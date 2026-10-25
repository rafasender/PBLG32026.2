# src/calibration.jl
# Calibração por paciente (CLAUDE.md §10): erro quadrático da PAM + regularização em log(θ),
# otimizador Nelder–Mead (poucos parâmetros) e sensibilidade local normalizada.

using Optimization, OptimizationOptimJL

"""
    settled_state(p, sc; T=1800.0)

Estado "compensado" no início da janela: parte do equilíbrio de volume `p.Vb0` e deixa o
barorreflexo, o interstício e o vasopressor assentarem por `T` s com as entradas do instante
inicial mantidas constantes. Evita um transiente artificial no começo da simulação.
"""
function settled_state(p, sc::Scenario; T=1800.0)
    t0 = sc.tstops isa Vector && !isempty(sc.tstops) ? sc.tstops[1] : 0.0
    q0 = sc.Qinf(t0); d0 = sc.dose(t0); h0 = sc.hct_inf isa Function ? sc.hct_inf(t0) : sc.hct_inf
    hr0 = sc.hr_forced === nothing ? nothing : (t -> sc.hr_forced(t0))
    sc0 = Scenario(Qinf = t -> q0, hct_inf = h0, dose = t -> d0, hr_forced = hr0, baro = sc.baro)
    sol = simulate_mean(p, sc0, (t0, t0 + T); reltol=1e-6, abstol=1e-6)
    return sol.u[end]
end

"""
    simulate_patient(p, pd; sc=scenario_from_patient(pd))

PAM simulada nos instantes de `pd.t` (modelo lento, estado inicial compensado).
"""
function simulate_patient(p, pd::PatientData; sc=scenario_from_patient(pd))
    u0 = settled_state(p, sc)
    sol = simulate_mean(p, sc, (pd.t[1], pd.t[end]); u0, saveat=pd.t, reltol=1e-6, abstol=1e-6)
    return [mean_hemodynamics(p, sol(t), sc, t).Psa for t in pd.t]
end

"Raiz do erro quadrático médio ignorando NaN."
function rmse(a, b)
    ok = findall(i -> !isnan(a[i]) && !isnan(b[i]), eachindex(a))
    return sqrt(sum(abs2, a[ok] .- b[ok]) / length(ok))
end

_theta_params(p, names, θ) = with_params(p; (n => v for (n, v) in zip(names, θ))...)

"""
    calibrate(p0, pd, names; λ=0.0, maxiters=600, sc=scenario_from_patient(pd))

Ajusta os parâmetros `names` (ex.: `[:Rsa, :Vb0, :EmaxL]`) minimizando
`MSE(PAM_sim, PAM_obs) + λ·Σ(log θ − log θ₀)²` em `log θ` (garante θ > 0 e escalas
comparáveis). `p0` fornece o valor inicial e o *prior*. Devolve `(p, θ, rmse, resultado)`.
"""
function calibrate(p0, pd::PatientData, names; λ=0.0, maxiters=600, sc=scenario_from_patient(pd))
    θ0 = [float(getfield(p0, n)) for n in names]
    function loss(x, _)
        θ = θ0 .* exp.(x)
        p = _theta_params(p0, names, θ)
        y = try
            simulate_patient(p, pd; sc)
        catch
            return 1e6
        end
        return rmse(y, pd.PAM)^2 + λ * sum(abs2, x)
    end
    prob = OptimizationProblem(OptimizationFunction(loss), zeros(length(names)))
    res = solve(prob, NelderMead(); maxiters)
    θ = θ0 .* exp.(res.u)
    p = _theta_params(p0, names, θ)
    return (; p, θ, rmse=rmse(simulate_patient(p, pd; sc), pd.PAM), result=res)
end

"""
    sensitivity_matrix(p, pd, names; h=0.01)

Sensibilidade local normalizada `S_ij = (∂PAM_i/∂θ_j)·(θ_j/PAM_i)` por diferenças centrais
(§10.3). A condição de `SᵀS` mede a identificabilidade prática: colunas quase colineares
(condição enorme) ou uma sensibilidade média muito menor que a dos demais indicam parâmetros
que a PAM sozinha não identifica (no ensaio sintético, `EmaxL`: |S| ≈ 0,015 contra 0,37 de `Vb0`).
Devolve `(S, cond_SᵀS)`.
"""
function sensitivity_matrix(p, pd::PatientData, names; h=0.01, sc=scenario_from_patient(pd))
    y0 = simulate_patient(p, pd; sc)
    S = zeros(length(pd.t), length(names))
    for (j, n) in enumerate(names)
        θ = float(getfield(p, n))
        yp = simulate_patient(with_params(p; n => θ * (1 + h)), pd; sc)
        ym = simulate_patient(with_params(p; n => θ * (1 - h)), pd; sc)
        S[:, j] = (yp .- ym) ./ (2 * h * θ) .* θ ./ y0
    end
    return S, cond(S' * S)
end
