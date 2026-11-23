# src/calibration.jl
# Calibração por paciente (CLAUDE.md §10): erro quadrático da PAM (e, opcionalmente, da FC) +
# regularização em log(θ), otimizador Nelder–Mead (poucos parâmetros) e sensibilidade local.
#
# Para calibrar um modelo que depois vai PREVER outros tratamentos, a FC não é imposta: o
# barorreflexo do modelo a gera, e a FC medida entra no erro (a mesma física nas duas etapas).

using Optimization, OptimizationOptimJL

"""
    settled_state(p, sc; T=1800.0)

Estado "compensado" no início da janela: parte do equilíbrio de volume `p.Vb0` e deixa o
barorreflexo, o interstício e o vasopressor assentarem por `T` s com as entradas do instante
inicial mantidas constantes. Evita um transiente artificial no começo da simulação.

Numa hipovolemia grave (perda ≳ 1 L) não existe equilíbrio **sem** os efetores, que é o ponto de
partida de `initial_state_mean`. Nesse caso parte-se da normovolemia (`Vb_ref`) e o déficit
`Vb_ref − Vb0` é retirado em 10 min, como uma hemorragia, com o barorreflexo atuando.
"""
function settled_state(p, sc::Scenario; T=1800.0)
    t0 = sc.tstops isa Vector && !isempty(sc.tstops) ? sc.tstops[1] : 0.0
    q0 = sc.Qinf(t0); d0 = sc.dose(t0); h0 = sc.hct_inf isa Function ? sc.hct_inf(t0) : sc.hct_inf
    hr0 = sc.hr_forced === nothing ? nothing : (_ -> sc.hr_forced(t0))
    u0, Td = try
        initial_state_mean(p), 0.0
    catch
        initial_state_mean(with_params(p; Vb0=p.Vb_ref)), 600.0
    end
    loss = Td > 0 ? (p.Vb_ref - p.Vb0) / Td : 0.0
    sc0 = Scenario(Qhem = t -> (t < t0 + Td ? loss : 0.0), Qinf = _ -> q0, hct_inf = h0,
                   dose = _ -> d0, hr_forced = hr0, baro = sc.baro, tstops = [t0 + Td])
    sol = simulate_mean(p, sc0, (t0, t0 + Td + T); u0, reltol=1e-6, abstol=1e-6)
    return sol.u[end]
end

"""
    simulate_patient_series(p, pd; sc=scenario_from_patient(pd))

Séries simuladas nos instantes de `pd.t` (modelo lento, estado inicial compensado): PAM, FC,
CO [L/min], Hct e DO₂ [mL/min].
"""
function simulate_patient_series(p, pd::PatientData; sc=scenario_from_patient(pd))
    u0 = settled_state(p, sc)
    sol = simulate_mean(p, sc, (pd.t[1], pd.t[end]); u0, saveat=pd.t, reltol=1e-6, abstol=1e-6)
    o = mean_observables(sol, p, sc, pd.t)
    return (; PAM=o.PAM, HR=o.HR, CO=o.CO, Hct=o.Hct, DO2=o.DO2)
end

"PAM simulada nos instantes de `pd.t` (ver [`simulate_patient_series`](@ref))."
simulate_patient(p, pd::PatientData; sc=scenario_from_patient(pd)) =
    simulate_patient_series(p, pd; sc).PAM

"Raiz do erro quadrático médio ignorando NaN (NaN se não houver nenhum par válido)."
function rmse(a, b)
    ok = findall(i -> !isnan(a[i]) && !isnan(b[i]), eachindex(a))
    isempty(ok) && return NaN
    return sqrt(sum(abs2, a[ok] .- b[ok]) / length(ok))
end

_theta_params(p, names, θ) = with_params(p; (n => v for (n, v) in zip(names, θ))...)

"""
    calibrate(p0, pd, names; λ=0.0, w_hr=0.0, maxiters=600, sc=scenario_from_patient(pd))

Ajusta os parâmetros `names` (ex.: `[:Rsa, :Vb0, :HR]`) minimizando

    RMSE(PAM)² + w_hr·RMSE(FC)² + λ·Σ(log θ − log θ₀)²

em `log θ` (garante θ > 0 e escalas comparáveis). `w_hr = 0,25` faz um erro de 10 bpm pesar
como 5 mmHg; com `w_hr = 0` só a PAM entra. `p0` fornece o valor inicial e o *prior*.

`limites = [(lo, hi), ...]` restringe cada parâmetro à faixa `[lo·θ₀, hi·θ₀]` (faixa fisiológica):
`log(θ/θ₀) = m + h·tanh(x)`, com `m` e `h` o centro e a meia-largura da faixa em escala log, de modo
que o otimizador continua sem restrições em `x`. Devolve `(p, θ, rmse, rmse_hr, resultado)`.
"""
function calibrate(p0, pd::PatientData, names; λ=0.0, w_hr=0.0, maxiters=600, limites=nothing,
                   sc=scenario_from_patient(pd))
    θ0 = [float(getfield(p0, n)) for n in names]
    if limites === nothing
        logr = x -> x
        x0 = zeros(length(names))
    else
        m = [log(lo * hi) / 2 for (lo, hi) in limites]
        h = [log(hi / lo) / 2 for (lo, hi) in limites]
        logr = x -> m .+ h .* tanh.(x)
        x0 = atanh.(clamp.(-m ./ h, -0.99, 0.99))         # começa em θ₀ (ou na borda mais próxima)
    end
    function loss(x, _)
        p = _theta_params(p0, names, θ0 .* exp.(logr(x)))
        y = try
            simulate_patient_series(p, pd; sc)
        catch
            return 1e6
        end
        e_hr = w_hr > 0 ? rmse(y.HR, pd.HR) : 0.0
        return rmse(y.PAM, pd.PAM)^2 + (isnan(e_hr) ? 0.0 : w_hr * e_hr^2) + λ * sum(abs2, logr(x))
    end
    prob = OptimizationProblem(OptimizationFunction(loss), x0)
    res = solve(prob, NelderMead(); maxiters)
    θ = θ0 .* exp.(logr(res.u))
    p = _theta_params(p0, names, θ)
    y = simulate_patient_series(p, pd; sc)
    return (; p, θ, rmse=rmse(y.PAM, pd.PAM), rmse_hr=rmse(y.HR, pd.HR), result=res)
end

"""
    sensitivity_matrix(p, pd, names; h=0.01)

Sensibilidade local normalizada `S_ij = (∂PAM_i/∂θ_j)·(θ_j/PAM_i)` por diferenças centrais
(§10.3). A condição de `SᵀS` mede a identificabilidade prática: colunas quase colineares
(condição enorme) ou uma sensibilidade média muito menor que a dos demais indicam parâmetros
que a PAM sozinha não identifica (é o caso de `EmaxL`: ver `docs/resultados_mimic.md`).
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
