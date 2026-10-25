# src/physiology.jl
# Peças elementares: elastância, válvulas, pressão venosa, sangue, balanço de fluidos,
# vasopressor e barorreflexo. Tudo que é compartilhado pelos modelos pulsátil e "lento".

# ---------------------------------------------------------------- funções elementares
"Sigmoide logística numericamente estável (sem overflow de exp)."
stable_sigmoid(x) = x >= 0 ? inv(1 + exp(-x)) : (e = exp(x); e / (1 + e))

"""
    activation(τ, T)

Ativação normalizada `e(τ) ∈ [0,1]` do ventrículo (CLAUDE.md §4.4): contração e relaxamento em
cosseno; a sístole encurta pouco com a taquicardia (`Ts1 = 0,3·√T`).
"""
function activation(τ, T)
    Ts1 = 0.3 * sqrt(T)
    Ts2 = 0.5 * Ts1
    if τ < Ts1
        return 0.5 * (1 - cos(π * τ / Ts1))
    elseif τ < Ts1 + Ts2
        return 0.5 * (1 + cos(π * (τ - Ts1) / Ts2))
    else
        return zero(τ)
    end
end

"Elastância variante no tempo `E(t) = Emin + (Emax − Emin)·e` (equação E16)."
elastance(Emin, Emax, e) = Emin + (Emax - Emin) * e

"""
    valve_flow(ΔP, R, k)

Válvula como diodo suave: `Q = (ΔP/R)·σ(k·ΔP)` (CLAUDE.md §4.5). Para `ΔP ≫ 1/k` é a lei de
Ohm; para `ΔP < 0` o fluxo é praticamente nulo, e a função é diferenciável.
"""
valve_flow(ΔP, R, k) = ΔP / R * stable_sigmoid(k * ΔP)

"""
    venous_pressure(stressed, C, s)

Pressão venosa a partir do volume estressado `V − V0`: linear (`stressed/C`) quando positiva e
saturando em 0 quando a veia colapsa (softplus de largura `s`). Sem isso uma hemorragia faria
a pressão venosa ficar fortemente negativa.
"""
function venous_pressure(stressed, C, s)
    x = stressed / C
    return x > 20 * s ? x : s * log1p(exp(x / s))
end

"Inverso de [`venous_pressure`](@ref): volume estressado que gera a pressão `P > 0`."
function stressed_volume_from_pressure(P, C, s)
    return P > 20 * s ? C * P : C * s * log(expm1(P / s))
end

# ------------------------------------------------------------------ estados "extras"
# Compartilhados pelos dois modelos (vêm sempre no fim do vetor de estados):
#   1 Vrbc  volume de hemácias [mL]           5 xR   efeito simpático na resistência
#   2 Vis   excesso de líquido intersticial   6 xV   efeito simpático em V0 venoso
#   3 Pbar  PAM filtrada pelo barorreceptor   7 xE   efeito simpático na contratilidade
#   4 xHR   efeito simpático na FC [bpm]      8 Ce   concentração no sítio de efeito
const N_EXTRA = 8

"""
    Scenario

Entradas externas de uma simulação (funções do tempo em mL/s ou unidade da droga):
- `Qhem(t)`: hemorragia; `Qinf(t)`: infusão; `hct_inf`: hematócrito do fluido infundido
  (0 = cristaloide, ~0,40 = sangue total; número ou função do tempo);
- `dose(t)`: dose de vasopressor (mesma unidade de `EC50`);
- `baro`: liga/desliga o barorreflexo; `hr_forced(t)`: impõe a FC (ex.: dado do MIMIC);
- `tstops`: instantes de descontinuidade (início/fim de infusões) para o solver.
"""
Base.@kwdef struct Scenario
    Qhem::Function = t -> 0.0
    Qinf::Function = t -> 0.0
    hct_inf::Union{Float64,Function} = 0.0
    dose::Function = t -> 0.0
    baro::Bool = true
    hr_forced::Union{Nothing,Function} = nothing
    tstops::Vector{Float64} = Float64[]
end

"Taxa constante `volume/dur` (mL/s) entre `t0` e `t0+dur`."
window_rate(volume, t0, dur) = t -> (t0 <= t < t0 + dur) ? volume / dur : 0.0

"""
    effectors(p, sc, t, ex, Vblood)

Traduz os estados extras em parâmetros efetivos do instante: hematócrito, fator de
viscosidade (`μ(Hct)/μ(Hct_ref)`, E18), FC, `Rsa`, `Rpa`, `V0sv`, `Emax` e o sinal do
barorreflexo `σ ∈ (−1,1)` (positivo quando a pressão está baixa).
"""
function effectors(p, sc, t, ex, Vblood)
    Vrbc, Vis, Pbar, xHR, xR, xV, xE, Ce = ex
    Hct = Vrbc / Vblood
    μr = blood_viscosity(Hct; mu_plasma=p.mu_plasma) /
         blood_viscosity(p.Hct_ref; mu_plasma=p.mu_plasma)
    f = Ce / (p.EC50 + Ce)                       # Emax, n = 1
    HR = sc.hr_forced === nothing ? clamp(p.HR + xHR, 30.0, 220.0) : sc.hr_forced(t)
    σ = sc.baro ? 1 - 2 * stable_sigmoid((Pbar - p.Pn) / p.ks) : zero(Pbar)
    return (; Hct, μr, f, σ, HR,
            Rsa=p.Rsa * μr * (1 + xR) * (1 + p.G_Rd * f),
            Rpa=p.Rpa * μr,
            V0sv=p.V0sv * (1 - xV) * (1 - p.G_Vd * f),
            EmaxL=p.EmaxL * (1 + xE) * (1 + p.G_Ed * f),
            EmaxR=p.EmaxR * (1 + xE) * (1 + p.G_Ed * f))
end

"""
    capillary_exchange(p, Vblood, Vrbc, Vis)

Fluxo líquido plasma → interstício `J` [mL/s] (§7.2, nível 1): só o excesso de **plasma**
(não as hemácias) é filtrado; o refil, no sentido inverso, satura em `Jmax`.
"""
function capillary_exchange(p, Vblood, Vrbc, Vis)
    Vpl_ref = p.Vb_ref * (1 - p.Hct_ref)
    x = ((Vblood - Vrbc) - Vpl_ref) - Vis / p.rho
    return x > 0 ? p.kf * x : p.Jmax * tanh(p.kf * x / p.Jmax)
end

"""
    extras_rhs!(du, u, o, p, sc, t, e, Psa, Vblood, J)

Equações dos 8 estados extras (hemácias, intersticial, barorreflexo, vasopressor), com
`o` o deslocamento do primeiro estado extra em `u`.
"""
function extras_rhs!(du, u, o, p, sc, t, e, Psa, Vblood, J)
    Qinf = sc.Qinf(t); Qhem = sc.Qhem(t)
    hinf = sc.hct_inf isa Function ? sc.hct_inf(t) : sc.hct_inf
    du[o+1] = hinf * Qinf - e.Hct * Qhem
    du[o+2] = J
    du[o+3] = (Psa - u[o+3]) / p.tau_p
    du[o+4] = (p.G_HR * e.σ - u[o+4]) / p.tau_HR
    du[o+5] = (p.G_R * e.σ - u[o+5]) / p.tau_R
    du[o+6] = (p.G_V * e.σ - u[o+6]) / p.tau_V
    du[o+7] = (p.G_E * e.σ - u[o+7]) / p.tau_E
    du[o+8] = (sc.dose(t) - u[o+8]) / p.tau_e
    return nothing
end

"Estados extras no basal: hemácias para `Hct_ref`, sem líquido extra, reflexo em repouso."
extras_baseline(p, Vb, Pbar) = [p.Hct_ref * Vb, 0.0, Pbar, 0.0, 0.0, 0.0, 0.0, 0.0]
