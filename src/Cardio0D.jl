"""
    Cardio0D

Módulo do projeto PBL 2026.2 - Grupo 3:
"Modelagem 0D do sistema cardiovascular para simulação de choque hipovolêmico
e comparação de terapias de reposição volêmica em UTI".

Frentes integradas:
- Frente A: Modelagem matemática (Windkessel / SciML)
- Frente B: Fisiologia, viscosidade μ(Hct) e cenários terapêuticos
- Frente C: Calibração de parâmetros com dados MIMIC-IV
"""
module Cardio0D

using DifferentialEquations
using LinearAlgebra
using CSV, DataFrames

export BloodProperties, blood_viscosity, blood_viscosity_krieger, blood_viscosity_pries
export Windkessel2Params, windkessel2_ode!, simulate_windkessel2
export calculate_map
# parâmetros e fisiologia
export CircParams, baseline_params, with_params, lpm_to_mls, mls_to_lpm, cmH2O_to_mmHg
export Scenario, window_rate, activation, elastance, valve_flow
# modelos
export mean_rhs!, simulate_mean, initial_state_mean, mean_observables, mean_hemodynamics
export operating_point, oxygen_delivery
export pulsatile_rhs!, simulate_pulsatile, initial_state_pulsatile, pulsatile_observables
# análise
export time_mean, cycle_metrics, last_beats, venous_return_curve, cardiac_output_curve
export guyton_intersection, pv_loop
# dados e calibração
export PatientData, load_patient, build_patient, synthetic_patient, write_patient_csv
export read_mimic_table, scenario_from_patient, step_function, linear_function
export calibrate, simulate_patient, sensitivity_matrix, rmse, settled_state
export MIMIC_ITEMIDS

"""
    BloodProperties(hct, mu_plasma)

Propriedades reológicas do sangue.
- `hct`: hematócrito (fração volumétrica, ex: 0.40 para 40%)
- `mu_plasma`: viscosidade do plasma em cP (mPa·s), valor típico ~1.2 cP
"""
Base.@kwdef struct BloodProperties{T<:Real}
    hct::T = 0.40
    mu_plasma::T = 1.20
end

"""
    blood_viscosity(bp::BloodProperties)
    blood_viscosity(hct::Real; mu_plasma=1.20)

Calcula a viscosidade sanguínea efetiva μ (em cP) em função do hematócrito (Hct).

Relação adotada no modelo (entregável 1.5): μ = μ_plasma·exp(2.5·Hct), com coeficiente de
Einstein 2.5. Para Hct = 40%, dá ~3.26 cP (sangue total em vasos grandes).

O modelo só usa a **razão** `μ(Hct)/μ(Hct_ref)` (E18), e essa razão coincide com a de
[`blood_viscosity_pries`](@ref) em vaso grande (D = 1000 µm) dentro de 3% para Hct entre 0,15
e 0,55 (ver `docs/viscosidade.md`). Em microvasos (D = 30–100 µm) o efeito Fåhræus–Lindqvist
torna a resposta ao Hct mais fraca: use `blood_viscosity_pries(hct; D)` para a análise de
sensibilidade.
"""
function blood_viscosity(hct::Real; mu_plasma::Real=1.20)
    return mu_plasma * exp(2.5 * hct)
end

function blood_viscosity(bp::BloodProperties)
    return blood_viscosity(bp.hct; mu_plasma=bp.mu_plasma)
end

"""
    Windkessel2Params

Parâmetros do modelo arterial Windkessel de 2 elementos:
- `R`: Resistência vascular sistêmica (RVS) [mmHg·s/mL]
- `C`: Complacência arterial total [mL/mmHg]
- `Pv`: Pressão venosa central de referência [mmHg] (default 4.0 mmHg)
"""
Base.@kwdef struct Windkessel2Params{T<:Real}
    R::T = 1.0       # Resistência vascular sistêmica baseline
    C::T = 1.2       # Complacência arterial baseline
    Pv::T = 4.0      # Pressão venosa
end

"""
    windkessel2_ode!(dP, P, params, t)

Equação diferencial ordinária do Windkessel de 2 elementos:
dPa/dt = (Qin(t) - (Pa - Pv)/R) / C
"""
function windkessel2_ode!(dP, P, p, t)
    # p pode ser uma tupla (params::Windkessel2Params, Qin_func)
    params, Qin_func = p
    Pa = P[1]
    Qin = Qin_func(t)
    
    # Derivada dPa/dt
    dP[1] = (Qin - (Pa - params.Pv) / params.R) / params.C
end

"""
    calculate_map(sol)

Calcula a Pressão Arterial Média (PAM) a partir da série temporal de pressão simulada.
"""
function calculate_map(times, pressures)
    # Integração trapezoidal simples
    N = length(times)
    N < 2 && return pressures[1]
    dt_total = times[end] - times[1]
    dt_total <= 0 && return pressures[end]
    
    integral = 0.0
    for i in 1:(N-1)
        dt = times[i+1] - times[i]
        integral += 0.5 * (pressures[i] + pressures[i+1]) * dt
    end
    return integral / dt_total
end

"""
    blood_viscosity_krieger(hct; mu_plasma=1.20, hct_max=0.62, intrinsic=2.5)

Alternativa tipo Krieger–Dougherty: `μ = μ_pl·(1 − Hct/Hct_max)^(−[η]·Hct_max)`.
**Não recomendada com os valores padrão**: superestima a viscosidade (≈ 6 cP em Hct 0,40 e
≈ 35 cP em 0,55, contra ≈ 3,4 e 5,0 cP de Pries) porque sangue exige `[η]·Hct_max` bem
menor que 1,55. Mantida só para comparação (`docs/viscosidade.md`).
"""
function blood_viscosity_krieger(hct::Real; mu_plasma::Real=1.20, hct_max::Real=0.62,
                                 intrinsic::Real=2.5)
    return mu_plasma * (1 - hct / hct_max)^(-intrinsic * hct_max)
end

"""
    blood_viscosity_pries(hct; D=100.0, mu_plasma=1.20, invivo=true)

Viscosidade aparente do sangue de Pries et al. em cP, para um vaso de diâmetro `D` [µm] e
hematócrito de descarga `hct`. Forma comum:

```
μ_rel = 1 + (μ_0.45 − 1)·((1 − Hct)^C − 1)/((1 − 0,45)^C − 1)
C     = (0,8 + exp(−0,075·D))·(−1 + 1/(1 + 1e-11·D¹²)) + 1/(1 + 1e-11·D¹²)
```

- `invivo=false` — tubos de vidro (Pries, Neuhaus & Gaehtgens 1992, Am J Physiol 263:H1770):
  `μ_0.45 = 220·exp(−1,3·D) + 3,2 − 2,44·exp(−0,06·D^0,645)`.
- `invivo=true` — microvasos (Pries et al. 1994, Circ Res 75:904):
  `μ_0.45 = 6·exp(−0,085·D) + 3,2 − 2,44·exp(−0,06·D^0,645)` e
  `μ_rel,vivo = [1 + (μ_0.45 − 1)·(…)·(D/(D−1,1))²]·(D/(D−1,1))²` (camada de superfície do
  endotélio).

`μ = μ_pl·μ_rel`. Para D grande (≳ 300 µm) as duas versões coincidem. Fonte das fórmulas:
tabela B.2 da Frente B (valores conferidos em `docs/viscosidade.md`).
"""
function blood_viscosity_pries(hct::Real; D::Real=100.0, mu_plasma::Real=1.20, invivo::Bool=true)
    μ045 = (invivo ? 6 * exp(-0.085 * D) : 220 * exp(-1.3 * D)) + 3.2 - 2.44 * exp(-0.06 * D^0.645)
    g = 1 / (1 + 1e-11 * D^12)
    C = (0.8 + exp(-0.075 * D)) * (-1 + g) + g
    shape = ((1 - hct)^C - 1) / ((1 - 0.45)^C - 1)
    if invivo
        f = (D / (D - 1.1))^2
        μrel = (1 + (μ045 - 1) * shape * f) * f
    else
        μrel = 1 + (μ045 - 1) * shape
    end
    return mu_plasma * μrel
end

include("params.jl")
include("physiology.jl")
include("mean_model.jl")
include("closed_loop.jl")
include("analysis.jl")
include("data_io.jl")
include("calibration.jl")

end # module Cardio0D
