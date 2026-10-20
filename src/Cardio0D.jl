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

export BloodProperties, blood_viscosity
export Windkessel2Params, windkessel2_ode!, simulate_windkessel2
export calculate_map

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
Baseado nas formulações reológicas da literatura hemodinâmica (Pries et al., 1992):
μ(Hct) aumenta exponencialmente/não linearmente com o hematócrito.
Para a faixa fisiológica padrão:
μ ≈ μ_plasma * exp(2.5 * Hct)
"""
function blood_viscosity(hct::Real; mu_plasma::Real=1.20)
    # Modelo exponencial validado para hemodiluição e policitemia moderada
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

end # module Cardio0D
