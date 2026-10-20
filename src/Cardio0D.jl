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

[PLACEHOLDER - Entregável B.2]:
Aproximação inicial para escoamento macroscópico baseada no coeficiente de Einstein (2.5)
para suspensões: μ ≈ μ_plasma * exp(2.5 * Hct).
Para Hct = 40%, fornece ~3.26 cP (faixa fisiológica razoável para sangue total em grandes vasos).

*Nota metodológica*: As equações empíricas completas de Pries et al. (1992) descrevem a
viscosidade aparente in vivo levando em consideração também o diâmetro do vaso (D) e o
efeito Fåhræus-Lindqvist. A formulação analítica detalhada com diâmetros e referências
será integrada na entrega B.2 pela Frente B.
"""
function blood_viscosity(hct::Real; mu_plasma::Real=1.20)
    # Placeholder analítico: relação exponencial aproximada (coeficiente de Einstein)
    # Sujeito a refinamento com a formulação completa da Frente B (Entregável B.2)
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
