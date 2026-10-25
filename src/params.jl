# src/params.jl
# Parâmetros da malha fechada, conversões de unidade e projeto do estado basal.
#
# Unidades: pressão [mmHg], volume [mL], tempo [s], fluxo [mL/s],
#           resistência [mmHg·s/mL], complacência [mL/mmHg], elastância [mmHg/mL].

lpm_to_mls(q) = q * 1000 / 60      # L/min -> mL/s
mls_to_lpm(q) = q * 60 / 1000      # mL/s  -> L/min
cmH2O_to_mmHg(x) = 0.7356 * x

"""
    CircParams

Parâmetros do modelo de circulação em malha fechada (CLAUDE.md §5, §7, §8, §12).
Crie sempre via [`baseline_params`](@ref) e altere com [`with_params`](@ref).

Grupos de campos:
- vasos sistêmicos e pulmonares (`R*`, `C*`, `V0*`), coração (`E*`, `Vd*`, `EDVmax*`);
- `dPesL/dPesR`: correção entre a pressão média e a pressão de fim de sístole, usada só no
  modelo "lento" (a pressão de pulso eleva a pós-carga efetiva do ventrículo);
- sangue (`Hct_ref`, `mu_plasma`), troca capilar (`kf`, `rho`, `Jmax`), vasopressor
  (`tau_e`, `EC50`, `G_*d`) e barorreflexo (`Pn`, `ks`, `tau_*`, `G_*`).
"""
Base.@kwdef struct CircParams{T<:Real}
    HR::T = 75.0
    # circulação sistêmica
    Rsa::T = 1.0;   Csa::T = 1.5;    V0sa::T = 600.0
    Rsv::T = 0.03;  Csv::T = 100.0;  V0sv::T = 2500.0
    # circulação pulmonar
    Rpa::T = 0.08;  Cpa::T = 4.0;    V0pa::T = 150.0
    Cpv::T = 10.0;  V0pv::T = 350.0
    # coração esquerdo / direito
    EmaxL::T = 2.5; EminL::T = 0.1;  VdL::T = 10.0
    EmaxR::T = 0.5; EminR::T = 0.05; VdR::T = 10.0
    EDVmaxL::T = 200.0; EDVmaxR::T = 210.0
    dPesL::T = 7.0; dPesR::T = 4.0
    # válvulas, tórax e volumes
    Rvalv::T = 0.005; kvalv::T = 10.0; Pth::T = -4.0
    sveno::T = 0.25                 # suavidade do colapso venoso [mmHg]
    Vheart::T = 175.0               # volume médio dos ventrículos (só no modelo lento)
    Vb_ref::T = 4800.0              # volume sanguíneo de referência (normovolemia)
    Vb0::T = 4800.0                 # volume sanguíneo no início da simulação
    # sangue e líquido intersticial
    Hct_ref::T = 0.42; mu_plasma::T = 1.2
    kf::T = 1 / (25 * 60)           # 1/s  (equilíbrio capilar ~25 min)
    rho::T = 3.0                    # V_intersticial / V_plasma no equilíbrio
    Jmax::T = 0.25                  # refil capilar máximo [mL/s]
    # vasopressor
    tau_e::T = 60.0; EC50::T = 0.05
    G_Rd::T = 1.0; G_Vd::T = 0.20; G_Ed::T = 0.25
    # barorreflexo
    Pn::T = 93.0; ks::T = 8.0; tau_p::T = 3.0
    G_HR::T = 45.0; G_R::T = 0.6; G_V::T = 0.25; G_E::T = 0.3
    tau_HR::T = 2.0; tau_R::T = 8.0; tau_V::T = 25.0; tau_E::T = 10.0
end

"""
    with_params(p; kwargs...)

Devolve uma cópia de `p` com os campos indicados trocados (os tipos são promovidos, então
funciona com números `Float64` ou duais).
"""
function with_params(p::CircParams; kw...)
    names = fieldnames(CircParams)
    args = Any[haskey(kw, f) ? kw[f] : getfield(p, f) for f in names]
    T = promote_type(typeof.(args)...)
    return CircParams{T}(convert.(T, args)...)
end

"""
    baseline_params(; kwargs...)

Projeta o estado basal a partir de **alvos hemodinâmicos** (adulto ~70 kg):
resistências pela lei de Ohm, e complacências/elastâncias do coração pelas equações
`EDV − Vd = C_D·(P_f − P_th)` e `ESV − Vd = C_S·(P_pós − P_th)` (CLAUDE.md §4).
O resultado reproduz exatamente os alvos no modelo "lento" (equilíbrio algébrico).
"""
function baseline_params(; HR=75.0, MAP=93.0, CVP=5.0, RAP=2.5, Ppa=15.0, LAP=8.0, CO=5.25,
                         Pth=-4.0, EDV_L=120.0, ESV_L=50.0, EDV_R=125.0, ESV_R=55.0, Vd=10.0,
                         dPesL=7.0, dPesR=4.0, Csa=1.5, Csv=100.0, Cpa=4.0, Cpv=10.0,
                         V0sa=600.0, V0sv=2500.0, V0pa=150.0, V0pv=350.0,
                         Rvalv=0.005, kvalv=10.0, overrides...)
    Q = lpm_to_mls(CO)
    Rsa = (MAP - CVP) / Q
    Rsv = (CVP - RAP) / Q
    Rpa = (Ppa - LAP) / Q
    CD_L = (EDV_L - Vd) / (LAP - Pth)
    CS_L = (ESV_L - Vd) / (MAP + dPesL - Pth)
    CD_R = (EDV_R - Vd) / (RAP - Pth)
    CS_R = (ESV_R - Vd) / (Ppa + dPesR - Pth)
    Vheart = (EDV_L + ESV_L) / 2 + (EDV_R + ESV_R) / 2
    Vves = (V0sa + Csa * MAP) + (V0sv + Csv * CVP) +
           (V0pa + Cpa * (Ppa - Pth)) + (V0pv + Cpv * (LAP - Pth))
    Vb = Vves + Vheart
    base = (; HR=float(HR), Rsa, Csa=float(Csa), V0sa=float(V0sa), Rsv, Csv=float(Csv),
            V0sv=float(V0sv), Rpa, Cpa=float(Cpa), V0pa=float(V0pa), Cpv=float(Cpv),
            V0pv=float(V0pv), EmaxL=1 / CS_L, EminL=1 / CD_L, VdL=float(Vd),
            EmaxR=1 / CS_R, EminR=1 / CD_R, VdR=float(Vd),
            EDVmaxL=EDV_L + 80.0, EDVmaxR=EDV_R + 85.0, dPesL=float(dPesL),
            dPesR=float(dPesR), Rvalv=float(Rvalv), kvalv=float(kvalv), Pth=float(Pth),
            Vheart, Vb_ref=Vb, Vb0=Vb, Pn=float(MAP))
    return CircParams(; merge(base, NamedTuple(overrides))...)
end
