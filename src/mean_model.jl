# src/mean_model.jl
# Modelo "lento" (valores médios por batimento) — CLAUDE.md §5.5 — e ponto de funcionamento.
#
# Estados: u = [Vsa, Vsv, Vpa, Vpv | Vrbc, Vis, Pbar, xHR, xR, xV, xE, Ce]
# O coração vira a fonte algébrica da aula: CO = HR · (C_D·P_f − C_S·P_pós) (equação E14).

const MEAN_N = 4          # estados hemodinâmicos
const MEAN_O = 4          # deslocamento dos estados extras

"""
    mean_hemodynamics(p, u, sc, t)

Pressões, fluxos e efetores instantâneos do modelo lento. É a única fonte das relações
constitutivas: o lado direito da EDO e a análise posterior usam esta mesma função.
"""
function mean_hemodynamics(p, u, sc, t)
    Vsa, Vsv, Vpa, Vpv = u[1], u[2], u[3], u[4]
    ex = (u[5], u[6], u[7], u[8], u[9], u[10], u[11], u[12])
    Vblood = Vsa + Vsv + Vpa + Vpv + p.Vheart
    e = effectors(p, sc, t, ex, Vblood)

    Psa = (Vsa - p.V0sa) / p.Csa
    Psv = venous_pressure(Vsv - e.V0sv, p.Csv, p.sveno)
    Ppa = (Vpa - p.V0pa) / p.Cpa + p.Pth
    Ppv = (Vpv - p.V0pv) / p.Cpv + p.Pth

    k = e.HR / 60
    CDL = 1 / p.EminL; CSL = 1 / e.EmaxL
    CDR = 1 / p.EminR; CSR = 1 / e.EmaxR
    afterL = Psa + p.dPesL - p.Pth
    afterR = Ppa + p.dPesR - p.Pth

    # coração esquerdo: SV = C_D·(P_f − P_th) − C_S·(P_pós − P_th), limitado entre 0 e EDVmax
    SVL = max(0.0, min(CDL * (Ppv - p.Pth), p.EDVmaxL - p.VdL) - CSL * afterL)
    # coração direito: a resistência venosa Rsv fica em série com o enchimento (RAP = Psv − Q·Rsv)
    SVR_free = (CDR * (Psv - p.Pth) - CSR * afterR) / (1 + k * CDR * p.Rsv)
    SVR = max(0.0, min(SVR_free, p.EDVmaxR - p.VdR - CSR * afterR))
    QL = k * SVL
    QR = k * SVR

    Qsys = (Psa - Psv) / e.Rsa
    Qlung = (Ppa - Ppv) / e.Rpa
    J = capillary_exchange(p, Vblood, ex[1], ex[2])
    RAP = Psv - QR * p.Rsv
    return (; Vblood, e, ex, Psa, Psv, Ppa, Ppv, RAP, QL, QR, SVL, SVR, Qsys, Qlung, J)
end

function mean_rhs!(du, u, par, t)
    p, sc = par
    h = mean_hemodynamics(p, u, sc, t)
    du[1] = h.QL - h.Qsys
    du[2] = h.Qsys - h.QR + sc.Qinf(t) - sc.Qhem(t) - h.J
    du[3] = h.QR - h.Qlung
    du[4] = h.Qlung - h.QL
    extras_rhs!(du, u, MEAN_O, p, sc, t, h.e, h.Psa, h.Vblood, h.J)
    return nothing
end

# ----------------------------------------------------- ponto de funcionamento algébrico
"Raiz de `f` em `[a,b]` por bisseção (exige troca de sinal)."
function bisect(f, a, b; tol=1e-10, maxit=200)
    fa, fb = f(a), f(b)
    fa == 0 && return a
    fb == 0 && return b
    fa * fb > 0 && error("bisect: sem troca de sinal em [$a, $b]")
    for _ in 1:maxit
        m = (a + b) / 2
        fm = f(m)
        (fm == 0 || (b - a) / 2 < tol) && return m
        if fa * fm < 0
            b, fb = m, fm
        else
            a, fa = m, fm
        end
    end
    return (a + b) / 2
end

"""
    pressures_for_flow(p, Q)

Dado o débito `Q` [mL/s], resolve o sistema linear 4×4 (Ohm nas duas circulações + as duas
equações do coração, §4.2 e §5.4) e devolve `(Psa, Psv, Ppa, Ppv)`; basal (sem efetores).
"""
function pressures_for_flow(p, Q)
    k = p.HR / 60
    CDL = 1 / p.EminL; CSL = 1 / p.EmaxL
    CDR = 1 / p.EminR; CSR = 1 / p.EmaxR
    # incógnitas: Psa, Psv, Ppa, Ppv
    A = [1.0      -1.0     0.0      0.0;
         0.0       0.0     1.0     -1.0;
         -k*CSL    0.0     0.0      k*CDL;
         0.0       k*CDR  -k*CSR    0.0]
    b = [Q * p.Rsa,
         Q * p.Rpa,
         Q + k * (CDL * p.Pth + CSL * (p.dPesL - p.Pth)),
         Q + k * (CDR * (Q * p.Rsv + p.Pth) + CSR * (p.dPesR - p.Pth))]
    x = A \ b
    return x[1], x[2], x[3], x[4]
end

"""
    volumes_from_pressures(p, Psa, Psv, Ppa, Ppv)

Inverte as relações P–V dos quatro compartimentos (parâmetros basais).
"""
function volumes_from_pressures(p, Psa, Psv, Ppa, Ppv)
    Vsa = p.V0sa + p.Csa * Psa
    Vsv = p.V0sv + stressed_volume_from_pressure(Psv, p.Csv, p.sveno)
    Vpa = p.V0pa + p.Cpa * (Ppa - p.Pth)
    Vpv = p.V0pv + p.Cpv * (Ppv - p.Pth)
    return Vsa, Vsv, Vpa, Vpv
end

"""
    operating_point(p; Vb=p.Vb0)

Equilíbrio (CO = VR) do modelo lento para um volume sanguíneo `Vb`, sem efetores e sem trocas
com o interstício. Varre o débito `Q` e resolve `V_total(Q) = Vb` por bisseção; é a
verificação analítica contra a qual a EDO é testada (CLAUDE.md §9, V5).
Devolve uma `NamedTuple` com volumes, pressões, `Q`, `CO` [L/min], `RAP` e `Pms`.
"""
function operating_point(p; Vb=p.Vb0)
    function Vtotal(Q)
        Psa, Psv, Ppa, Ppv = pressures_for_flow(p, Q)
        (Psv <= 0.02 || Ppv <= p.Pth) && return NaN
        return sum(volumes_from_pressures(p, Psa, Psv, Ppa, Ppv)) + p.Vheart
    end
    g(Q) = Vtotal(Q) - Vb
    Qs = range(1.0, 400.0; length=800)
    vals = [g(Q) for Q in Qs]
    idx = findfirst(i -> !isnan(vals[i]) && !isnan(vals[i+1]) && vals[i] * vals[i+1] <= 0,
                    1:length(Qs)-1)
    idx === nothing && error("operating_point: sem equilíbrio para Vb = $Vb mL (colapso hemodinâmico)")
    Q = bisect(g, Qs[idx], Qs[idx+1])
    Psa, Psv, Ppa, Ppv = pressures_for_flow(p, Q)
    Vsa, Vsv, Vpa, Vpv = volumes_from_pressures(p, Psa, Psv, Ppa, Ppv)
    Pms = ((Vsa + Vsv) - (p.V0sa + p.V0sv)) / (p.Csa + p.Csv)
    return (; Q, CO=mls_to_lpm(Q), Psa, Psv, Ppa, Ppv, RAP=Psv - Q * p.Rsv, Pms,
            Vsa, Vsv, Vpa, Vpv, Vb)
end

"Estado inicial do modelo lento: equilíbrio para `p.Vb0` + extras no basal."
function initial_state_mean(p)
    op = operating_point(p; Vb=p.Vb0)
    return vcat([op.Vsa, op.Vsv, op.Vpa, op.Vpv], extras_baseline(p, p.Vb0, op.Psa))
end

"""
    simulate_mean(p, sc, tspan; u0=initial_state_mean(p), solver=Tsit5(), kwargs...)

Resolve o modelo lento. `kwargs` vão para `solve` (`saveat`, `reltol`, ...).
"""
function simulate_mean(p, sc::Scenario, tspan; u0=initial_state_mean(p), solver=Tsit5(),
                       reltol=1e-8, abstol=1e-8, kwargs...)
    prob = ODEProblem(mean_rhs!, u0, tspan, (p, sc))
    return solve(prob, solver; reltol, abstol, tstops=sc.tstops, kwargs...)
end

"""
    mean_observables(sol, p, sc, ts)

Séries derivadas nos instantes `ts`: pressões, CO [L/min], FC, Hct, viscosidade relativa,
volume sanguíneo, líquido intersticial e DO₂ [mL/min].
"""
function mean_observables(sol, p, sc, ts)
    hs = [mean_hemodynamics(p, sol(t), sc, t) for t in ts]
    return observables_from(hs, p, ts)
end

function observables_from(hs, p, ts)
    CO = [mls_to_lpm(h.QL) for h in hs]
    Hct = [h.e.Hct for h in hs]
    return (; t=collect(ts),
            PAM=[h.Psa for h in hs], CVP=[h.Psv for h in hs], RAP=[h.RAP for h in hs],
            Ppa=[h.Ppa for h in hs], LAP=[h.Ppv for h in hs],
            CO, HR=[h.e.HR for h in hs], SV=[h.SVL for h in hs],
            Hct, mu_rel=[h.e.μr for h in hs], Vblood=[h.Vblood for h in hs],
            Vis=[h.ex[2] for h in hs], Pbar=[h.ex[3] for h in hs],
            DO2=[oxygen_delivery(c, hc) for (c, hc) in zip(CO, Hct)])
end

"""
    oxygen_delivery(CO_lpm, Hct; SaO2=0.98)

Oferta de O₂ [mL/min]: `DO₂ = CO·10·1,34·Hb·SaO₂`, com `Hb ≈ 33,3·Hct` g/dL (§11).
"""
oxygen_delivery(CO_lpm, Hct; SaO2=0.98) = CO_lpm * 10 * 1.34 * (33.3 * Hct) * SaO2
