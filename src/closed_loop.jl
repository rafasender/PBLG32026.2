# src/closed_loop.jl
# Modelo pulsátil em malha fechada (CLAUDE.md §5.1–5.3): 6 volumes + fase cardíaca + extras.
#
# Estados: u = [VL, Vsa, Vsv, VR, Vpa, Vpv, φ | Vrbc, Vis, Pbar, xHR, xR, xV, xE, Ce]
#   φ = fase do ciclo (em batimentos); dφ/dt = HR/60 permite FC variável sem descontinuidade.

const PULS_N = 7
const PULS_O = 7

"""
    pulsatile_hemodynamics(p, u, sc, t)

Relações constitutivas e fluxos instantâneos da malha fechada (usada pelo lado direito da
EDO e pela análise posterior).
"""
function pulsatile_hemodynamics(p, u, sc, t)
    VL, Vsa, Vsv, VR, Vpa, Vpv, φ = u[1], u[2], u[3], u[4], u[5], u[6], u[7]
    ex = (u[8], u[9], u[10], u[11], u[12], u[13], u[14], u[15])
    Vblood = VL + Vsa + Vsv + VR + Vpa + Vpv
    e = effectors(p, sc, t, ex, Vblood)

    frac = φ - floor(φ)                 # `mod` daria NaN na derivada (AD) em φ inteiro
    act = activation(frac * 60 / e.HR, 60 / e.HR)
    EL = elastance(p.EminL, e.EmaxL, act)
    ER = elastance(p.EminR, e.EmaxR, act)

    PL = EL * (VL - p.VdL) + p.Pth
    PR = ER * (VR - p.VdR) + p.Pth
    Psa = (Vsa - p.V0sa) / p.Csa
    Psv = venous_pressure(Vsv - e.V0sv, p.Csv, p.sveno)
    Ppa = (Vpa - p.V0pa) / p.Cpa + p.Pth
    Ppv = (Vpv - p.V0pv) / p.Cpv + p.Pth

    Qao = valve_flow(PL - Psa, p.Rvalv, p.kvalv)
    Qmi = valve_flow(Ppv - PL, p.Rvalv, p.kvalv)
    Qpul = valve_flow(PR - Ppa, p.Rvalv, p.kvalv)
    Qtri = valve_flow(Psv - PR, p.Rvalv + p.Rsv, p.kvalv)   # Rsv em série com a tricúspide
    Qsys = (Psa - Psv) / e.Rsa
    Qlung = (Ppa - Ppv) / e.Rpa
    J = capillary_exchange(p, Vblood, ex[1], ex[2])
    return (; Vblood, e, ex, act, PL, PR, Psa, Psv, Ppa, Ppv, Qao, Qmi, Qpul, Qtri,
            Qsys, Qlung, J)
end

function pulsatile_rhs!(du, u, par, t)
    p, sc = par
    h = pulsatile_hemodynamics(p, u, sc, t)
    du[1] = h.Qmi - h.Qao
    du[2] = h.Qao - h.Qsys
    du[3] = h.Qsys - h.Qtri + sc.Qinf(t) - sc.Qhem(t) - h.J
    du[4] = h.Qtri - h.Qpul
    du[5] = h.Qpul - h.Qlung
    du[6] = h.Qlung - h.Qmi
    du[7] = h.e.HR / 60
    extras_rhs!(du, u, PULS_O, p, sc, t, h.e, h.Psa, h.Vblood, h.J)
    return nothing
end

"""
    initial_state_pulsatile(p)

Vasos no equilíbrio do modelo lento; ventrículos no fim da diástole; fase 0.
"""
function initial_state_pulsatile(p)
    op = operating_point(p; Vb=p.Vb0)
    VL = p.VdL + (op.Ppv - p.Pth) / p.EminL
    VR = p.VdR + (op.RAP - p.Pth) / p.EminR
    Vb = op.Vsa + op.Vsv + op.Vpa + op.Vpv + VL + VR
    return vcat([VL, op.Vsa, op.Vsv, VR, op.Vpa, op.Vpv, 0.0], extras_baseline(p, Vb, op.Psa))
end

"""
    simulate_pulsatile(p, sc, tspan; u0=initial_state_pulsatile(p), kwargs...)

Resolve a malha fechada com um método para sistemas rígidos (válvulas) — `Rodas5P()`.
"""
function simulate_pulsatile(p, sc::Scenario, tspan; u0=initial_state_pulsatile(p),
                            solver=Rodas5P(), reltol=1e-6, abstol=1e-6, kwargs...)
    prob = ODEProblem(pulsatile_rhs!, u0, tspan, (p, sc))
    return solve(prob, solver; reltol, abstol, tstops=sc.tstops, kwargs...)
end

"""
    pulsatile_observables(sol, p, sc, ts)

Séries (pressões, volumes, fluxos aórtico e pulmonar) nos instantes `ts`.
"""
function pulsatile_observables(sol, p, sc, ts)
    hs = [pulsatile_hemodynamics(p, sol(t), sc, t) for t in ts]
    us = [sol(t) for t in ts]
    return (; t=collect(ts),
            PL=[h.PL for h in hs], PR=[h.PR for h in hs],
            Psa=[h.Psa for h in hs], Psv=[h.Psv for h in hs],
            Ppa=[h.Ppa for h in hs], Ppv=[h.Ppv for h in hs],
            VL=[u[1] for u in us], VR=[u[4] for u in us],
            Qao=[h.Qao for h in hs], Qpul=[h.Qpul for h in hs],
            Qsys=[h.Qsys for h in hs], HR=[h.e.HR for h in hs],
            Hct=[h.e.Hct for h in hs], Vblood=[h.Vblood for h in hs])
end
