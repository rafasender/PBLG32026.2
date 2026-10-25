# src/analysis.jl
# Métricas por ciclo, PV loop e curvas de Guyton (retorno venoso × função cardíaca).

"Média temporal por regra do trapézio de `y(t)` em toda a janela `ts`."
function time_mean(ts, y)
    length(ts) < 2 && return float(y[1])
    s = 0.0
    for i in 1:length(ts)-1
        s += 0.5 * (y[i] + y[i+1]) * (ts[i+1] - ts[i])
    end
    return s / (ts[end] - ts[1])
end

"""
    cycle_metrics(obs)

Métricas de uma janela pulsátil (de preferência vários batimentos em regime):
PAM (média temporal), PAS, PAD, CO [L/min] (média do fluxo aórtico), SV [mL],
EDV/ESV/EF do VE, RAP/CVP médios.
"""
function cycle_metrics(obs)
    t = obs.t
    PAM = time_mean(t, obs.Psa)
    CO = mls_to_lpm(time_mean(t, obs.Qao))
    HR = time_mean(t, obs.HR)
    EDV = maximum(obs.VL); ESV = minimum(obs.VL)
    return (; PAM, PAS=maximum(obs.Psa), PAD=minimum(obs.Psa), CO, HR,
            SV=CO * 1000 / HR, EDV, ESV, EF=(EDV - ESV) / EDV,
            CVP=time_mean(t, obs.Psv), Ppa=time_mean(t, obs.Ppa),
            LAP=time_mean(t, obs.Ppv), PAM_pulm=time_mean(t, obs.Ppa))
end

"Janela dos últimos `n` batimentos de uma solução pulsátil (grade de `dt` s)."
function last_beats(sol, p; n=3, dt=0.005)
    T = 60 / p.HR
    tend = sol.t[end]
    return range(tend - n * T, tend; step=dt)
end

# --------------------------------------------------------------------- curvas de Guyton
"""
    venous_return_curve(p, op, RAP)

Retorno venoso [mL/s] em função de `RAP` (equação E10) com a `Pms` do ponto de operação `op`:
`VR = (Pms − RAP) / (Rv + TPR·AC/(AC+VC))`. Para `RAP > Pms` o fluxo seria reverso: devolve 0.
"""
function venous_return_curve(p, op, RAP; Rsa=p.Rsa)
    RVR = p.Rsv + Rsa * p.Csa / (p.Csa + p.Csv)
    return max(0.0, (op.Pms - RAP) / RVR)
end

"""
    cardiac_output_curve(p, op, RAP; Ppa=op.Ppa)

Função cardíaca do coração direito (E14): `CO = HR·(C_D·(RAP−P_th) − C_S·(P_pa+dPes−P_th))`,
limitada entre 0 e o platô de `EDVmax` (§4.3). Pós-carga (`Ppa`) fixa no ponto de operação.
"""
function cardiac_output_curve(p, op, RAP; Ppa=op.Ppa, EmaxR=p.EmaxR)
    k = p.HR / 60
    CDR = 1 / p.EminR; CSR = 1 / EmaxR
    after = CSR * (Ppa + p.dPesR - p.Pth)
    SV = min(CDR * (RAP - p.Pth), p.EDVmaxR - p.VdR) - after
    return k * max(0.0, SV)
end

"""
    guyton_intersection(p, op; kwargs...)

Interseção das curvas de VR e de função cardíaca no plano (RAP, fluxo): o ponto de
funcionamento do §5.4. Devolve `(RAP, Q, CO)`.
"""
function guyton_intersection(p, op; Rsa=p.Rsa, Ppa=op.Ppa, EmaxR=p.EmaxR)
    f(RAP) = venous_return_curve(p, op, RAP; Rsa) - cardiac_output_curve(p, op, RAP; Ppa, EmaxR)
    RAP = bisect(f, p.Pth + 0.5, op.Pms)
    Q = venous_return_curve(p, op, RAP; Rsa)
    return (; RAP, Q, CO=mls_to_lpm(Q))
end

"""
    pv_loop(obs, p; beat=-1)

Pares `(V, P)` do VE em um batimento em regime (por padrão o último), útil para o gráfico.
"""
function pv_loop(obs, p; dt=0.005)
    T = 60 / p.HR
    t1 = obs.t[end]
    sel = findall(t -> t >= t1 - T, obs.t)
    return obs.VL[sel], obs.PL[sel]
end
