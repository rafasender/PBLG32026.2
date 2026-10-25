# scripts/run_model.jl
# Gera as figuras e a tabela de resultados do relatório (docs/figuras/*.png, docs/resultados.md).
#
#   julia --project=. scripts/run_model.jl

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))

using Cardio0D
using Plots
using Printf

gr()
default(size=(900, 520), linewidth=2, titlefontsize=11, guidefontsize=10, legendfontsize=8,
        left_margin=4Plots.mm, bottom_margin=4Plots.mm)

const FIGDIR = normpath(joinpath(@__DIR__, "..", "docs", "figuras"))
mkpath(FIGDIR)
const OUT = IOBuffer()
out(args...) = (println(OUT, args...); println(args...))

p = baseline_params()

# ------------------------------------------------------------------ 1. basal pulsátil
out("## 1. Basal pulsátil (40 s, regime periódico)\n")
sc0 = Scenario(baro=false)
sol0 = simulate_pulsatile(p, sc0, (0.0, 40.0))
ts = last_beats(sol0, p; n=3)
obs0 = pulsatile_observables(sol0, p, sc0, ts)
m0 = cycle_metrics(obs0)
op0 = operating_point(p)
out("| Variável | Pulsátil | Modelo lento (alvo) |\n|---|---|---|")
out(@sprintf("| PAM (mmHg) | %.1f | %.1f |", m0.PAM, op0.Psa))
out(@sprintf("| PAS/PAD (mmHg) | %.0f/%.0f | — |", m0.PAS, m0.PAD))
out(@sprintf("| CO (L/min) | %.2f | %.2f |", m0.CO, op0.CO))
out(@sprintf("| SV (mL) | %.1f | %.1f |", m0.SV, op0.Q * 60 / p.HR))
out(@sprintf("| EDV/ESV VE (mL) | %.0f/%.0f | 120/50 |", m0.EDV, m0.ESV))
out(@sprintf("| FEVE | %.2f | 0.58 |", m0.EF))
out(@sprintf("| PVC média (mmHg) | %.1f | %.1f |", m0.CVP, op0.Psv))
out(@sprintf("| P. arterial pulmonar média (mmHg) | %.1f | %.1f |\n", m0.Ppa, op0.Ppa))

tt = obs0.t .- obs0.t[1]
p1 = plot(tt, obs0.PL, label="Ventrículo esq.", xlabel="t (s)", ylabel="Pressão (mmHg)",
          title="Pressões (últimos 3 batimentos)")
plot!(p1, tt, obs0.Psa, label="Aorta/artérias sist.")
plot!(p1, tt, obs0.Ppv, label="Átrio esq. (veias pulm.)")
plot!(p1, tt, obs0.PR, label="Ventrículo dir.", linestyle=:dash)
plot!(p1, tt, obs0.Ppa, label="Art. pulmonar", linestyle=:dash)
p2 = plot(tt, obs0.Qao, label="Aórtico", xlabel="t (s)", ylabel="Fluxo (mL/s)",
          title="Fluxos valvares")
plot!(p2, tt, obs0.Qpul, label="Pulmonar", linestyle=:dash)
V, P = pv_loop(obs0, p)
tR = obs0.t[end] - 60 / p.HR
selR = findall(t -> t >= tR, obs0.t)
p3 = plot(V, P, label="VE", xlabel="Volume (mL)", ylabel="Pressão (mmHg)", title="Loops P–V",
          legend=:topright)
plot!(p3, obs0.VR[selR], obs0.PR[selR], label="VD")
savefig(plot(p1, p2, p3; layout=(1, 3), size=(1500, 460)), joinpath(FIGDIR, "01_basal_pulsatil.png"))

# ----------------------------------------------------------------------- 2. Guyton
out("## 2. Curvas de Guyton: basal × hemorragia × terapias\n")
pk = with_params(p; kf=0.0)
cases = [("Basal", p, operating_point(p)),
         ("Hemorragia 400 mL", pk, operating_point(pk; Vb=pk.Vb0 - 400)),
         ("Hemorragia 400 mL + 400 mL de volume", pk, operating_point(pk; Vb=pk.Vb0)),
         ("Hemorragia 400 mL + vasoconstrição (TPR×1,5)", pk, operating_point(pk; Vb=pk.Vb0 - 400))]
rap = range(-3, 12; length=200)
pg = plot(xlabel="Pressão do átrio direito RAP (mmHg)", ylabel="Fluxo (mL/s)",
          title="Retorno venoso (linhas contínuas) e função cardíaca (tracejada)", xlims=(-3, 12),
          ylims=(0, 180))
cols = [:black, :red, :green, :orange]
out("| Cenário | Pms (mmHg) | RAP* (mmHg) | CO* (L/min) |\n|---|---|---|---|")
for (i, (nm, pp, op)) in enumerate(cases)
    Rsa = i == 4 ? 1.5 * pp.Rsa : pp.Rsa
    plot!(pg, rap, [venous_return_curve(pp, op, r; Rsa) for r in rap], label="VR: $nm", color=cols[i])
    g = guyton_intersection(pp, op; Rsa)
    scatter!(pg, [g.RAP], [g.Q], label="", color=cols[i], markersize=6)
    out(@sprintf("| %s | %.2f | %.2f | %.2f |", nm, op.Pms, g.RAP, g.CO))
end
op0g = operating_point(p)
plot!(pg, rap, [cardiac_output_curve(p, op0g, r) for r in rap], label="Função cardíaca (VD)",
      color=:blue, linestyle=:dash)
savefig(pg, joinpath(FIGDIR, "02_guyton.png"))
out("")

# --------------------------------------------------------- 3. hemorragia, barorreflexo
out("## 3. Hemorragia classe III (1500 mL em 10 min), barorreflexo ligado × desligado\n")
hem = window_rate(1500.0, 60.0, 600.0)
res3 = Dict()
for baro in (false, true)
    sc = Scenario(Qhem=hem, baro=baro, tstops=[60.0, 660.0])
    sol = simulate_mean(p, sc, (0.0, 3600.0))
    res3[baro] = mean_observables(sol, p, sc, 0:10:3600)
end
t_min = res3[true].t ./ 60
pa = plot(t_min, res3[true].PAM, label="Barorreflexo ON", xlabel="Tempo (min)", ylabel="PAM (mmHg)", title="PAM")
plot!(pa, t_min, res3[false].PAM, label="Barorreflexo OFF", linestyle=:dash)
hline!(pa, [65], label="PAM 65", color=:gray, linestyle=:dot)
pb = plot(t_min, res3[true].CO, label="ON", xlabel="Tempo (min)", ylabel="CO (L/min)", title="Débito cardíaco")
plot!(pb, t_min, res3[false].CO, label="OFF", linestyle=:dash)
pc = plot(t_min, res3[true].HR, label="ON", xlabel="Tempo (min)", ylabel="FC (bpm)", title="Frequência cardíaca")
plot!(pc, t_min, res3[false].HR, label="OFF", linestyle=:dash)
savefig(plot(pa, pb, pc; layout=(1, 3), size=(1500, 440)), joinpath(FIGDIR, "03_hemorragia_barorreflexo.png"))
out("| | PAM mín (mmHg) | CO mín (L/min) | FC máx (bpm) | PAM aos 60 min |\n|---|---|---|---|---|")
for baro in (false, true)
    r = res3[baro]
    out(@sprintf("| Barorreflexo %s | %.0f | %.2f | %.0f | %.0f |", baro ? "ON" : "OFF",
                 minimum(r.PAM), minimum(r.CO), maximum(r.HR), r.PAM[end]))
end
out("")

# classes ATLS (barorreflexo ligado)
out("### Classes de hemorragia (barorreflexo ligado, PAM aos 20 min)\n")
out("| Classe | Perda (mL) | PAM (mmHg) | CO (L/min) | FC (bpm) |\n|---|---|---|---|---|")
pcl = plot(xlabel="Tempo (min)", ylabel="PAM (mmHg)", title="PAM por classe de hemorragia (barorreflexo ON)")
for (cl, vol) in (("I", 600.0), ("II", 1000.0), ("III", 1500.0), ("IV", 2000.0))
    sc = Scenario(Qhem=window_rate(vol, 60.0, 600.0), tstops=[60.0, 660.0])
    sol = try
        simulate_mean(p, sc, (0.0, 3600.0))
    catch
        nothing
    end
    sol === nothing && continue
    o = mean_observables(sol, p, sc, 0:10:3600)
    plot!(pcl, o.t ./ 60, o.PAM, label="Classe $cl ($(Int(vol)) mL)")
    i = findfirst(>=(1200), o.t)
    out(@sprintf("| %s | %d | %.0f | %.2f | %.0f |", cl, vol, o.PAM[i], o.CO[i], o.HR[i]))
end
hline!(pcl, [65], label="PAM 65", color=:gray, linestyle=:dot)
savefig(pcl, joinpath(FIGDIR, "04_classes_hemorragia.png"))
out("")

# ---------------------------------------------------------------------- 4. terapias
out("## 4. Terapias após hemorragia de 1500 mL (infusão a partir de t = 20 min, 15 min de duração)\n")
terapias = [("Sem terapia", 0.0, 0.0, 0.0),
            ("Ringer Lactato 1500 mL", 1500.0, 0.0, 0.0),
            ("Sangue total 1000 mL", 1000.0, 0.40, 0.0),
            ("Noradrenalina 0,1 µg/kg/min", 0.0, 0.0, 0.1),
            ("RL 1000 mL + noradrenalina", 1000.0, 0.0, 0.1)]
res4 = Dict()
for (nm, vol, hct, dose) in terapias
    sc = Scenario(Qhem=hem, Qinf=window_rate(vol, 1200.0, 900.0), hct_inf=hct,
                  dose=t -> (t >= 1200.0 ? dose : 0.0), tstops=[60.0, 660.0, 1200.0, 2100.0])
    sol = simulate_mean(p, sc, (0.0, 7200.0))
    res4[nm] = mean_observables(sol, p, sc, 0:30:7200)
end
"Minutos após o início da terapia até a PAM atingir `thr` (NaN se nunca atingir)."
function first_time_above(o, thr; after=1200.0)
    i = findfirst(k -> o.t[k] >= after && o.PAM[k] >= thr, eachindex(o.t))
    return i === nothing ? NaN : (o.t[i] - after) / 60
end
out("| Terapia | Tempo até PAM ≥ 85 (min) | PAM 40 min (mmHg) | PAM 100 min | CO 100 min (L/min) | Hct 100 min | μ rel. | DO₂ 100 min (mL/min) | Vol. intersticial 100 min (mL) |\n|---|---|---|---|---|---|---|---|---|")
for (nm, _...) in terapias
    o = res4[nm]
    i40 = findfirst(>=(2400), o.t); i100 = findfirst(>=(6000), o.t)
    out(@sprintf("| %s | %s | %.0f | %.0f | %.2f | %.3f | %.2f | %.0f | %.0f |", nm,
                 isnan(first_time_above(o, 85)) ? "não atinge" : @sprintf("%.1f", first_time_above(o, 85)),
                 o.PAM[i40], o.PAM[i100], o.CO[i100], o.Hct[i100], o.mu_rel[i100], o.DO2[i100], o.Vis[i100]))
end
out("")
cl4 = Dict("Sem terapia" => :gray, "Ringer Lactato 1500 mL" => :blue, "Sangue total 1000 mL" => :red,
           "Noradrenalina 0,1 µg/kg/min" => :green, "RL 1000 mL + noradrenalina" => :purple)
f4 = [plot(xlabel="Tempo (min)", ylabel=yl, title=tl) for (yl, tl) in
      (("PAM (mmHg)", "PAM"), ("CO (L/min)", "Débito cardíaco"), ("Hct", "Hematócrito"),
       ("DO₂ (mL/min)", "Oferta de O₂"))]
for (nm, _...) in terapias
    o = res4[nm]
    for (k, y) in enumerate((o.PAM, o.CO, o.Hct, o.DO2))
        plot!(f4[k], o.t ./ 60, y, label=k == 1 ? nm : "", color=cl4[nm])
    end
end
hline!(f4[1], [85], label="", color=:gray, linestyle=:dot)
savefig(plot(f4...; layout=(2, 2), size=(1200, 800)), joinpath(FIGDIR, "05_terapias.png"))

# ------------------------------------------------------------ 5. PEEP (efeito da intubação)
out("## 5. Efeito da pressão intratorácica (PEEP) no retorno venoso\n")
out("| P_th (mmHg) | PAM (mmHg) | CO (L/min) | RAP (mmHg) |\n|---|---|---|---|")
for pth in (-4.0, -2.0, 0.0, 3.7)
    pp = with_params(p; Pth=pth)
    op = operating_point(pp)
    out(@sprintf("| %.1f | %.1f | %.2f | %.2f |", pth, op.Psa, op.CO, op.RAP))
end
out("")

write(joinpath(@__DIR__, "..", "docs", "resultados.md"), String(take!(OUT)))
println("\n[OK] Figuras em docs/figuras/ e tabelas em docs/resultados.md")
