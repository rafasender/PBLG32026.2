# scripts/run_model.jl
# Gera as figuras e a tabela de resultados do relatório (docs/figuras/*.png, docs/resultados.md).
# Todas usam o MODELO com parâmetros de literatura (Frente B, adulto de 70 kg): nenhum dado de
# paciente entra aqui. Os dados do MIMIC entram em scripts/pipeline_mimic.jl (figuras 07 e 08).
#
#   julia --project=. scripts/run_model.jl

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))

using Cardio0D
using Printf
include(joinpath(@__DIR__, "plot_style.jl"))
const FIGDIR = pasta_figuras(normpath(joinpath(@__DIR__, "..")))

const OUT = IOBuffer()
out(args...) = (println(OUT, args...); println(args...))
br(x, d=1) = replace(@sprintf("%.*f", d, x), "." => ",")      # número com vírgula decimal

"Figura = cabeçalho explicativo em cima + painéis embaixo (no relatório, só os painéis)."
function figura(cab, paineis...; layout_paineis, largura, altura, h_cab=0.16)
    corpo = plot(paineis...; layout=layout_paineis)
    RELATORIO && return plot(corpo; size=(largura, round(Int, altura * (1 - h_cab))),
                             left_margin=6Plots.mm, bottom_margin=5Plots.mm, right_margin=4Plots.mm)
    return plot(cab, corpo; layout=grid(2, 1, heights=[h_cab, 1 - h_cab]), size=(largura, altura),
                left_margin=6Plots.mm,
                bottom_margin=5Plots.mm, right_margin=4Plots.mm)
end

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
p1 = plot(xlabel="tempo (s)", ylabel="pressão (mmHg)", title="Pressão nas câmaras em 3 batimentos",
          legend=:outerright)
# faixas de ejeção (válvula aórtica aberta)
ej = obs0.Qao .> 1.0
ejecoes = []
for i in findall(k -> ej[k] && (k == 1 || !ej[k-1]), eachindex(ej))
    j = findnext(!, ej, i); j = j === nothing ? length(ej) : j
    fase!(p1, tt[i], tt[j], FUNDO.tratamento)
    push!(ejecoes, (tt[i], tt[j], "ejeção"))
end
plot!(p1, tt, obs0.PL, color=SERIE[1], label="ventrículo esquerdo")
plot!(p1, tt, obs0.Psa, color=SERIE[2], label="aorta")
plot!(p1, tt, obs0.Ppv, color=SERIE[3], label="átrio esquerdo")
plot!(p1, tt, obs0.Ppa, color=SERIE[4], label="artéria pulmonar")
nomes_fases!(p1, ejecoes)

V, P = pv_loop(obs0, p)
selR = findall(t -> t >= obs0.t[end] - 60 / p.HR, obs0.t)
p3 = plot(xlabel="volume do ventrículo (mL)", ylabel="pressão (mmHg)",
          title="Ciclo pressão–volume (PV loop)", legend=:topright)
plot!(p3, V, P, color=SERIE[1], label="ventrículo esquerdo")
plot!(p3, obs0.VR[selR], obs0.PR[selR], color=SERIE[5], label="ventrículo direito")
EDV, ESV, Pm = maximum(V), minimum(V), maximum(P)
annotate!(p3, (EDV + ESV) / 2, -9, text("1. enchimento (mitral aberta)", :center, 8, TINTA.media))
annotate!(p3, EDV + 2, Pm * 0.45, text("2. contração\nisovolumétrica", :left, 8, TINTA.media))
annotate!(p3, (EDV + ESV) / 2, Pm + 7, text("3. ejeção (aórtica aberta)", :center, 8, TINTA.media))
annotate!(p3, ESV - 2, Pm * 0.45, text("4. relaxamento\nisovolumétrico", :right, 8, TINTA.media))
annotate!(p3, (EDV + ESV) / 2, Pm * 0.55, text(@sprintf("volume ejetado ≈ %.0f mL", EDV - ESV), :center, 9, TINTA.forte))
xlims!(p3, (ESV - 30, EDV + 35)); ylims!(p3, (-14, Pm + 15))

c1 = cabecalho("Figura 1 — Coração saudável em repouso (modelo pulsátil)",
    ["Dados: nenhum paciente; parâmetros de literatura (Frente B, adulto de 70 kg). 40 s simulados, mostrados os 3 últimos batimentos.",
     "Esquerda: pressão em cada câmara; as faixas azuis marcam a ejeção (válvula aórtica aberta), quando o ventrículo e a aorta sobem juntos.",
     @sprintf("Direita: o ciclo de um batimento. Largura = volume ejetado. Resultado: PA %s/%s mmHg, PAM %s, débito %s L/min, FE %s%%.",
              br(m0.PAS, 0), br(m0.PAD, 0), br(m0.PAM, 0), br(m0.CO, 2), br(100m0.EF, 0))]; cols=185)
savefig(figura(c1, p1, p3; layout_paineis=(1, 2), largura=1400, altura=600, h_cab=0.15),
        joinpath(FIGDIR, "01_basal_pulsatil.png"))

# ----------------------------------------------------------------------- 2. Guyton
out("## 2. Curvas de Guyton: basal × hemorragia × terapias\n")
pk = with_params(p; kf=0.0)
cases = [("Basal", p, operating_point(p), TINTA.forte, :solid),
         ("Perdeu 400 mL", pk, operating_point(pk; Vb=pk.Vb0 - 400), SERIE[2], :solid),
         ("Perdeu 400 mL e recebeu 400 mL", pk, operating_point(pk; Vb=pk.Vb0), SERIE[1], :dash),
         ("Perdeu 400 mL + vasoconstrição (resistência ×1,5)", pk, operating_point(pk; Vb=pk.Vb0 - 400), SERIE[4], :solid)]
rap = range(-3, 12; length=200)
pg = plot(xlabel="pressão no átrio direito, RAP (mmHg)", ylabel="fluxo (mL/s)",
          title="Onde o coração opera: retorno venoso × função cardíaca", xlims=(-3, 12.5),
          ylims=(0, 185), legend=:topright)
out("| Cenário | Pms (mmHg) | RAP* (mmHg) | CO* (L/min) |\n|---|---|---|---|")
pts = []
for (i, (nm, pp, op, cor, ls)) in enumerate(cases)
    Rsa = i == 4 ? 1.5 * pp.Rsa : pp.Rsa
    plot!(pg, rap, [venous_return_curve(pp, op, r; Rsa) for r in rap], color=cor, ls=ls,
          label="retorno venoso: " * nm)
    g = guyton_intersection(pp, op; Rsa)
    push!(pts, (g, cor, nm))
    out(@sprintf("| %s | %.2f | %.2f | %.2f |", nm, op.Pms, g.RAP, g.CO))
end
plot!(pg, rap, [cardiac_output_curve(p, op0, r) for r in rap], color=TINTA.fraca, ls=:dash, lw=2.5,
      label="função cardíaca (coração direito)")
for (k, (g, cor, nm)) in enumerate(pts[[1, 2, 4]])
    scatter!(pg, [g.RAP], [g.Q], color=cor, ms=8, markerstrokecolor=:white, markerstrokewidth=2, label="")
    lado = k == 3 ? :right : :left
    annotate!(pg, g.RAP + (k == 3 ? -0.3 : 0.3), g.Q + (k == 3 ? -4 : 4),
              text("débito " * br(g.CO, 2) * " L/min", lado, 8, TINTA.forte))
end
c2 = cabecalho("Figura 2 — Curvas de Guyton: por que perder sangue derruba o débito",
    ["Dados: nenhum paciente; equações de equilíbrio do modelo (retorno venoso e função cardíaca) com parâmetros de literatura. Não é simulação no tempo.",
     "Como ler: o coração opera no cruzamento da reta do retorno venoso (cheia) com a função cardíaca (tracejada cinza). O ponto marca o débito resultante.",
     "Perder 400 mL desce a reta (débito 5,25 → 4,08 L/min); repor 400 mL a devolve (tracejada azul, sobre a preta); contrair os vasos inclina a reta e piora o débito."]; cols=165)
savefig(figura(c2, pg; layout_paineis=(1, 1), largura=1100, altura=760, h_cab=0.22),
        joinpath(FIGDIR, "02_guyton.png"))
out("")

# --------------------------------------------------------- 3. hemorragia com barorreflexo
out("## 3. Hemorragia classe III (1500 mL em 10 min), com barorreflexo\n")
hem = window_rate(1500.0, 60.0, 600.0)
sc3 = Scenario(Qhem=hem, tstops=[60.0, 660.0])
r3 = mean_observables(simulate_mean(p, sc3, (0.0, 3600.0)), p, sc3, 0:10:3600)
t_min = r3.t ./ 60
fases3 = [(1, 11, "perdendo
sangue"), (11, 60, "sem tratamento: o corpo compensa")]
f3 = map([(r3.PAM, "PAM (mmHg)", "Pressão arterial média"), (r3.CO, "débito (L/min)", "Débito cardíaco"),
          (r3.HR, "FC (bpm)", "Frequência cardíaca")]) do (y, yl, tl)
    q = plot(xlabel="tempo (min)", ylabel=yl, title=tl, legend=false)
    fase!(q, 1, 11, FUNDO.perda)
    plot!(q, t_min, y, color=SERIE[1])
    nomes_fases!(q, fases3)
end
hline!(f3[1], [65], color=TINTA.fraca, ls=:dot, lw=1.5)
annotate!(f3[1], 58, 66.5, text("meta PAM 65", :right, 8, TINTA.media))
c3 = cabecalho("Figura 3 — Hemorragia de 1500 mL (classe III) e a compensação do corpo",
    ["Dados: nenhum paciente real. Paciente hipotético padrão do modelo (adulto de 70 kg, 4,7 L de sangue, valores normais de literatura), com barorreflexo ligado.",
     "Cenário escolhido: perde 1500 mL (32% do sangue; a classe III do ATLS vai de 30% a 40%) em 10 min (faixa vermelha) e depois não recebe tratamento. Os valores abaixo são o que a simulação calcula.",
     @sprintf("Resultado: a PAM cai até %s mmHg e o débito até %s L/min; a FC sobe a %s bpm. Sem tratamento, a PAM se recupera devagar (%s mmHg em 60 min), com o reflexo e o líquido que volta do interstício.",
              br(minimum(r3.PAM), 0), br(minimum(r3.CO), 2), br(maximum(r3.HR), 0), br(r3.PAM[end], 0))]; cols=195)
savefig(figura(c3, f3...; layout_paineis=(1, 3), largura=1500, altura=580, h_cab=0.17),
        joinpath(FIGDIR, "03_hemorragia_barorreflexo.png"))
out("| PAM mín (mmHg) | CO mín (L/min) | FC máx (bpm) | PAM aos 60 min |\n|---|---|---|---|")
out(@sprintf("| %.0f | %.2f | %.0f | %.0f |", minimum(r3.PAM), minimum(r3.CO), maximum(r3.HR), r3.PAM[end]))
out("")

# classes ATLS (barorreflexo ligado)
out("### Classes de hemorragia (barorreflexo ligado, PAM aos 20 min)\n")
out("| Classe | Perda (mL) | PAM (mmHg) | CO (L/min) | FC (bpm) |\n|---|---|---|---|---|")
azul = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]      # quanto mais escuro, maior a perda
pcl = plot(xlabel="tempo (min)", ylabel="PAM (mmHg)", title="PAM por classe de hemorragia",
           legend=false, xlims=(0, 70))
fase!(pcl, 1, 11, FUNDO.perda)
for (k, (cl, vol)) in enumerate((("I", 600.0), ("II", 1000.0), ("III", 1500.0), ("IV", 2000.0)))
    sc = Scenario(Qhem=window_rate(vol, 60.0, 600.0), tstops=[60.0, 660.0])
    sol = try
        simulate_mean(p, sc, (0.0, 3600.0))
    catch
        nothing
    end
    sol === nothing && continue
    o = mean_observables(sol, p, sc, 0:10:3600)
    plot!(pcl, o.t ./ 60, o.PAM, color=azul[k])
    rotulo!(pcl, 60, o.PAM[end], "classe $cl ($(Int(vol)) mL)")
    i = findfirst(>=(1200), o.t)
    out(@sprintf("| %s | %d | %.0f | %.2f | %.0f |", cl, vol, o.PAM[i], o.CO[i], o.HR[i]))
end
hline!(pcl, [65], color=TINTA.fraca, ls=:dot, lw=1.5)
annotate!(pcl, 59, 66.5, text("meta PAM 65", :right, 8, TINTA.media))
nomes_fases!(pcl, [(1, 11, "perdendo sangue"), (11, 60, "sem tratamento")])
c4 = cabecalho("Figura 4 — Quanto maior a perda, menor a pressão (classes de hemorragia do ATLS)",
    ["Dados: nenhum paciente; modelo de valores médios com parâmetros de literatura e barorreflexo ligado. Cada curva perde um volume diferente em 10 min.",
     "Como ler: cores mais escuras = perdas maiores. Limitação: na classe IV (2000 mL) o modelo ainda segura a PAM em ~70 mmHg, porque o reflexo",
     "do modelo é de uma pessoa saudável e não entra em falência. Na clínica a classe IV descompensa."]; cols=165)
savefig(figura(c4, pcl; layout_paineis=(1, 1), largura=1100, altura=680, h_cab=0.21),
        joinpath(FIGDIR, "04_classes_hemorragia.png"))
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
curto = Dict("Sem terapia" => "sem terapia", "Ringer Lactato 1500 mL" => "Ringer 1500 mL",
             "Sangue total 1000 mL" => "sangue total 1000 mL", "Noradrenalina 0,1 µg/kg/min" => "noradrenalina",
             "RL 1000 mL + noradrenalina" => "Ringer 1000 mL + nora")
fases4 = [(1, 11, "perdendo\nsangue"), (20, 35, "infusão\nde volume"), (35, 120, "após a infusão")]
f4 = map([(:PAM, "PAM (mmHg)", "Pressão arterial média"), (:CO, "débito (L/min)", "Débito cardíaco"),
          (:Hct, "hematócrito", "Hematócrito (diluição do sangue)"), (:DO2, "DO₂ (mL/min)", "Oferta de O₂ aos tecidos")]) do (campo, yl, tl)
    q = plot(xlabel="tempo (min)", ylabel=yl, title=tl, legend=false, xlims=(0, 150))
    fase!(q, 1, 11, FUNDO.perda); fase!(q, 20, 35, FUNDO.tratamento)
    vline!(q, [20], color=TINTA.media, ls=:dot, lw=1.2, label="")
    for (nm, _...) in terapias
        o = res4[nm]; y = getfield(o, campo)
        plot!(q, o.t ./ 60, y, color=COR_TRAT[nm], label=campo == :CO ? curto[nm] : "")
        campo in (:PAM, :DO2) && rotulo!(q, 120, y[end], curto[nm])
    end
    nomes_fases!(q, fases4)
end
# legenda no painel do débito (as curvas de "sem terapia" e "noradrenalina" se cruzam ali)
plot!(f4[2], legend=:topright)
annotate!(f4[3], 122, 0.351, text("  sem terapia e
  noradrenalina
  (mesma curva)", :left, 8, TINTA.forte))
hline!(f4[1], [65], color=TINTA.fraca, ls=:dot, lw=1.5)
c5 = cabecalho("Figura 5 — Comparação de tratamentos depois de uma hemorragia de 1500 mL",
    ["Dados: nenhum paciente; modelo de valores médios com parâmetros de literatura e barorreflexo ligado. O mesmo paciente simulado recebe cada tratamento.",
     "Linha do tempo: perde 1500 mL entre 1 e 11 min (faixa vermelha); o tratamento começa aos 20 min (linha pontilhada). Ringer e sangue entram em 15 min (faixa azul);",
     "a noradrenalina segue contínua até o fim. Como ler: compare as cores no fim de cada painel. A noradrenalina sobe a PAM mais rápido, mas o débito e a oferta",
     "de O₂ ficam os mais baixos; o Ringer dilui o sangue (hematócrito cai); o sangue total entrega mais O₂. Pressão alta não significa tecido bem oxigenado."]; cols=190)
savefig(figura(c5, f4...; layout_paineis=(2, 2), largura=1450, altura=1000, h_cab=0.14),
        joinpath(FIGDIR, "05_terapias.png"))

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
