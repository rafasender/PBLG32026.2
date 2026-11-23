# scripts/figuras_mimic.jl
# Figuras 07 (calibração e teste) e 08 (previsão dos tratamentos) dos pacientes do MIMIC-IV,
# refeitas a partir dos parâmetros já calibrados, SEM recalibrar. É chamado no fim do
# pipeline_mimic.jl, mas também pode rodar sozinho:
#
#   julia --project=. scripts/figuras_mimic.jl [arquivo_de_calibracao.csv]
#
# Lê data/processed/mimic_calibracao.csv (Rsa, Vb0, HR, Hct_ref por paciente) e
# data/processed/mimic_PXX.csv (série de cada paciente). Ambos são locais e ignorados pelo git.

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))

using Cardio0D
using CSV, DataFrames, Printf
include(joinpath(@__DIR__, "plot_style.jl"))

const ROOT = normpath(joinpath(@__DIR__, ".."))
const PROC = joinpath(ROOT, "data", "processed")
const FIGDIR = pasta_figuras(ROOT)
const FIT_H, HOURS, T_TRAT = 6.0, 12.0, 3600.0        # os mesmos do pipeline_mimic.jl
br(x, d=1) = replace(@sprintf("%.*f", d, x), "." => ",")

arq = isempty(ARGS) ? joinpath(PROC, "mimic_calibracao.csv") : ARGS[1]
cal = CSV.read(arq, DataFrame)
p_lit = baseline_params()

"Trechos contínuos `(início, fim)` [h] em que `v > 0`."
function trechos(t, v)
    out = Tuple{Float64,Float64}[]
    i = 1
    while i <= length(v)
        if v[i] > 0
            j = i
            while j < length(v) && v[j+1] > 0
                j += 1
            end
            push!(out, (t[i] / 3600, (t[j] + 60) / 3600))
            i = j + 1
        else
            i += 1
        end
    end
    return out
end

pacientes = []
for r in eachrow(cal)
    pd = load_patient(joinpath(PROC, "mimic_$(r.lab).csv"); id=r.lab)
    p = with_params(p_lit; Hct_ref=r.Hct_ref, Rsa=r.Rsa, Vb0=r.Vb0, HR=r.HR)
    y = simulate_patient_series(p, pd)
    test = findall(>(FIT_H * 3600), pd.t)
    ultima = filter(!isnan, pd.PAM[pd.t .<= FIT_H * 3600])[end]
    prev = [simulate_patient_series(p, pd; sc=treatment_scenario(tr, T_TRAT)) for tr in TREATMENTS]
    fim = [(; PAM=s.PAM[end], DO2=s.DO2[end]) for s in prev]
    push!(pacientes, (; lab=r.lab, pd, y, prev, fim, melhor=best_treatment(fim),
                      e_teste=rmse(y.PAM[test], pd.PAM[test]),
                      e_palpite=rmse(fill(ultima, length(test)), pd.PAM[test])))
    println("    ", r.lab, ": figuras simuladas")
end
N = length(pacientes)

# ------------------------------------------------------------------ figura 07: ajuste e teste
blocos = map(pacientes) do q
    h = q.pd.t ./ 3600
    a = plot(ylabel="PAM (mmHg)", legend=(q === pacientes[1] ? :bottomleft : false),
             title=@sprintf("%s: erro no teste %s mmHg (palpite \"fica igual\": %s)", q.lab,
                            br(q.e_teste), br(q.e_palpite)), titlefontsize=10, xlims=(0, HOURS))
    fase!(a, 0, FIT_H, FUNDO.ajuste)
    vline!(a, [T_TRAT / 3600], color=TINTA.media, ls=:dot, lw=1.2, label="")
    plot!(a, h, q.y.PAM, color=SERIE[1], label="modelo")
    scatter!(a, h, q.pd.PAM, color=TINTA.forte, ms=4, markerstrokewidth=0, label="PAM medida")
    nomes_fases!(a, [(0, FIT_H, "ajuste: o modelo aprende"), (FIT_H, HOURS, "teste: o modelo prevê sozinho")])
    # faixa de tratamentos recebidos
    vols = (sum(q.pd.rl_mlh), sum(q.pd.blood_mlh)) .* 60 ./ 3600
    b = plot(xlabel="horas desde 1 h antes do 1º tratamento", legend=false, xlims=(0, HOURS),
             ylims=(0.4, 3.6), yticks=([1, 2, 3], ["noradrenalina", @sprintf("hemácias %.0f mL", vols[2]),
                                                   @sprintf("cristaloide %.0f mL", vols[1])]),
             grid=false, ytickfontsize=7)
    fase!(b, 0, FIT_H, FUNDO.ajuste)
    for (k, serie, cor) in ((3, q.pd.rl_mlh, SERIE[1]), (2, q.pd.blood_mlh, SERIE[2]), (1, q.pd.dose, SERIE[3]))
        for (t0, t1) in trechos(q.pd.t, serie)
            plot!(b, [t0, t1], [k, k], color=cor, lw=7)
        end
    end
    plot(a, b; layout=grid(2, 1, heights=[0.74, 0.26]), left_margin=14Plots.mm)
end
linhas = ceil(Int, N / 2)
txt7 = [@sprintf("Dados: MIMIC-IV, %d internações de UTI com CID de choque hipovolêmico ou traumático (sem sepse) que receberam noradrenalina, cristaloide ou hemácias. ", N) *
     "PAM arterial invasiva e FC vêm do chartevents, os tratamentos do inputevents e o hematócrito do labevents. Cada painel é um paciente (anônimo).",
     "Como ler: pontos = PAM medida (cerca de 1 por hora); linha azul = modelo. Na faixa cinza (6 h) o modelo é ajustado ao paciente (3 números: resistência, " *
     "volume de sangue e FC intrínseca). Depois ele prevê sozinho, sabendo só os tratamentos (barras embaixo). A linha pontilhada marca o 1º tratamento.",
     "O erro no teste é comparado com um palpite ingênuo, \"a PAM fica igual à última medida\": o modelo só é útil para prever se errar menos que isso."]
c7 = cabecalho("Figura 7 — O modelo calibrado com pacientes reais (MIMIC-IV, choque hipovolêmico)", txt7; cols=175)
h7 = altura_cabecalho(txt7; cols=175)
corpo7 = plot(blocos...; layout=(linhas, 2))
alt7 = h7 + (RELATORIO ? 300 : 400) * linhas
savefig(RELATORIO ? plot(corpo7; size=(1500, alt7 - h7), left_margin=6Plots.mm, bottom_margin=4Plots.mm) :
        plot(c7, corpo7; layout=grid(2, 1, heights=[h7 / alt7, 1 - h7 / alt7]), size=(1500, alt7),
             left_margin=6Plots.mm, bottom_margin=4Plots.mm),
        joinpath(FIGDIR, "07_calibracao_mimic.png"))

# ------------------------------------------------------------------ figura 08: previsão
ex = pacientes[1]
h = ex.pd.t ./ 3600
fases8 = [(2, 5.4, "← tratamento (1 h)")]
ta, tb = map([(:PAM, "PAM (mmHg)", "PAM prevista para $(ex.lab)"),
              (:DO2, "DO₂ (mL/min)", "Oferta de O₂ prevista para $(ex.lab)")]) do (campo, yl, tl)
    q = plot(xlabel="horas", ylabel=yl, title=tl, legend=false, xlims=(0, HOURS + 3.2))
    fase!(q, 1, 2, FUNDO.tratamento)
    for (tr, s) in zip(TREATMENTS, ex.prev)
        y = getfield(s, campo)
        plot!(q, h, y, color=COR_TRAT[tr.nome])
        rotulo!(q, HOURS, y[end], tr.nome)
    end
    campo == :PAM && hline!(q, [65], color=TINTA.fraca, ls=:dot, lw=1.5)
    nomes_fases!(q, fases8)
end
xs = 1:N
pa, pb = map([(:PAM, "PAM ao fim (mmHg)", "PAM prevista em 12 h, todos os pacientes"),
              (:DO2, "DO₂ ao fim (mL/min)", "Oferta de O₂ prevista em 12 h (círculo preto = melhor)")]) do (campo, yl, tl)
    pl = plot(xlabel="paciente", ylabel=yl, title=tl, xticks=(xs, [q.lab for q in pacientes]),
              legend=(campo == :PAM ? :outerright : false), xlims=(0.4, N + 0.6))
    melhores = [getfield(q.fim[q.melhor], campo) for q in pacientes]
    dx = [(q.melhor - 2.5) * 0.14 for q in pacientes]
    scatter!(pl, xs .+ dx, melhores, color=:white, markerstrokecolor=TINTA.forte,
             markerstrokewidth=2, ms=13, label=(campo == :PAM ? "melhor (critério)" : ""))
    for (j, tr) in enumerate(TREATMENTS)
        v = [getfield(q.fim[j], campo) for q in pacientes]
        scatter!(pl, xs .+ (j - 2.5) * 0.14, v, color=COR_TRAT[tr.nome], ms=7, markerstrokewidth=0,
                 label=tr.nome)
    end
    campo == :PAM && hline!(pl, [65], color=TINTA.fraca, ls=:dot, lw=1.5, label="PAM 65")
    pl
end
cont = [count(q -> q.melhor == j, pacientes) for j in eachindex(TREATMENTS)]
txt8 = ["Dados: os mesmos pacientes da Figura 7, cada um com seus parâmetros calibrados. Aqui cada tratamento é simulado SOZINHO (não é o que o paciente recebeu): " *
     "cristaloide 1000 mL em 1 h, hemácias 600 mL em 1 h ou noradrenalina 0,1 µg/kg/min contínua, todos começando 1 h após o início da janela (faixa azul).",
     "Em cima: o primeiro paciente ao longo do tempo. Embaixo: o valor ao fim de 12 h para todos. Critério de melhor: maior oferta de O₂ entre os tratamentos " *
     "que mantêm PAM ≥ 65 mmHg. Contagem: " * join(["$(tr.nome) $(cont[j])" for (j, tr) in enumerate(TREATMENTS)], "; ") * ".",
     "Atenção: o que o paciente teria feito com outro tratamento nunca é observado. A confiança nessas previsões vem do teste da Figura 7."]
c8 = cabecalho("Figura 8 — Previsão: qual tratamento responderia melhor em cada paciente", txt8; cols=175)
h8 = altura_cabecalho(txt8; cols=175)
corpo8 = plot(ta, tb, pa, pb; layout=grid(2, 2, widths=[0.5, 0.5]))
alt8 = h8 + 950
savefig(RELATORIO ? plot(corpo8; size=(1500, alt8 - h8), left_margin=6Plots.mm, bottom_margin=5Plots.mm, right_margin=6Plots.mm) :
        plot(c8, corpo8; layout=grid(2, 1, heights=[h8 / alt8, 1 - h8 / alt8]), size=(1500, alt8),
             left_margin=6Plots.mm, bottom_margin=5Plots.mm, right_margin=6Plots.mm),
        joinpath(FIGDIR, "08_previsao_tratamentos.png"))
println("[OK] figuras 07 e 08 em ", relpath(FIGDIR, ROOT))
