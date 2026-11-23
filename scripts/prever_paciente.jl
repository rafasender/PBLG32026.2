# scripts/prever_paciente.jl
# Previsão para um paciente informado no terminal: qual tratamento responde melhor?
#
#   julia --project=. scripts/prever_paciente.jl               # pergunta PAM, FC e Hct
#   julia --project=. scripts/prever_paciente.jl 62 118 30     # ou passa direto (Hct em % ou fração)
#
# 1. Parte dos parâmetros calibrados no MIMIC (mediana da coorte, docs/parametros_mimic.csv) ou,
#    se o arquivo não existir, dos parâmetros de literatura.
# 2. Ajusta ao paciente o volume sanguíneo e a FC intrínseca que explicam a PAM e a FC atuais
#    (com o barorreflexo ligado): é a "hipovolemia estimada".
# 3. Simula cada tratamento sozinho por 6 h e aponta o melhor (maior DO₂ com PAM ≥ 65).
# Idade, sexo e peso não entram: o modelo não tem parâmetros que dependam deles (adulto de 70 kg).

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))

using Cardio0D
using CSV, DataFrames, Printf
include(joinpath(@__DIR__, "plot_style.jl"))

const ROOT = normpath(joinpath(@__DIR__, ".."))
const HORAS = 6.0                    # horizonte da previsão
const T_START = 600.0                # tratamento começa 10 min após o início

"Lê um número do terminal (aceita vírgula); Enter ou fim da entrada usa o padrão."
function pergunta(texto, padrao, lo, hi)
    while true
        print(texto, " [", padrao, "]: ")
        s = eof(stdin) ? "" : strip(readline())
        isempty(s) && return float(padrao)
        x = tryparse(Float64, replace(s, "," => "."))
        x !== nothing && lo <= x <= hi && return x
        println("  valor inválido: use um número entre ", lo, " e ", hi)
    end
end

# ------------------------------------------------------------------ 1. dados do paciente
println("\nPrevisão de resposta ao tratamento — modelo 0D calibrado com o MIMIC-IV\n")
if length(ARGS) >= 3
    pam, fc, hct = parse.(Float64, replace.(ARGS[1:3], "," => "."))
else
    pam = pergunta("PAM atual (mmHg)", 60, 30, 130)
    fc  = pergunta("FC atual (bpm)", 115, 40, 200)
    hct = pergunta("Hematócrito (% ou fração)", 32, 0.1, 70)
end
hct = hct > 1 ? hct / 100 : hct

# ------------------------------------------------------------------ 2. parâmetros de partida
p0 = baseline_params()
arq = joinpath(ROOT, "docs", "parametros_mimic.csv")
if isfile(arq)
    t = CSV.read(arq, DataFrame)
    razao(nome) = t.razao_mediana[findfirst(==(nome), t.parametro)]
    p0 = with_params(p0; Rsa=p0.Rsa * razao("Rsa"))
    println("\nResistência vascular: mediana calibrada nos pacientes do MIMIC (", relpath(arq, ROOT), ").")
else
    println("\nAviso: ", relpath(arq, ROOT), " não encontrado; usando parâmetros de literatura.")
end
p0 = with_params(p0; Hct_ref=clamp(hct, 0.15, 0.60))

# ------------------------------------------------------------------ 3. estado atual do paciente
print("Ajustando o modelo ao estado atual ... ")
agora = PatientData("terminal", [0.0, 60.0], [pam, pam], [fc, fc], [hct, hct], zeros(2), zeros(2), zeros(2))
fit = calibrate(p0, agora, [:Vb0, :HR]; w_hr=0.25, maxiters=80)
y0 = simulate_patient_series(fit.p, agora)
# o Hct informado é o de AGORA (já diluído pelo líquido que volta do interstício após a perda):
# corrige o Hct de antes da perda e ajusta de novo
p0 = with_params(p0; Hct_ref=clamp(p0.Hct_ref * hct / y0.Hct[1], 0.15, 0.65))
fit = calibrate(p0, agora, [:Vb0, :HR]; w_hr=0.25, maxiters=80)
p = fit.p
y0 = simulate_patient_series(p, agora)
println("ok")
sc_agora = scenario_from_patient(agora)
v_agora = mean_hemodynamics(p, settled_state(p, sc_agora), sc_agora, 0.0).Vblood
deficit = p.Vb_ref - v_agora
@printf("  modelo: PAM %.0f mmHg, FC %.0f bpm, Hct %.2f (informado: %.0f, %.0f e %.2f)\n",
        y0.PAM[end], y0.HR[end], y0.Hct[1], pam, fc, hct)
@printf("  volume de sangue agora %.0f mL: déficit de %.0f mL (%.0f%% da volemia normal)\n",
        v_agora, deficit, 100 * deficit / p.Vb_ref)
abs(y0.PAM[end] - pam) > 3 &&
    println("  atenção: o modelo não reproduz bem esse estado (erro > 3 mmHg); a previsão é menos confiável.")

# ------------------------------------------------------------------ 4. previsão dos tratamentos
println("\nSimulando ", length(TREATMENTS), " tratamentos por ", Int(HORAS), " h ...")
grade = PatientData("grade", collect(0.0:60.0:HORAS*3600), fill(NaN, Int(HORAS * 60) + 1),
                    fill(NaN, Int(HORAS * 60) + 1), fill(NaN, Int(HORAS * 60) + 1),
                    zeros(Int(HORAS * 60) + 1), zeros(Int(HORAS * 60) + 1), zeros(Int(HORAS * 60) + 1))
series = [simulate_patient_series(p, grade; sc=treatment_scenario(tr, T_START)) for tr in TREATMENTS]
fim = [(; PAM=s.PAM[end], HR=s.HR[end], Hct=s.Hct[end], DO2=s.DO2[end], CO=s.CO[end]) for s in series]
k1 = findfirst(>=(T_START + 3600), grade.t)
b = best_treatment(fim)

@printf("\n%-22s %10s %10s %8s %7s %11s %10s\n", "Tratamento", "PAM 1 h", "PAM 6 h", "FC 6 h", "Hct", "CO (L/min)", "DO₂")
for (tr, s, f) in zip(TREATMENTS, series, fim)
    @printf("%-22s %10.0f %10.0f %8.0f %7.2f %11.1f %10.0f%s\n", tr.nome, s.PAM[k1], f.PAM, f.HR, f.Hct,
            f.CO, f.DO2, tr === TREATMENTS[b] ? "  <- melhor" : "")
end
println("\nCritério: maior oferta de O₂ (DO₂, mL/min) entre os tratamentos que mantêm PAM ≥ 65 mmHg.")
println("Doses: cristaloide e hemácias em 1 h; noradrenalina em µg/kg/min contínua. Início em 10 min.")

# ------------------------------------------------------------------ 5. figura
h = grade.t ./ 3600
fim_t = T_START / 3600 + 1
pn = map([(:PAM, "PAM (mmHg)", "Pressão arterial média"), (:HR, "FC (bpm)", "Frequência cardíaca"),
          (:Hct, "hematócrito", "Hematócrito"), (:DO2, "DO₂ (mL/min)", "Oferta de O₂ aos tecidos")]) do (campo, yl, tl)
    q = plot(xlabel="horas", ylabel=yl, title=tl, legend=false, xlims=(0, HORAS + 1.8))
    fase!(q, T_START / 3600, fim_t, FUNDO.tratamento)
    for (tr, s) in zip(TREATMENTS, series)
        y = getfield(s, campo)
        plot!(q, h, y, color=COR_TRAT[tr.nome])
        campo in (:PAM, :DO2) && rotulo!(q, HORAS, y[end], tr.nome)
    end
    campo == :PAM && hline!(q, [65], color=TINTA.fraca, ls=:dot, lw=1.5)
    nomes_fases!(q, [(fim_t + 0.1, fim_t + 2.2, "← tratamento (1 h)")])
end
txt = [@sprintf("Paciente informado: PAM %.0f mmHg, FC %.0f bpm, hematócrito %.2f. Volume de sangue agora, estimado pelo modelo: %.0f mL (déficit de %.0f mL).",
                pam, fc, hct, v_agora, deficit),
       "Cada cor é um tratamento simulado sozinho a partir do mesmo estado (faixa azul = 1 h de infusão; noradrenalina contínua). " *
       "Melhor pelo critério (maior oferta de O₂ com PAM ≥ 65): " * TREATMENTS[b].nome * "."]
cab = cabecalho("Previsão de resposta ao tratamento", txt; cols=150)
hc = altura_cabecalho(txt; cols=150)
fig = joinpath(ROOT, "output", "previsao_paciente.png")
mkpath(dirname(fig))
savefig(plot(cab, plot(pn...; layout=(2, 2)); layout=grid(2, 1, heights=[hc / (hc + 760), 760 / (hc + 760)]),
             size=(1250, hc + 760), left_margin=5Plots.mm, bottom_margin=4Plots.mm), fig)
println("Figura: ", relpath(fig, ROOT))
