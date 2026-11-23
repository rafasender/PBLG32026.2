# scripts/pipeline_mimic.jl
# Pipeline com dados REAIS do MIMIC-IV, em três etapas por paciente:
#   1. calibrar  o modelo nas primeiras FIT_H horas da janela (PAM e FC medidas);
#   2. validar   prevendo as horas seguintes com os tratamentos que o paciente recebeu de fato;
#   3. prever    a resposta a cada tratamento isolado (cristaloide, hemácias, noradrenalina).
#
#   julia --project=. scripts/pipeline_mimic.jl
#
# Arquivos esperados (NÃO versionados, ver .gitignore) em data/raw/mimic/:
#   icustays.csv.gz, diagnoses_icd.csv.gz, inputevents.csv.gz, labevents.csv.gz
# e o chartevents em data/raw/mimic/chartevents.csv.gz ou data/Dados Mimic/chartevents.csv-001.gz.
#
# Coorte = interseção de quatro tabelas:
#   diagnoses_icd  choque hipovolêmico/traumático (HYPOVOLEMIC_ICD), excluindo sepse (SEPSIS_ICD)
#   inputevents    recebeu noradrenalina OU cristaloide (RL, NaCl 0,9%, ...) OU hemácias
#   chartevents    PAM arterial invasiva e FC nas duas metades da janela
#   labevents      Hct, quando houver (não é critério)
# Janela: [início do primeiro desses tratamentos − 1 h, + 12 h]; o tratamento tem de aparecer na
# janela (≥ 250 mL de volume ou noradrenalina).
#
# O barorreflexo fica LIGADO e a FC é gerada pelo modelo: é a mesma física na calibração e na
# previsão (em que não existe FC medida). A FC medida entra no erro da calibração.
# Resultados agregados vão para docs/resultados_mimic.md (internações anônimas: P01, P02, ...).

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))

using Cardio0D
using CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW  = joinpath(ROOT, "data", "raw", "mimic")
const PROC = joinpath(ROOT, "data", "processed")
const FIGDIR = joinpath(ROOT, "docs", "figuras")
mkpath(PROC); mkpath(FIGDIR)

# ------------------------------------------------------------------ configuração
const PRE_H, HOURS, FIT_H = 1.0, 12.0, 6.0   # janela [t_trat − 1 h, + 12 h]; ajuste nas 6 h iniciais
# O MIMIC-IV registra sinais vitais ~de hora em hora: exigir PAM invasiva e FC nas duas metades.
const MIN_MAP_HALF, MIN_HR_HALF, MIN_VOL = 5, 3, 250.0
const NMAX = 10                      # nº de pacientes calibrados
const DT = 60.0
const NAMES = [:Rsa, :Vb0, :HR]      # resistência, volume sanguíneo inicial, FC intrínseca
const W_HR = 0.25                    # 10 bpm de erro na FC pesam como 5 mmHg na PAM
const LAMBDA = 5.0                   # regularização: afastar um parâmetro ×2,7 custa 5 mmHg²
const MAXITERS = 100
# faixas fisiológicas (razão em relação à literatura): resistência ÷3 a ×3 (vasoplegia a
# vasoconstrição intensa), volume de sangue 60–120% (perda de 40% a sobrecarga), FC intrínseca 50–120 bpm
const LIMITES = [(1 / 3, 3.0), (0.6, 1.2), (50 / 75, 120 / 75)]

# Tratamentos comparados na previsão (TREATMENTS, src/treatments.jl), todos iniciados 1 h após o
# início da janela: o instante em que o paciente recebeu o primeiro tratamento real.
const T_TRAT = PRE_H * 3600

chart_path = first(filter(isfile, [joinpath(RAW, "chartevents.csv.gz"),
                                   joinpath(ROOT, "data", "Dados Mimic", "chartevents.csv-001.gz")]))
ids = MIMIC_ITEMIDS
const NE = ids.input.norepinephrine
const VOLUME_ITEMS = (CRYSTALLOID_ITEMS..., BLOOD_ITEMS...)
const THERAPY_ITEMS = (NE, VOLUME_ITEMS...)
tempo(f) = (t = time(); r = f(); (r, time() - t))
funil = String[]                     # contagens da coorte, para o relatório
registra(s) = (println("    ", s); push!(funil, s))

# ------------------------------------------------------------------ 1. diagnósticos
println("1/6 diagnoses_icd + icustays ...")
icu = CSV.read(joinpath(RAW, "icustays.csv.gz"), DataFrame; select=[:subject_id, :hadm_id, :stay_id])
dx = CSV.read(joinpath(RAW, "diagnoses_icd.csv.gz"), DataFrame;
              select=[:hadm_id, :icd_code, :icd_version], types=Dict(:icd_code => String))
h_hypo = hadm_with_icd(dx, HYPOVOLEMIC_ICD)
h_coorte = setdiff(h_hypo, hadm_with_icd(dx, SEPSIS_ICD))
stays_dx = Set(r.stay_id for r in eachrow(icu) if r.hadm_id in h_coorte)
registra(@sprintf("%d internações hospitalares com choque hipovolêmico/traumático; %d sem sepse; %d internações de UTI",
                  length(h_hypo), length(h_coorte), length(stays_dx)))

# ------------------------------------------------------------------ 2. tratamentos
println("2/6 inputevents ...")
inp, s = tempo() do
    d = CSV.read(joinpath(RAW, "inputevents.csv.gz"), DataFrame;
                 select=[:stay_id, :starttime, :endtime, :itemid, :amount, :rate, :statusdescription],
                 types=Dict(:starttime => DateTime, :endtime => DateTime))
    filter!(r -> r.statusdescription != "Rewritten" && r.stay_id in stays_dx, d)
    return d
end
tr = filter(r -> r.itemid in THERAPY_ITEMS, inp)
t_first = combine(groupby(tr, :stay_id), :starttime => minimum => :t_first)
cand = Set(t_first.stay_id)
n_grupo(items) = length(unique(tr.stay_id[[x in items for x in tr.itemid]]))
registra(@sprintf("%d com algum tratamento: noradrenalina %d, cristaloide %d, hemácias %d (%.0f s)",
                  length(cand), n_grupo((NE,)), n_grupo(CRYSTALLOID_ITEMS), n_grupo(BLOOD_ITEMS), s))

# ------------------------------------------------------------------ 3. chartevents (streaming)
cache = joinpath(PROC, "mimic_chart_hipovolemico.csv")
chart_all = if isfile(cache)
    println("3/6 chartevents: usando cache ", relpath(cache, ROOT))
    CSV.read(cache, DataFrame; types=Dict(:charttime => DateTime))
else
    println("3/6 chartevents em streaming (2–3 GB comprimidos, alguns minutos) ...")
    vit = (ids.chart.HR, ids.chart.MAP_invasive, ids.chart.MAP_nibp)
    d, s = tempo() do
        stream_mimic_table(chart_path,
            [:stay_id => Int, :charttime => DateTime, :itemid => Int, :valuenum => Float64],
            r -> r.itemid in vit && r.stay_id in cand; verbose=true)
    end
    @printf("    %d linhas (%.0f s)\n", nrow(d), s)
    CSV.write(cache, d)
    d
end
n_rem = filter_physiological!(chart_all)
println("    artefatos removidos (limites fisiológicos): ", n_rem)
registra(@sprintf("%d com sinais vitais no arquivo de chartevents disponível", length(unique(chart_all.stay_id))))

# ------------------------------------------------------------------ 4. seleção das janelas
println("4/6 seleção de janelas ...")
gb_chart = groupby(chart_all, :stay_id)
gb_inp = groupby(inp, :stay_id)
tfirst_of = Dict(zip(t_first.stay_id, t_first.t_first))
elig = NamedTuple[]
for stay in sort(unique(chart_all.stay_id))
    t0 = tfirst_of[stay] - Minute(round(Int, PRE_H * 60)); t1 = t0 + Minute(round(Int, HOURS * 60))
    tm = t0 + Minute(round(Int, FIT_H * 60))
    ch = filter(r -> t0 <= r.charttime <= t1, gb_chart[(stay_id=stay,)])
    n_in(item, a, b) = count(r -> r.itemid == item && a <= r.charttime < b, eachrow(ch))
    ok = n_in(ids.chart.MAP_invasive, t0, tm) >= MIN_MAP_HALF && n_in(ids.chart.MAP_invasive, tm, t1 + Second(1)) >= MIN_MAP_HALF &&
         n_in(ids.chart.HR, t0, tm) >= MIN_HR_HALF && n_in(ids.chart.HR, tm, t1 + Second(1)) >= MIN_HR_HALF
    ok || continue
    w = filter(r -> t0 <= r.starttime <= t1, gb_inp[(stay_id=stay,)])
    vol = sum(skipmissing(w.amount[[x in VOLUME_ITEMS for x in w.itemid]]); init=0.0)
    ne = any(==(NE), w.itemid)
    (vol >= MIN_VOL || ne) || continue
    push!(elig, (; stay, t0, t1, vol, ne))
end
isempty(elig) && error("nenhuma internação atende aos critérios da coorte")
registra(@sprintf("%d com PAM invasiva (≥ %d) e FC (≥ %d) em cada metade da janela e tratamento na janela (≥ %d mL ou noradrenalina)",
                  length(elig), MIN_MAP_HALF, MIN_HR_HALF, Int(MIN_VOL)))
sel = elig[1:min(NMAX, length(elig))]
registra(@sprintf("%d calibradas (as de menor stay_id), %d delas com noradrenalina na janela",
                  length(sel), count(x -> x.ne, sel)))

# ------------------------------------------------------------------ 5. labevents (só Hct, só as internações escolhidas)
println("5/6 labevents (Hct) em streaming ...")
hadm_of = Dict(r.stay_id => r.hadm_id for r in eachrow(icu))
hadms = Set(hadm_of[x.stay] for x in sel)
labs_all, s = tempo() do
    stream_mimic_table(joinpath(RAW, "labevents.csv.gz"),
        [:hadm_id => Int, :charttime => DateTime, :itemid => Int, :valuenum => Float64],
        r -> r.itemid == ids.lab.hct && r.hadm_id in hadms; verbose=true)
end
@printf("    %d medidas de Hct (%.0f s)\n", nrow(labs_all), s)

# ------------------------------------------------------------------ 6. calibrar, validar, prever
println("6/6 calibração (", Int(FIT_H), " h), validação (", Int(HOURS - FIT_H), " h) e previsão ...")
p_lit = baseline_params()

"Previsão de um tratamento isolado: PAM, FC, Hct e DO₂ ao fim da janela."
function prever(p, pd, tr)
    y = simulate_patient_series(p, pd; sc=treatment_scenario(tr, T_TRAT))
    return (; PAM=y.PAM[end], HR=y.HR[end], Hct=y.Hct[end], DO2=y.DO2[end])
end

"Calibra, valida e prevê um paciente (independente dos demais: roda em paralelo)."
function processa(k, x)
    lab = @sprintf("P%02d", k)
    ch = filter(r -> x.t0 <= r.charttime <= x.t1, gb_chart[(stay_id=x.stay,)])
    ip = DataFrame(gb_inp[(stay_id=x.stay,)])
    lb = filter(r -> r.hadm_id == hadm_of[x.stay] && x.t0 - Hour(6) <= r.charttime <= x.t1, labs_all)
    pd = build_patient(x.stay, ch, ip, lb; dt=DT, t0=x.t0, hours=HOURS)
    write_patient_csv(joinpath(PROC, "mimic_$(lab).csv"), pd)     # id interno só no nome local, ignorado pelo git
    # parâmetros de partida: literatura + Hct medido do próprio paciente (quando houver)
    hct = filter(!isnan, pd.Hct)
    p0 = isempty(hct) ? p_lit : with_params(p_lit; Hct_ref=clamp(hct[1], 0.15, 0.60))
    fitd = truncate_patient(pd, FIT_H * 3600)
    test = findall(>(FIT_H * 3600), pd.t)
    fit = calibrate(p0, fitd, NAMES; λ=LAMBDA, w_hr=W_HR, maxiters=MAXITERS, limites=LIMITES)
    # validação: o modelo calibrado prevê as horas seguintes com os tratamentos reais
    y = simulate_patient_series(fit.p, pd)
    y_lit = simulate_patient_series(p0, pd)
    last_obs = filter(!isnan, fitd.PAM)[end]                     # palpite ingênuo: "fica igual"
    r = (; lab, n_fit=count(!isnan, fitd.PAM), n_test=count(!isnan, pd.PAM[test]),
         vol_cr=sum(pd.rl_mlh) * DT / 3600, vol_bl=sum(pd.blood_mlh) * DT / 3600,
         dose=maximum(pd.dose), pam_min=minimum(filter(!isnan, pd.PAM)),
         r_lit=rmse(y_lit.PAM[test], pd.PAM[test]), r_fit=fit.rmse,
         r_test=rmse(y.PAM[test], pd.PAM[test]),
         r_naive=rmse(fill(last_obs, length(test)), pd.PAM[test]),
         r_hr=rmse(y.HR[test], pd.HR[test]),
         Rsa=fit.θ[1] / p0.Rsa, Vb0=fit.θ[2] / p0.Vb0, HR0=fit.θ[3],
         co=sum(y.CO) / length(y.CO))
    pr = [prever(fit.p, pd, tr) for tr in TREATMENTS]
    @printf("    %s: RMSE ajuste %.1f | teste %.1f (palpite %.1f, literatura %.1f) mmHg | melhor: %s\n",
            lab, r.r_fit, r.r_test, r.r_naive, r.r_lit, TREATMENTS[best_treatment(pr)].nome)
    return (; r, pr, fit, fitd,
            linha=(; lab, Rsa=fit.θ[1], Vb0=fit.θ[2], HR=fit.θ[3], Hct_ref=p0.Hct_ref))
end

# os pacientes são independentes: um por thread (julia -t N). Os índices agrupados já foram
# construídos na etapa 4, então a leitura concorrente de gb_chart/gb_inp é segura.
println("    ", Threads.nthreads(), " threads")
res = Vector{Any}(nothing, length(sel))
Threads.@threads for k in eachindex(sel)
    res[k] = try
        processa(k, sel[k])
    catch e
        println("    ", @sprintf("P%02d", k), " falhou: ", e); nothing
    end
end
ok_res = filter(!isnothing, res)
isempty(ok_res) && error("nenhum paciente calibrado")
rows = [q.r for q in ok_res]
preds = [q.pr for q in ok_res]
fits = [q.linha for q in ok_res]
# identificabilidade (1º paciente calibrado): sensibilidade média |S| de cada parâmetro
S, _ = sensitivity_matrix(ok_res[1].fit.p, ok_res[1].fitd, [:Rsa, :Vb0, :HR, :EmaxL])
sens = vec(sum(abs.(S); dims=1)) ./ size(S, 1)

# parâmetros calibrados (local, ignorado pelo git): as figuras 07 e 08 são refeitas a partir deles
CSV.write(joinpath(PROC, "mimic_calibracao.csv"), DataFrame(fits))

# ------------------------------------------------------------------ relatório
mediana(v) = (w = sort(filter(!isnan, v)); isempty(w) ? NaN : w[cld(length(w), 2)])
io = IOBuffer()
println(io, "# Modelo calibrado com dados reais do MIMIC-IV\n")
println(io, "Gerado por `scripts/pipeline_mimic.jl`. Internações anonimizadas (P01…); nenhum dado individual é versionado.\n")
println(io, "## 1. Coorte (interseção das tabelas)\n")
foreach(s -> println(io, "- ", s), funil)
println(io, "\nChoque hipovolêmico/traumático pelos CID `", join(["$(c) (CID-$(v))" for (v, c) in HYPOVOLEMIC_ICD], ", "),
            "`, sem sepse. Tratamentos: noradrenalina, cristaloides (RL, NaCl 0,9%, 0,45%, 3%, cristaloide de CC/RPA) ",
            "e concentrado de hemácias (o item \"sangue total\" não tem registros no MIMIC-IV). ",
            "Só a parte `-001` do `chartevents` estava disponível.\n")
println(io, "**Janela:** [início do primeiro tratamento − $(PRE_H) h, + $(HOURS) h]. **Ajuste** nas primeiras $(FIT_H) h, ",
            "**teste** nas $(HOURS - FIT_H) h seguintes, com os tratamentos que o paciente recebeu de fato.\n")
println(io, "## 2. Calibração e validação\n")
println(io, "Parâmetros ajustados por paciente: `Rsa` (resistência), `Vb0` (volume sanguíneo no início) e `HR` (FC intrínseca), ",
            "com o barorreflexo ligado e a FC gerada pelo modelo, restritos a faixas fisiológicas (resistência ÷3 a ×3, volume 60–120% do normal, ",
            "FC intrínseca 50–120 bpm). Erro = RMSE(PAM)² + $(W_HR)·RMSE(FC)² + $(LAMBDA)·Σ(log θ/θ₀)².\n")
println(io, "| Paciente | PAM ajuste/teste (n) | Cristaloide (mL) | Hemácias (mL) | NE máx | PAM mín | RMSE ajuste | **RMSE teste** | Palpite \"fica igual\" | Modelo sem calibrar | RMSE FC teste (bpm) | Rsa/Rsa₀ | Vb0/Vb0₀ | FC intrínseca | CO médio (L/min) |")
println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in rows
    @printf(io, "| %s | %d/%d | %.0f | %.0f | %.2f | %.0f | %.1f | **%.1f** | %.1f | %.1f | %.1f | %.2f | %.2f | %.0f | %.1f |\n",
            r.lab, r.n_fit, r.n_test, r.vol_cr, r.vol_bl, r.dose, r.pam_min, r.r_fit, r.r_test, r.r_naive,
            r.r_lit, r.r_hr, r.Rsa, r.Vb0, r.HR0, r.co)
end
n_win = count(r -> r.r_test < r.r_naive, rows)
@printf(io, "\n**Medianas (mmHg):** ajuste %.1f; teste %.1f; palpite \"fica igual\" %.1f; modelo sem calibrar %.1f. ",
        mediana([r.r_fit for r in rows]), mediana([r.r_test for r in rows]),
        mediana([r.r_naive for r in rows]), mediana([r.r_lit for r in rows]))
@printf(io, "O modelo calibrado erra menos que o palpite \"fica igual\" em **%d de %d** pacientes.\n", n_win, length(rows))
sens === nothing || @printf(io, "\n**Identificabilidade (P01):** sensibilidade média |S| da PAM: Rsa %.3f, Vb0 %.3f, FC %.3f, EmaxL %.3f.\n",
                            sens[1], sens[2], sens[3], sens[4])
println(io, "\nO \"palpite fica igual\" repete a última PAM medida no ajuste: um modelo só é útil para prever se errar menos que isso. ",
            "Figura: `docs/figuras/07_calibracao_mimic.png` (linha tracejada = fim do ajuste).\n")
println(io, "## 3. Previsão: qual tratamento para cada paciente\n")
println(io, "Com os parâmetros calibrados, cada tratamento é simulado **sozinho**, começando 1 h após o início da janela: ",
            join(["$(tr.nome)" * (tr.vol > 0 ? " (em $(Int(tr.dur ÷ 60)) min)" : "") for tr in TREATMENTS[2:end]], ", "),
            " µg/kg/min (contínua). Valores ao fim da janela (12 h). ",
            "**Critério:** maior DO₂ entre os tratamentos que mantêm PAM ≥ 65 mmHg; se nenhum mantém, maior PAM.\n")
println(io, "| Paciente | " * join(["$(tr.nome) PAM / DO₂" for tr in TREATMENTS], " | ") * " | Melhor |")
println(io, "|---|" * repeat("---|", length(TREATMENTS)) * "---|")
for (r, pr) in zip(rows, preds)
    println(io, "| ", r.lab, " | ", join([@sprintf("%.0f / %.0f", q.PAM, q.DO2) for q in pr], " | "),
            " | ", TREATMENTS[best_treatment(pr)].nome, " |")
end
cont = Dict(tr.nome => 0 for tr in TREATMENTS)
foreach(pr -> cont[TREATMENTS[best_treatment(pr)].nome] += 1, preds)
println(io, "\n**Contagem do melhor tratamento:** ", join(["$(tr.nome) $(cont[tr.nome])" for tr in TREATMENTS], "; "), ".\n")
println(io, "Figura: `docs/figuras/08_previsao_tratamentos.png`. **Atenção:** a previsão do tratamento que o paciente ",
            "não recebeu não pode ser conferida com dados (nunca vemos o mesmo paciente com outro tratamento). ",
            "O que sustenta a confiança nela é a validação da seção 2.")
write(joinpath(ROOT, "docs", "resultados_mimic.md"), String(take!(io)))
# medianas calibradas (razão em relação à literatura): usadas por scripts/prever_paciente.jl
CSV.write(joinpath(ROOT, "docs", "parametros_mimic.csv"),
          DataFrame(parametro=["Rsa", "Vb0", "HR"], n=length(rows),
                    razao_mediana=[mediana([r.Rsa for r in rows]), mediana([r.Vb0 for r in rows]),
                                   mediana([r.HR0 for r in rows]) / p_lit.HR]))
println("\n[OK] docs/resultados_mimic.md, docs/parametros_mimic.csv, figuras 07 e 08")

# figuras 07 e 08 (podem ser refeitas sem recalibrar: julia --project=. scripts/figuras_mimic.jl)
run(`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__, "figuras_mimic.jl"))`)
