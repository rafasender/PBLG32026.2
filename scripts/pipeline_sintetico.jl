# scripts/pipeline_sintetico.jl
# Entregável 1.6 — pipeline de dados com dados SINTÉTICOS (CSV.jl + DataFrames.jl).
#
#   tabelas brutas no esquema do MIMIC-IV  ->  leitura (só as colunas necessárias)
#   -> filtragem (internação, itemid, valores fisiológicos) -> série em grade regular
#   -> CSV processado -> releitura.
#
# As tabelas sintéticas imitam `chartevents`, `inputevents` e `labevents` (colunas extras,
# itens irrelevantes e artefatos de medida incluídos de propósito, para o filtro ter o que
# filtrar). Quando o acesso ao MIMIC-IV sair, basta apontar `RAW` para os CSVs reais.
#
#   julia --project=. scripts/pipeline_sintetico.jl

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))

using Cardio0D
using CSV, DataFrames
using Dates: DateTime, Second
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW  = joinpath(ROOT, "data", "raw", "sintetico")
const PROC = joinpath(ROOT, "data", "processed")
mkpath(RAW); mkpath(PROC)

ids = MIMIC_ITEMIDS
t0 = DateTime(2150, 3, 1, 8, 0, 0)
at(s) = t0 + Second(round(Int, s))

# ------------------------------------------------------------------ 1. tabelas brutas
# Verdade: o paciente sintético do modelo (hipovolemia + RL + sangue + noradrenalina).
p0 = baseline_params()
ptrue = with_params(p0; Rsa=1.25 * p0.Rsa, Vb0=0.86 * p0.Vb_ref, EmaxL=0.85 * p0.EmaxL)
truth, _ = synthetic_patient(ptrue; id="1001", hours=3.0, dt=60.0, noise=2.0, seed=3)

STAY, OTHER = 1001, 2002     # a internação 2002 é "de outro paciente": deve ser descartada

# chartevents: PAM invasiva (a cada 1 min) e FC (a cada 5 min) + um item irrelevante (SpO2)
chart_rows = NamedTuple[]
for (i, t) in enumerate(truth.t)
    push!(chart_rows, (subject_id=1, hadm_id=11, stay_id=STAY, charttime=at(t),
                       itemid=ids.chart.MAP_invasive, valuenum=truth.PAM[i], warning=0))
    i % 5 == 1 && push!(chart_rows, (subject_id=1, hadm_id=11, stay_id=STAY, charttime=at(t),
                       itemid=ids.chart.HR, valuenum=truth.HR[i], warning=0))
    i % 5 == 1 && push!(chart_rows, (subject_id=1, hadm_id=11, stay_id=STAY, charttime=at(t),
                       itemid=220277, valuenum=97.0, warning=0))             # SpO2: irrelevante
    push!(chart_rows, (subject_id=2, hadm_id=22, stay_id=OTHER, charttime=at(t),
                       itemid=ids.chart.MAP_invasive, valuenum=55.0, warning=0))  # outro paciente
end
# artefatos: PAM = 0 (linha arterial lavada) e 300 (flush); devem ser removidos
push!(chart_rows, (subject_id=1, hadm_id=11, stay_id=STAY, charttime=at(600),
                   itemid=ids.chart.MAP_invasive, valuenum=0.0, warning=1))
push!(chart_rows, (subject_id=1, hadm_id=11, stay_id=STAY, charttime=at(1200),
                   itemid=ids.chart.MAP_invasive, valuenum=300.0, warning=1))
CSV.write(joinpath(RAW, "chartevents.csv"), DataFrame(chart_rows))

# inputevents: converte as taxas (mL/h) do paciente sintético em eventos [start, end)
input_rows = NamedTuple[]
function add_events!(rows, vec, itemid; is_rate=false)
    i = 1
    while i <= length(vec)
        if vec[i] > 0
            j = i
            while j + 1 <= length(vec) && vec[j+1] == vec[i]; j += 1; end
            dur_h = (truth.t[j] + 60 - truth.t[i]) / 3600
            push!(rows, (subject_id=1, hadm_id=11, stay_id=STAY, starttime=at(truth.t[i]),
                         endtime=at(truth.t[j] + 60), itemid=itemid,
                         amount=is_rate ? 0.0 : vec[i] * dur_h, rate=is_rate ? vec[i] : 0.0))
            i = j + 1
        else
            i += 1
        end
    end
end
add_events!(input_rows, truth.rl_mlh, ids.input.ringer)
add_events!(input_rows, truth.blood_mlh, ids.input.prbc)
add_events!(input_rows, truth.dose, ids.input.norepinephrine; is_rate=true)
push!(input_rows, (subject_id=1, hadm_id=11, stay_id=STAY, starttime=at(0), endtime=at(600),
                   itemid=225799, amount=50.0, rate=0.0))                    # item irrelevante
CSV.write(joinpath(RAW, "inputevents.csv"), DataFrame(input_rows))

# labevents: sem stay_id (como no MIMIC); Hct em % e hemoglobina como item irrelevante
lab_rows = [(subject_id=1, hadm_id=11, charttime=at(t), itemid=ids.lab.hct,
             valuenum=100 * truth.Hct[i]) for (i, t) in enumerate(truth.t) if i % 30 == 1]
push!(lab_rows, (subject_id=1, hadm_id=11, charttime=at(0), itemid=ids.lab.hgb, valuenum=13.5))
CSV.write(joinpath(RAW, "labevents.csv"), DataFrame(lab_rows))

# ------------------------------------------------------------------ 2. leitura + filtragem
chart = read_mimic_table(joinpath(RAW, "chartevents.csv"),
                         [:stay_id, :charttime, :itemid, :valuenum];
                         itemids=[ids.chart.MAP_invasive, ids.chart.MAP_nibp, ids.chart.HR],
                         stay_ids=[STAY])
inputs = read_mimic_table(joinpath(RAW, "inputevents.csv"),
                          [:stay_id, :starttime, :endtime, :itemid, :amount, :rate];
                          itemids=collect(values(ids.input)), stay_ids=[STAY])
labs = read_mimic_table(joinpath(RAW, "labevents.csv"), [:hadm_id, :charttime, :itemid, :valuenum];
                        itemids=[ids.lab.hct])
filter!(r -> r.hadm_id == 11, labs)       # labevents não tem stay_id: filtra por hadm_id

n_antes = nrow(chart)
# limites fisiológicos (descarta artefatos)
filter!(r -> r.itemid == ids.chart.HR ? 20 <= r.valuenum <= 250 : 20 <= r.valuenum <= 200, chart)
n_removidos = n_antes - nrow(chart)

# ------------------------------------------------------------------ 3. série regular + CSV
pd = build_patient(STAY, chart, inputs, labs; dt=60.0)
out_csv = joinpath(PROC, "sintetico_pipeline.csv")
write_patient_csv(out_csv, pd)
pd2 = load_patient(out_csv)

# ------------------------------------------------------------------ 4. conferências
ok = Bool[]
check(name, cond) = (push!(ok, cond); println(cond ? "  [OK]    " : "  [FALHA] ", name))

println("\nPipeline sintético — conferências")
check("só a internação $STAY e itens de interesse (colunas lidas: $(join(names(chart), ", ")))",
      all(chart.stay_id .== STAY) && all(in.(chart.itemid, Ref([ids.chart.MAP_invasive, ids.chart.HR]))))
check("artefatos removidos ($n_removidos linhas fora dos limites fisiológicos)", n_removidos == 2)
check("série em grade regular de 60 s ($(length(pd.t)) pontos)",
      all(diff(pd.t) .== 60.0) && length(pd.t) == length(truth.t))
okpam = filter(!isnan, pd.PAM)
check("PAM dentro da faixa fisiológica ($(round(minimum(okpam), digits=1))–$(round(maximum(okpam), digits=1)) mmHg)",
      all(20 .<= okpam .<= 200))
check("PAM idêntica à verdade sintética (ruído incluído)",
      all(isapprox.(pd.PAM[.!isnan.(pd.PAM)], truth.PAM[.!isnan.(pd.PAM)]; atol=1e-9)))
vol(x) = sum(x) * 60 / 3600
check("volume de RL conservado ($(round(vol(pd.rl_mlh), digits=1)) mL)",
      isapprox(vol(pd.rl_mlh), vol(truth.rl_mlh); rtol=1e-6))
check("volume de sangue conservado ($(round(vol(pd.blood_mlh), digits=1)) mL)",
      isapprox(vol(pd.blood_mlh), vol(truth.blood_mlh); rtol=1e-6))
check("dose de noradrenalina recuperada (máx $(maximum(pd.dose)) µg/kg/min)",
      maximum(pd.dose) == maximum(truth.dose))
check("Hct em fração (0–1) e sem NaN após forward-fill",
      all(0.1 .< filter(!isnan, pd.Hct) .< 0.7) && count(isnan, pd.Hct) <= 1)
check("CSV processado relido igual ao original",
      pd2.t == pd.t && isequal(pd2.PAM, pd.PAM) && pd2.dose == pd.dose)
sc = scenario_from_patient(pd2)
check("Scenario construído a partir do CSV (Qinf(t) em mL/s ≥ 0)",
      sc.Qinf(2000.0) >= 0 && sc.hr_forced(0.0) > 0)

# ---- pipeline -> modelo: simula com o CSV processado e compara a PAM -------------------
pam_true = simulate_patient(ptrue, pd2)      # parâmetros verdadeiros: deve cair no nível do ruído
pam_base = simulate_patient(p0, pd2)         # paciente "médio": erro maior (não sabe o déficit)
r_true, r_base = rmse(pam_true, pd2.PAM), rmse(pam_base, pd2.PAM)
println(@sprintf("\nPAM simulada × CSV: RMSE com parâmetros verdadeiros = %.2f mmHg (ruído 2,0); com baseline = %.2f mmHg",
                 r_true, r_base))
check("modelo alimentado pelo CSV reproduz a PAM no nível do ruído (RMSE < 3 mmHg)", r_true < 3.0)
# só a ordem é verificada: a margem é pequena (a PAM é pouco sensível ao déficit com os
# parâmetros da B.2), o que é um resultado a reportar, não um limite a ajustar
check("o modelo distingue os pacientes (RMSE baseline > RMSE verdadeiro)", r_base > r_true)

println("\nCSV processado: ", relpath(out_csv, ROOT))
println(first(CSV.read(out_csv, DataFrame), 5))
all(ok) ? println("\n[SUCESSO] Pipeline sintético OK.") : (println("\n[FALHA] ver acima."); exit(1))
