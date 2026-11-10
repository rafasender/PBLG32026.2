# scripts/pipeline_mimic.jl
# Pipeline com dados REAIS do MIMIC-IV: coorte -> leitura em streaming -> janelas por paciente
# -> CSV processado -> modelo -> calibração (Rsa, Vb0) por paciente.
#
#   julia --project=. scripts/pipeline_mimic.jl
#
# Arquivos esperados (NÃO versionados, ver .gitignore) em data/raw/mimic/:
#   icustays.csv.gz, inputevents.csv.gz, labevents.csv.gz
# e o chartevents em data/raw/mimic/chartevents.csv.gz ou data/Dados Mimic/chartevents.csv-001.gz.
#
# Coorte (proxy de choque com reposição; o MIMIC não marca "hemorragia"):
#   internações com noradrenalina que, na janela [início da NE − 1 h, + 12 h], têm PAM arterial
#   invasiva (≥ 12 medidas, ~1 por hora), FC (≥ 6) e ≥ 250 mL de RL/hemácias/sangue total.
# Resultados agregados vão para docs/resultados_mimic.md (internações anônimas: P01, P02, ...).

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))

using Cardio0D
using CSV, DataFrames, Dates, Printf, Plots

gr()
const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW  = joinpath(ROOT, "data", "raw", "mimic")
const PROC = joinpath(ROOT, "data", "processed")
const FIGDIR = joinpath(ROOT, "docs", "figuras")
mkpath(PROC); mkpath(FIGDIR)

# ------------------------------------------------------------------ configuração
const PRE_H, HOURS = 1.0, 12.0      # janela: [t_NE − PRE_H, t_NE − PRE_H + HOURS]
# O MIMIC-IV registra sinais vitais ~de hora em hora: exigir ≥ 1 PAM invasiva por hora de janela.
const MIN_MAP, MIN_HR, MIN_VOL = 12, 6, 250.0
const NMAX = 10                      # nº de pacientes calibrados
const DT = 60.0

chart_path = first(filter(isfile, [joinpath(RAW, "chartevents.csv.gz"),
                                   joinpath(ROOT, "data", "Dados Mimic", "chartevents.csv-001.gz")]))
ids = MIMIC_ITEMIDS
tempo(f) = (t = time(); r = f(); (r, time() - t))

# ------------------------------------------------------------------ 1. coorte candidata
println("1/5 inputevents + icustays ...")
inp, s = tempo() do
    d = CSV.read(joinpath(RAW, "inputevents.csv.gz"), DataFrame;
                 select=[:stay_id, :starttime, :endtime, :itemid, :amount, :rate, :statusdescription],
                 types=Dict(:starttime => DateTime, :endtime => DateTime))
    filter!(r -> r.statusdescription != "Rewritten", d)
    return d
end
@printf("    %d linhas, %d internações (%.0f s)\n", nrow(inp), length(unique(inp.stay_id)), s)
icu = CSV.read(joinpath(RAW, "icustays.csv.gz"), DataFrame; select=[:subject_id, :hadm_id, :stay_id])

ne = filter(r -> r.itemid == ids.input.norepinephrine, inp)
t_ne = combine(groupby(ne, :stay_id), :starttime => minimum => :t_ne)
cand = Set(t_ne.stay_id)
@printf("    %d internações com noradrenalina\n", length(cand))

# ------------------------------------------------------------------ 2. chartevents (streaming)
cache = joinpath(PROC, "mimic_chart_coorte.csv")
chart_all = if isfile(cache)
    println("2/5 chartevents: usando cache ", relpath(cache, ROOT))
    CSV.read(cache, DataFrame; types=Dict(:charttime => DateTime))
else
    println("2/5 chartevents em streaming (2–3 GB comprimidos, alguns minutos) ...")
    vit = (ids.chart.HR, ids.chart.MAP_invasive, ids.chart.MAP_nibp)
    d, s = tempo() do
        stream_mimic_table(chart_path,
            [:stay_id => Int, :charttime => DateTime, :itemid => Int, :valuenum => Float64],
            r -> r.itemid in vit && r.stay_id in cand; verbose=true)
    end
    @printf("    %d linhas de %d internações (%.0f s)\n", nrow(d), length(unique(d.stay_id)), s)
    CSV.write(cache, d)
    d
end
n_rem = filter_physiological!(chart_all)
println("    artefatos removidos (limites fisiológicos): ", n_rem)
println("    internações do arquivo de chartevents com NE: ", length(unique(chart_all.stay_id)),
        "  (o arquivo é só uma parte do chartevents completo)")

# ------------------------------------------------------------------ 3. seleção das janelas
println("3/5 seleção de janelas ...")
therapy = (ids.input.ringer, ids.input.prbc, ids.input.prc_old, ids.input.whole_blood)
sel = NamedTuple[]
gb_chart = groupby(chart_all, :stay_id)
gb_inp = groupby(inp, :stay_id)
tne_of = Dict(zip(t_ne.stay_id, t_ne.t_ne))
for stay in sort(unique(chart_all.stay_id))
    tne = tne_of[stay]
    t0 = tne - Minute(round(Int, PRE_H * 60)); t1 = t0 + Minute(round(Int, HOURS * 60))
    ch = filter(r -> t0 <= r.charttime <= t1, gb_chart[(stay_id=stay,)])
    nmap = count(==(ids.chart.MAP_invasive), ch.itemid)
    nhr = count(==(ids.chart.HR), ch.itemid)
    (nmap >= MIN_MAP && nhr >= MIN_HR) || continue
    vol = sum(skipmissing(filter(r -> r.itemid in therapy && t0 <= r.starttime <= t1,
                                 gb_inp[(stay_id=stay,)]).amount); init=0.0)
    vol >= MIN_VOL || continue
    push!(sel, (; stay, t0, t1, nmap, nhr, vol))
    length(sel) >= NMAX && break
end
isempty(sel) && error("nenhuma internação atende aos critérios da coorte")
@printf("    %d internações selecionadas (as de menor stay_id que cumprem os critérios)\n", length(sel))

# ------------------------------------------------------------------ 4. labevents (só Hct, só as internações escolhidas)
println("4/5 labevents (Hct) em streaming ...")
hadm_of = Dict(r.stay_id => r.hadm_id for r in eachrow(icu))
hadms = Set(hadm_of[x.stay] for x in sel)
labs_all, s = tempo() do
    stream_mimic_table(joinpath(RAW, "labevents.csv.gz"),
        [:hadm_id => Int, :charttime => DateTime, :itemid => Int, :valuenum => Float64],
        r -> r.itemid == ids.lab.hct && r.hadm_id in hadms; verbose=true)
end
@printf("    %d medidas de Hct (%.0f s)\n", nrow(labs_all), s)

# ------------------------------------------------------------------ 5. pacientes + calibração
println("5/5 pacientes e calibração ...")
p0 = baseline_params()
rows = NamedTuple[]
panels = Any[]
for (k, x) in enumerate(sel)
    lab = @sprintf("P%02d", k)
    ch = filter(r -> x.t0 <= r.charttime <= x.t1, gb_chart[(stay_id=x.stay,)])
    ip = DataFrame(gb_inp[(stay_id=x.stay,)])
    lb = filter(r -> r.hadm_id == hadm_of[x.stay] && x.t0 - Hour(6) <= r.charttime <= x.t1, labs_all)
    pd = build_patient(x.stay, ch, ip, lb; dt=DT, t0=x.t0, hours=HOURS)
    write_patient_csv(joinpath(PROC, "mimic_$(lab).csv"), pd)     # id interno só no nome local, ignorado pelo git
    pam = filter(!isnan, pd.PAM)
    r_base = try rmse(simulate_patient(p0, pd), pd.PAM) catch; NaN end
    fit = try calibrate(p0, pd, [:Rsa, :Vb0]; maxiters=120) catch e; (println("    ", lab, " falhou: ", e); nothing) end
    vol_rl = sum(pd.rl_mlh) * DT / 3600; vol_bl = sum(pd.blood_mlh) * DT / 3600
    push!(rows, (; lab, n=length(pam), vol_rl, vol_bl, dose=maximum(pd.dose),
                 pam_ini=pam[1], pam_min=minimum(pam), r_base,
                 r_fit=fit === nothing ? NaN : fit.rmse,
                 Rsa=fit === nothing ? NaN : fit.θ[1] / p0.Rsa,
                 Vb0=fit === nothing ? NaN : fit.θ[2] / p0.Vb0))
    @printf("    %s: %d PAM, RL %.0f mL, sangue %.0f mL, NE máx %.2f | RMSE %.1f -> %.1f mmHg\n",
            lab, length(pam), vol_rl, vol_bl, maximum(pd.dose), r_base, rows[end].r_fit)
    if fit !== nothing
        y0 = simulate_patient(p0, pd); y1 = simulate_patient(fit.p, pd)
        pl = plot(pd.t ./ 60, pd.PAM, seriestype=:scatter, ms=2, color=:gray, label="PAM medida",
                  title=lab, xlabel="min", ylabel="PAM (mmHg)", legend=:bottomright, titlefontsize=9)
        plot!(pl, pd.t ./ 60, y0, color=:orange, ls=:dash, label="basal")
        plot!(pl, pd.t ./ 60, y1, color=:blue, lw=2, label="calibrado")
        push!(panels, pl)
    end
end

isempty(panels) || savefig(plot(panels...; layout=(ceil(Int, length(panels) / 2), 2),
                                size=(1100, 260 * ceil(Int, length(panels) / 2))),
                           joinpath(FIGDIR, "07_calibracao_mimic.png"))

ok = filter(r -> !isnan(r.r_fit), rows)
io = IOBuffer()
println(io, "# Pipeline com dados reais do MIMIC-IV\n")
println(io, "Gerado por `scripts/pipeline_mimic.jl`. Internações anonimizadas (P01…); nenhum dado individual é versionado.\n")
println(io, "**Janela:** [início da noradrenalina − $(PRE_H) h, + $(HOURS) h], PAM arterial invasiva (dado do MIMIC ~horário) ",
            "numa grade de $(Int(DT)) s (NaN onde não há medida). **Critérios:** ≥ $MIN_MAP medidas de PAM, ≥ $MIN_HR de FC, ≥ $(Int(MIN_VOL)) mL ",
            "de RL/hemácias/sangue total. Proxy de choque com reposição (o MIMIC não marca hemorragia).\n")
println(io, "| Paciente | PAM (n) | RL (mL) | Sangue (mL) | NE máx (µg/kg/min) | PAM inicial | PAM mín | RMSE basal | RMSE calibrado | Rsa/Rsa₀ | Vb0/Vb0₀ |")
println(io, "|---|---|---|---|---|---|---|---|---|---|---|")
for r in rows
    @printf(io, "| %s | %d | %.0f | %.0f | %.2f | %.0f | %.0f | %.1f | %.1f | %.2f | %.2f |\n",
            r.lab, r.n, r.vol_rl, r.vol_bl, r.dose, r.pam_ini, r.pam_min, r.r_base, r.r_fit, r.Rsa, r.Vb0)
end
if !isempty(ok)
    @printf(io, "\n**Mediana do RMSE:** basal %.1f mmHg → calibrado %.1f mmHg (%d pacientes calibrados).\n",
            sort([r.r_base for r in ok])[cld(length(ok), 2)], sort([r.r_fit for r in ok])[cld(length(ok), 2)], length(ok))
end
println(io, "\nO RMSE calibrado é o erro **na própria janela de ajuste** (sem validação em dados independentes); ",
            "`Rsa/Rsa₀` e `Vb0/Vb0₀` são razões em relação ao paciente basal. Figura: `docs/figuras/07_calibracao_mimic.png`.")
write(joinpath(ROOT, "docs", "resultados_mimic.md"), String(take!(io)))
println("\n[OK] docs/resultados_mimic.md e docs/figuras/07_calibracao_mimic.png")
