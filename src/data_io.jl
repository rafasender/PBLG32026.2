# src/data_io.jl
# Entrada de dados: do MIMIC-IV (tabelas brutas) a um `PatientData` em grade regular, e daí a
# um `Scenario` (funções do tempo) que alimenta a EDO. Também gera pacientes sintéticos para
# testar todo o fluxo antes do acesso à base.
#
# Formato processado (um CSV por internação em data/processed/, tempo em minutos):
#   t_min, PAM, HR, Hct, rate_rl_mlh, rate_blood_mlh, dose
# `rate_*` em mL/h (valem do instante da linha até a próxima), `dose` em µg/kg/min.

using Dates: DateTime, Millisecond

"""
IDs de item do MIMIC-IV usados pelo projeto. **Conferir em `d_items`/`d_labitems` da versão
instalada** (os códigos mudam entre versões e hospitais).
"""
const MIMIC_ITEMIDS = (
    chart = (HR=220045, MAP_invasive=220052, MAP_nibp=220181, SBP=220050, DBP=220051),
    input = (ringer=225828, saline=225158, prbc=225168, norepinephrine=221906),
    lab   = (hct=51221, hgb=51222),
)

"""
    PatientData

Série de um paciente em grade regular (tempo em segundos):
`PAM` observada [mmHg] (NaN = sem medida), `HR` [bpm], `Hct` (fração, NaN = sem dado),
taxas de infusão `rl_mlh`/`blood_mlh` [mL/h] e `dose` de vasopressor [µg/kg/min].
"""
struct PatientData
    id::String
    t::Vector{Float64}
    PAM::Vector{Float64}
    HR::Vector{Float64}
    Hct::Vector{Float64}
    rl_mlh::Vector{Float64}
    blood_mlh::Vector{Float64}
    dose::Vector{Float64}
end

# ----------------------------------------------------------- séries -> funções do tempo
"Função degrau contínua à direita: `ys[i]` vale de `ts[i]` até `ts[i+1]` (0 antes do início)."
function step_function(ts, ys)
    return t -> (i = searchsortedlast(ts, t); i < 1 ? 0.0 : float(ys[i]))
end

"Interpolação linear com extrapolação constante (ignora NaN, que são removidos)."
function linear_function(ts, ys)
    ok = findall(!isnan, ys)
    tt = ts[ok]; yy = ys[ok]
    return function (t)
        isempty(tt) && return NaN
        t <= tt[1] && return yy[1]
        t >= tt[end] && return yy[end]
        i = searchsortedlast(tt, t)
        w = (t - tt[i]) / (tt[i+1] - tt[i])
        return yy[i] * (1 - w) + yy[i+1] * w
    end
end

"""
    scenario_from_patient(pd; hct_blood=0.40, baro=true)

Converte os dados em `Scenario`: infusão total = RL + sangue (mL/s), hematócrito do fluido
ponderado pela mistura, dose de vasopressor e **FC imposta** (`hr_forced`) com os dados.
"""
function scenario_from_patient(pd::PatientData; hct_blood=0.40, baro=true)
    rl = step_function(pd.t, pd.rl_mlh ./ 3600)
    bl = step_function(pd.t, pd.blood_mlh ./ 3600)
    hr = linear_function(pd.t, pd.HR)
    return Scenario(Qinf = t -> rl(t) + bl(t),
                    hct_inf = t -> (q = rl(t) + bl(t); q > 0 ? hct_blood * bl(t) / q : 0.0),
                    dose = step_function(pd.t, pd.dose),
                    hr_forced = hr, baro = baro, tstops = collect(pd.t))
end

# ------------------------------------------------------------------------ CSV processado
"Grava o paciente no formato processado (tempo em minutos)."
function write_patient_csv(path, pd::PatientData)
    df = DataFrame(t_min = pd.t ./ 60, PAM = pd.PAM, HR = pd.HR, Hct = pd.Hct,
                   rate_rl_mlh = pd.rl_mlh, rate_blood_mlh = pd.blood_mlh, dose = pd.dose)
    CSV.write(path, df)
    return path
end

"Lê um CSV processado (colunas de [`write_patient_csv`](@ref)); células vazias viram NaN."
function load_patient(path; id=splitext(basename(path))[1])
    df = CSV.read(path, DataFrame; missingstring=["", "NA", "NaN"])
    nan(col) = [ismissing(x) ? NaN : float(x) for x in df[!, col]]
    return PatientData(String(id), nan(:t_min) .* 60, nan(:PAM), nan(:HR), nan(:Hct),
                       nan(:rate_rl_mlh), nan(:rate_blood_mlh), nan(:dose))
end

# ------------------------------------------------------------------------- MIMIC bruto
"""
    read_mimic_table(path, cols; itemids=nothing, stay_ids=nothing)

Lê **só as colunas necessárias** de uma tabela do MIMIC (`.csv` ou `.csv.gz`) e filtra por
`itemid`/`stay_id`. O `chartevents` tem vários GB: nunca leia todas as colunas.
"""
function read_mimic_table(path, cols; itemids=nothing, stay_ids=nothing)
    df = DataFrame(CSV.File(path; select=cols))
    itemids === nothing || filter!(r -> r.itemid in itemids, df)
    stay_ids === nothing || filter!(r -> r.stay_id in stay_ids, df)
    return df
end

_secs(t, t0) = Dates_value(t - t0) / 1000
Dates_value(x) = Int(Millisecond(x).value)

"Média da coluna `val` das linhas de `tt` em `[a, b)`; NaN se não houver medida."
function _bin_mean(tt, val, a, b)
    idx = findall(t -> a <= t < b, tt)
    return isempty(idx) ? NaN : sum(val[idx]) / length(idx)
end

"Preenche NaN com o último valor observado (forward-fill); NaN iniciais ficam NaN."
function _ffill(x)
    y = copy(x); last = NaN
    for i in eachindex(y)
        isnan(y[i]) ? (y[i] = last) : (last = y[i])
    end
    return y
end

"""
    build_patient(stay_id, chart, inputs, labs; dt=60.0, t0=nothing, hours=nothing)

Monta o [`PatientData`](@ref) de uma internação a partir de tabelas no esquema do MIMIC-IV:
- `chart`: `stay_id, charttime, itemid, valuenum` (PAM, FC);
- `inputs`: `stay_id, starttime, endtime, itemid, amount, rate` (RL/soro, PRBC, noradrenalina);
- `labs`: `charttime, itemid, valuenum` **já restrito à internação** (o `labevents` não tem
  `stay_id`: filtre por `hadm_id` e janela da UTI antes).

Tempo zero = `t0` (padrão: primeira medida de PAM). Volumes são distribuídos na grade de `dt`
segundos preservando o total (mL/h). PAM invasiva tem prioridade; sem ela usa a não invasiva.
"""
function build_patient(stay_id, chart, inputs, labs; dt=60.0, t0=nothing, hours=nothing)
    ids = MIMIC_ITEMIDS
    ch = filter(r -> r.stay_id == stay_id, chart)
    mapi = filter(r -> r.itemid == ids.chart.MAP_invasive, ch)
    mapn = filter(r -> r.itemid == ids.chart.MAP_nibp, ch)
    mp = nrow(mapi) > 0 ? mapi : mapn
    nrow(mp) == 0 && error("stay_id $stay_id sem PAM")
    hrdf = filter(r -> r.itemid == ids.chart.HR, ch)
    inp = filter(r -> r.stay_id == stay_id, inputs)

    t0 = t0 === nothing ? minimum(mp.charttime) : t0
    tend = maximum(mp.charttime)
    nbin = Int(floor(_secs(tend, t0) / dt)) + 1
    hours === nothing || (nbin = min(nbin, Int(hours * 3600 ÷ dt)))
    grid = [(i - 1) * dt for i in 1:nbin]

    tmap = [_secs(t, t0) for t in mp.charttime]
    thr = [_secs(t, t0) for t in hrdf.charttime]
    tlab = [_secs(t, t0) for t in labs.charttime]
    PAM = [_bin_mean(tmap, mp.valuenum, g, g + dt) for g in grid]
    HR = _ffill([_bin_mean(thr, hrdf.valuenum, g, g + dt) for g in grid])
    Hct = _ffill([_bin_mean(tlab[findall(==(ids.lab.hct), labs.itemid)],
                            labs.valuenum[findall(==(ids.lab.hct), labs.itemid)], g, g + dt) / 100
                  for g in grid])

    # volumes (mL) distribuídos por sobreposição com cada intervalo da grade
    function volume_rate(itemids)
        rate = zeros(nbin)
        for r in eachrow(filter(r -> r.itemid in itemids, inp))
            s = _secs(r.starttime, t0); e = _secs(r.endtime, t0)
            e <= s && continue
            for (i, g) in enumerate(grid)
                ov = max(0.0, min(e, g + dt) - max(s, g))
                rate[i] += r.amount * ov / (e - s) / dt * 3600      # mL/h
            end
        end
        return rate
    end
    rl = volume_rate((ids.input.ringer, ids.input.saline))
    bl = volume_rate((ids.input.prbc,))
    dose = zeros(nbin)
    for r in eachrow(filter(r -> r.itemid == ids.input.norepinephrine, inp))
        s = _secs(r.starttime, t0); e = _secs(r.endtime, t0)
        for (i, g) in enumerate(grid)
            s <= g < e && (dose[i] += r.rate)
        end
    end
    return PatientData(string(stay_id), grid, PAM, HR, Hct, rl, bl, dose)
end

# ---------------------------------------------------------------------- paciente sintético
"""
    synthetic_patient(p_true; id, hours=3, dt=60.0, noise=2.0, seed=1, rl=..., blood=..., dose=...)

Gera dados com a **verdade conhecida** para testar a calibração: simula o modelo lento com
`p_true` (volume inicial deficitário via `Vb0`), com infusões de RL/sangue e noradrenalina em
degraus (funções `(t_s) -> valor`), e amostra a PAM com ruído gaussiano (`noise` mmHg).
Devolve `(PatientData, sol)`.
"""
function synthetic_patient(p_true; id="sintetico", hours=3.0, dt=60.0, noise=2.0, seed=1,
                           rl=t -> (1800 <= t < 3000 ? 1500.0 : 0.0),
                           blood=t -> (5400 <= t < 6600 ? 800.0 : 0.0),
                           dose=t -> (3600 <= t < 7200 ? 0.1 : 0.0), HR0=95.0)
    grid = collect(0.0:dt:hours*3600)
    pd0 = PatientData(id, grid, fill(NaN, length(grid)), fill(HR0, length(grid)),
                      fill(NaN, length(grid)), rl.(grid), blood.(grid), dose.(grid))
    sc = scenario_from_patient(pd0)
    u0 = settled_state(p_true, sc)
    sol = simulate_mean(p_true, sc, (grid[1], grid[end]); u0, saveat=grid, reltol=1e-6, abstol=1e-6)
    h = [mean_hemodynamics(p_true, sol(t), sc, t) for t in grid]
    state = Random_LCG(seed)
    PAM = [x.Psa + noise * _randn!(state) for x in h]
    Hct = [x.e.Hct for x in h]
    pd = PatientData(id, grid, PAM, pd0.HR, Hct, pd0.rl_mlh, pd0.blood_mlh, pd0.dose)
    return pd, sol
end

# gerador pseudoaleatório mínimo (evita depender do pacote Random no módulo; reprodutível)
mutable struct Random_LCG
    s::UInt64
end
Random_LCG(seed::Integer) = Random_LCG(UInt64(seed) * 6364136223846793005 + 1442695040888963407)
function _rand!(r::Random_LCG)
    r.s = r.s * 6364136223846793005 + 1442695040888963407
    return ((r.s >> 11) + 1) / (2.0^53 + 2)
end
# Box–Muller
_randn!(r::Random_LCG) = sqrt(-2 * log(_rand!(r))) * cos(2π * _rand!(r))
