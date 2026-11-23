# src/data_io.jl
# Entrada de dados: do MIMIC-IV (tabelas brutas) a um `PatientData` em grade regular, e daí a
# um `Scenario` (funções do tempo) que alimenta a EDO. Também define a coorte (CID de choque
# hipovolêmico) e os grupos de tratamento.
#
# Formato processado (um CSV por internação em data/processed/, tempo em minutos):
#   t_min, PAM, HR, Hct, rate_rl_mlh, rate_blood_mlh, dose
# `rate_rl_mlh` = todos os cristaloides (RL, NaCl 0,9% e outros), `rate_blood_mlh` = hemácias;
# ambos em mL/h (valem do instante da linha até a próxima), `dose` em µg/kg/min.

using Dates: DateTime, DateFormat, Millisecond
using CodecZlib: GzipDecompressorStream

"""
IDs de item do MIMIC-IV usados pelo projeto, conferidos em `d_items`/`d_labitems` da versão
instalada. `whole_blood` (221013) existe no dicionário mas tem **zero** registros no
`inputevents` ("Not In Use"): o hemoderivado com volume é o concentrado de hemácias (`prbc`).
"""
const MIMIC_ITEMIDS = (
    chart = (HR=220045, MAP_invasive=220052, MAP_nibp=220181, SBP=220050, DBP=220051),
    input = (ringer=225828, saline=225158, saline_half=225159, saline_hyper=225161,
             crystalloid_or=226364, crystalloid_pacu=226375,
             prbc=225168, prc_old=220996, whole_blood=221013, norepinephrine=221906),
    lab   = (hct=51221, hgb=51222),
)

"Cristaloides: Ringer lactato, NaCl 0,9%, 0,45%, 3% e cristaloide do centro cirúrgico/RPA."
const CRYSTALLOID_ITEMS = (225828, 225158, 225159, 225161, 226364, 226375)
"Hemoderivados com volume: concentrado de hemácias (atual e antigo) e sangue total."
const BLOOD_ITEMS = (225168, 220996, 221013)

"""
Códigos de diagnóstico `(versão, prefixo)` da coorte. Choque hipovolêmico: R57.1 (CID-10);
choque traumático: T79.4 (CID-10) e 958.4 (CID-9). A CID-9 não tem código próprio de choque
hipovolêmico: ele fica em 785.59 ("outro choque sem menção de trauma").
"""
const HYPOVOLEMIC_ICD = ((10, "R571"), (10, "T794"), (9, "78559"), (9, "9584"))
"Sepse e choque séptico (excluídos da coorte): A40/A41, R65.2x, 038, 995.91/995.92, 785.52."
const SEPSIS_ICD = ((10, "A40"), (10, "A41"), (10, "R652"),
                    (9, "038"), (9, "99591"), (9, "99592"), (9, "78552"))

"""
    hadm_with_icd(diagnoses, codes)

Conjunto de `hadm_id` com algum diagnóstico (`icd_version`, `icd_code`) que comece por um dos
prefixos de `codes`. O diagnóstico é da alta hospitalar: diz *quem*, não *quando*.
"""
function hadm_with_icd(diagnoses, codes)
    hit(r) = any(r.icd_version == v && startswith(string(r.icd_code), c) for (v, c) in codes)
    return Set(r.hadm_id for r in eachrow(diagnoses) if hit(r))
end

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
    scenario_from_patient(pd; hct_blood=0.60, baro=true, force_hr=false)

Converte os dados em `Scenario`: infusão total = cristaloide + hemácias (mL/s), hematócrito do
fluido ponderado pela mistura e dose de vasopressor. `hct_blood = 0,60` é o Hct de um
concentrado de hemácias (0,55–0,65, CLAUDE.md §7.1), o único hemoderivado com volume
registrado no MIMIC-IV.

Por padrão a FC **não** é imposta: o barorreflexo do modelo a gera, e a FC medida serve para
calibrar (é o mesmo modelo que depois prevê outros tratamentos, em que não há FC medida).
`force_hr=true` impõe a FC dos dados (interpolada).
"""
function scenario_from_patient(pd::PatientData; hct_blood=0.60, baro=true, force_hr=false)
    rl = step_function(pd.t, pd.rl_mlh ./ 3600)
    bl = step_function(pd.t, pd.blood_mlh ./ 3600)
    return Scenario(Qinf = t -> rl(t) + bl(t),
                    hct_inf = t -> (q = rl(t) + bl(t); q > 0 ? hct_blood * bl(t) / q : 0.0),
                    dose = step_function(pd.t, pd.dose),
                    hr_forced = force_hr ? linear_function(pd.t, pd.HR) : nothing,
                    baro = baro, tstops = collect(pd.t))
end

"""
    truncate_patient(pd, tmax)

Parte inicial de `pd` (instantes `t ≤ tmax`, em s): a janela de ajuste quando a calibração é
validada no restante.
"""
function truncate_patient(pd::PatientData, tmax)
    k = findall(<=(tmax), pd.t)
    return PatientData(pd.id, pd.t[k], pd.PAM[k], pd.HR[k], pd.Hct[k], pd.rl_mlh[k],
                       pd.blood_mlh[k], pd.dose[k])
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

"""
    stream_mimic_table(path, types, keep; verbose=false)

Lê um `.csv` ou `.csv.gz` em **blocos de 16 MB** (memória constante, sem descompactar o arquivo
nem carregá-lo inteiro) e guarda só as linhas em que `keep(row)` é verdadeiro.

- `types`: `Vector{Pair{Symbol,DataType}}` com as colunas a manter e seus tipos
  (`Int`, `Float64` ou `DateTime`, no formato `yyyy-mm-dd HH:MM:SS` do MIMIC).
- `keep(row)`: recebe um `NamedTuple` só com as colunas `Int` (ex.: `stay_id`, `itemid`), lidas
  direto dos bytes; as demais colunas só são convertidas para as linhas aprovadas.
- Linhas com valor ausente ou ilegível em alguma coluna escolhida são descartadas.

Pensado para o `chartevents` e o `labevents` do MIMIC-IV (vários GB).
"""
function stream_mimic_table(path, types, keep; verbose=false)
    cols = first.(types); T = last.(types)
    io = endswith(path, ".gz") ? GzipDecompressorStream(open(path)) : open(path)
    try
        header = split(strip(readline(io)), ',')
        idx = [findfirst(==(string(c)), header) for c in cols]
        any(isnothing, idx) && error("colunas ausentes em $path: ", cols[findall(isnothing, idx)])
        maxidx = maximum(idx)
        intpos = [i for i in eachindex(cols) if T[i] === Int]
        keepnames = Tuple(cols[intpos])
        out = [Any[] for _ in cols]
        starts = zeros(Int, maxidx); ends = zeros(Int, maxidx)
        buf = Vector{UInt8}(undef, 1 << 24)
        carry = UInt8[]
        done_mb = 0; next_report = 512

        function process_line!(chunk, a, b)
            b >= a && chunk[b] == 0x0d && (b -= 1)
            b < a && return
            n = _field_bounds!(starts, ends, chunk, a, b, maxidx)
            n < maxidx && return
            ivals = ntuple(k -> _parse_int(chunk, starts[idx[intpos[k]]], ends[idx[intpos[k]]]), length(intpos))
            any(isnothing, ivals) && return
            keep(NamedTuple{keepnames}(ivals)) || return
            vals = Vector{Any}(undef, length(cols))
            for i in eachindex(cols)
                s, e = starts[idx[i]], ends[idx[i]]
                s > e && return
                txt = String(chunk[s:e])
                v = T[i] === Int ? _parse_int(chunk, s, e) :
                    T[i] === Float64 ? tryparse(Float64, txt) :
                    tryparse(DateTime, txt, _MIMIC_DTFORMAT)
                v === nothing && return
                vals[i] = v
            end
            for i in eachindex(cols)
                push!(out[i], vals[i])
            end
        end

        while !eof(io)
            nb = readbytes!(io, buf)
            chunk = isempty(carry) ? buf[1:nb] : vcat(carry, view(buf, 1:nb))
            lastnl = findlast(==(0x0a), chunk)
            if lastnl === nothing
                carry = chunk
                continue
            end
            i = 1
            while i <= lastnl
                j = findnext(==(0x0a), chunk, i)
                process_line!(chunk, i, j - 1)
                i = j + 1
            end
            carry = chunk[lastnl+1:end]
            done_mb += nb >> 20
            if verbose && done_mb >= next_report
                println("    ... ", done_mb, " MB descompactados, ", length(out[1]), " linhas mantidas"); flush(stdout)
                next_report += 512
            end
        end
        isempty(carry) || process_line!(carry, 1, length(carry))
        return DataFrame([cols[i] => identity.(out[i]) for i in eachindex(cols)])
    finally
        close(io)
    end
end

const _MIMIC_DTFORMAT = DateFormat("yyyy-mm-dd HH:MM:SS")

"Preenche `starts`/`ends` com os limites dos campos de uma linha CSV (respeita aspas); devolve quantos campos leu (até `maxn`)."
function _field_bounds!(starts, ends, buf, a, b, maxn)
    n = 0; s = a; inq = false
    @inbounds for i in a:b
        c = buf[i]
        if c == 0x22
            inq = !inq
        elseif c == 0x2c && !inq
            n += 1; starts[n] = s; ends[n] = i - 1; s = i + 1
            n >= maxn && return n
        end
    end
    n += 1; starts[n] = s; ends[n] = b
    return n
end

"Inteiro (possivelmente negativo) escrito em `buf[a:b]`; `nothing` se vazio ou não numérico."
function _parse_int(buf, a, b)
    a > b && return nothing
    i = a; neg = false
    if buf[i] == 0x2d
        neg = true; i += 1
    end
    i > b && return nothing
    v = 0
    @inbounds while i <= b
        d = buf[i] - 0x30
        d > 0x09 && return nothing
        v = v * 10 + Int(d); i += 1
    end
    return neg ? -v : v
end

"Limites fisiológicos para descartar artefatos de medida (linha arterial lavada, flush etc.)."
const PHYSIO_LIMITS = (MAP=(20.0, 200.0), HR=(20.0, 250.0))

"""
    filter_physiological!(chart)

Remove de `chart` (colunas `itemid`, `valuenum`) as medidas de PAM fora de 20–200 mmHg e de FC
fora de 20–250 bpm (`PHYSIO_LIMITS`). Devolve o número de linhas removidas.
"""
function filter_physiological!(chart)
    n0 = nrow(chart)
    filter!(chart) do r
        lo, hi = r.itemid == MIMIC_ITEMIDS.chart.HR ? PHYSIO_LIMITS.HR : PHYSIO_LIMITS.MAP
        lo <= r.valuenum <= hi
    end
    return n0 - nrow(chart)
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
    # linhas "Rewritten" são versões substituídas de uma ordem: contá-las duplica o volume
    hasproperty(inp, :statusdescription) && filter!(r -> r.statusdescription != "Rewritten", inp)

    t0 = t0 === nothing ? minimum(mp.charttime) : t0
    tend = maximum(mp.charttime)
    nbin = Int(floor(_secs(tend, t0) / dt)) + 1
    hours === nothing || (nbin = min(nbin, Int(hours * 3600 ÷ dt)))
    grid = [(i - 1) * dt for i in 1:nbin]

    tmap = [_secs(t, t0) for t in mp.charttime]
    thr = [_secs(t, t0) for t in hrdf.charttime]
    tlab = [_secs(t, t0) for t in labs.charttime]
    PAM = [_bin_mean(tmap, mp.valuenum, g, g + dt) for g in grid]
    HR = [_bin_mean(thr, hrdf.valuenum, g, g + dt) for g in grid]      # NaN onde não há medida
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
    rl = volume_rate(CRYSTALLOID_ITEMS)
    bl = volume_rate(BLOOD_ITEMS)
    dose = zeros(nbin)
    for r in eachrow(filter(r -> r.itemid == ids.input.norepinephrine, inp))
        s = _secs(r.starttime, t0); e = _secs(r.endtime, t0)
        for (i, g) in enumerate(grid)
            s <= g < e && (dose[i] += r.rate)
        end
    end
    return PatientData(string(stay_id), grid, PAM, HR, Hct, rl, bl, dose)
end
