# src/treatments.jl
# Os quatro tratamentos comparados na previsão e o critério para escolher o melhor.
# Usados pelo pipeline do MIMIC e pela previsão de um paciente informado no terminal.

"""
Tratamentos comparados, cada um aplicado **sozinho**: volume [mL] infundido em `dur` [s],
hematócrito do fluido e dose de noradrenalina [µg/kg/min] contínua. Hemácias = concentrado
(Hct 0,60), já que o MIMIC-IV não registra sangue total.
"""
const TREATMENTS = [(nome="Sem tratamento",      vol=0.0,    dur=3600.0, hct=0.0,  dose=0.0),
                    (nome="Cristaloide 1000 mL", vol=1000.0, dur=3600.0, hct=0.0,  dose=0.0),
                    (nome="Hemácias 600 mL",     vol=600.0,  dur=3600.0, hct=0.60, dose=0.0),
                    (nome="Noradrenalina 0,1",   vol=0.0,    dur=3600.0, hct=0.0,  dose=0.1)]

"Cenário de um tratamento iniciado em `t_start` [s], com a FC gerada pelo barorreflexo."
function treatment_scenario(tr, t_start)
    return Scenario(Qinf=window_rate(tr.vol, t_start, tr.dur), hct_inf=tr.hct,
                    dose=t -> (t >= t_start ? tr.dose : 0.0),
                    tstops=[0.0, t_start, t_start + tr.dur])
end

"""
    best_treatment(res)

Índice do melhor tratamento em `res` (um resultado por tratamento, com campos `PAM` e `DO2`
ao fim da simulação). Critério: **maior oferta de O₂ entre os que mantêm PAM ≥ 65 mmHg**; se
nenhum mantém, a maior PAM.
"""
function best_treatment(res)
    ok = findall(r -> r.PAM >= 65, res)
    return isempty(ok) ? argmax([r.PAM for r in res]) : ok[argmax([res[i].DO2 for i in ok])]
end
