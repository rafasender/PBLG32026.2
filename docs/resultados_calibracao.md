## Paciente sintético (3 h, PAM a cada 60 s, ruído gaussiano de 2 mmHg)

Verdade: `Rsa = 1.257`, `Vb0 = 3965.0` mL, `EmaxL = 2.21`. Infusões: RL 1500 mL (30–50 min), noradrenalina 0,1 µg/kg/min (60–120 min), sangue 800 mL (90–110 min).

RMSE do modelo basal sem calibrar: **3.94 mmHg**.

### Sensibilidade local normalizada (média de |S| na janela)

| Parâmetro | média \|S\| |
|---|---|
| Rsa | 0.115 |
| Vb0 | 0.370 |
| EmaxL | 0.015 |

cond(SᵀS) = 3.68e+05  (a condição sozinha não denuncia o problema: o que aponta o parâmetro não identificável é a sensibilidade média muito menor, aqui `EmaxL`).

### Resultados da calibração

| Parâmetros ajustados | λ | RMSE (mmHg) | estimativa | verdade | erro relativo |
|---|---|---|---|---|---|
| Rsa, Vb0 | 0.0 | 1.85 | 1.25 / 3.93e+03 | 1.26 / 3.97e+03 | -0.8% / -0.8% |
| Rsa, Vb0, EmaxL | 0.0 | 1.85 | 1.08 / 3.93e+03 / 2.93e+03 | 1.26 / 3.97e+03 / 2.21 | -14.4% / -1.0% / +132426.1% |
| Rsa, Vb0, EmaxL | 0.5 | 1.86 | 1.13 / 4.03e+03 / 2.73 | 1.26 / 3.97e+03 / 2.21 | -9.8% / +1.6% / +23.5% |

