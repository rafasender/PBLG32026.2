## Paciente sintético (3 h, PAM a cada 60 s, ruído gaussiano de 2 mmHg)

Verdade: `Rsa = 1.257`, `Vb0 = 4045.0` mL, `EmaxL = 2.21`. Infusões: RL 1500 mL (30–50 min), noradrenalina 0,1 µg/kg/min (60–120 min), sangue 800 mL (90–110 min).

RMSE do modelo basal sem calibrar: **2.84 mmHg**.

### Sensibilidade local normalizada (média de |S| na janela)

| Parâmetro | média \|S\| |
|---|---|
| Rsa | 0.131 |
| Vb0 | 0.306 |
| EmaxL | 0.015 |

cond(SᵀS) = 2.16e+05  (a condição sozinha não denuncia o problema: o que aponta o parâmetro não identificável é a sensibilidade média muito menor, aqui `EmaxL`).

### Resultados da calibração

| Parâmetros ajustados | λ | RMSE (mmHg) | estimativa | verdade | erro relativo |
|---|---|---|---|---|---|
| Rsa, Vb0 | 0.0 | 1.85 | 1.25 / 4.01e+03 | 1.26 / 4.05e+03 | -0.9% / -1.0% |
| Rsa, Vb0, EmaxL | 0.0 | 1.85 | 1.1 / 4.01e+03 / 6.67e+03 | 1.26 / 4.05e+03 / 2.21 | -12.7% / -0.9% / +301734.4% |
| Rsa, Vb0, EmaxL | 0.5 | 1.86 | 1.13 / 4.16e+03 / 2.72 | 1.26 / 4.05e+03 / 2.21 | -10.3% / +2.9% / +23.2% |

