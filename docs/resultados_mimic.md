# Pipeline com dados reais do MIMIC-IV

Gerado por `scripts/pipeline_mimic.jl`. Internações anonimizadas (P01…); nenhum dado individual é versionado.

**Janela:** [início da noradrenalina − 1.0 h, + 12.0 h], PAM arterial invasiva (dado do MIMIC ~horário) numa grade de 60 s (NaN onde não há medida). **Critérios:** ≥ 12 medidas de PAM, ≥ 6 de FC, ≥ 250 mL de RL/hemácias/sangue total. Proxy de choque com reposição (o MIMIC não marca hemorragia).

| Paciente | PAM (n) | RL (mL) | Sangue (mL) | NE máx (µg/kg/min) | PAM inicial | PAM mín | RMSE basal | RMSE calibrado | Rsa/Rsa₀ | Vb0/Vb0₀ |
|---|---|---|---|---|---|---|---|---|---|---|
| P01 | 25 | 1515 | 0 | 0.08 | 84 | 60 | 20.8 | 9.3 | 0.20 | 1.35 |
| P02 | 14 | 7466 | 0 | 0.20 | 73 | 53 | 45.4 | 12.2 | 0.06 | 2.05 |
| P03 | 12 | 576 | 181 | 0.10 | 78 | 61 | 26.6 | 12.9 | 0.16 | 1.28 |
| P04 | 13 | 2180 | 0 | 0.10 | 71 | 63 | 35.9 | 7.4 | 0.07 | 2.06 |
| P05 | 18 | 2658 | 700 | 0.05 | 118 | 47 | 36.2 | 19.3 | 0.13 | 1.37 |
| P06 | 32 | 5329 | 1050 | 0.16 | 57 | 23 | 75.1 | 11.9 | 0.01 | 2.03 |
| P07 | 37 | 2042 | 700 | 0.03 | 81 | 55 | 25.8 | 10.1 | 0.17 | 1.26 |
| P08 | 24 | 1274 | 275 | 0.10 | 78 | 50 | 25.5 | 9.1 | 0.25 | 1.00 |
| P09 | 24 | 2012 | 350 | 0.15 | 94 | 69 | 16.5 | 10.4 | 0.28 | 1.24 |
| P10 | 12 | 504 | 283 | 0.05 | 88 | 59 | 24.7 | 6.4 | 0.18 | 1.22 |

**Mediana do RMSE:** basal 25.8 mmHg → calibrado 10.1 mmHg (10 pacientes calibrados).

O RMSE calibrado é o erro **na própria janela de ajuste** (sem validação em dados independentes); `Rsa/Rsa₀` e `Vb0/Vb0₀` são razões em relação ao paciente basal. Figura: `docs/figuras/07_calibracao_mimic.png`.
