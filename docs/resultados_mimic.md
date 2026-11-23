# Modelo calibrado com dados reais do MIMIC-IV

Gerado por `scripts/pipeline_mimic.jl`. Internações anonimizadas (P01…); nenhum dado individual é versionado.

## 1. Coorte (interseção das tabelas)

- 1171 internações hospitalares com choque hipovolêmico/traumático; 937 sem sepse; 956 internações de UTI
- 936 com algum tratamento: noradrenalina 440, cristaloide 927, hemácias 547 (33 s)
- 936 com sinais vitais no arquivo de chartevents disponível
- 220 com PAM invasiva (≥ 5) e FC (≥ 3) em cada metade da janela e tratamento na janela (≥ 250 mL ou noradrenalina)
- 10 calibradas (as de menor stay_id), 5 delas com noradrenalina na janela

Choque hipovolêmico/traumático pelos CID `R571 (CID-10), T794 (CID-10), 78559 (CID-9), 9584 (CID-9)`, sem sepse. Tratamentos: noradrenalina, cristaloides (RL, NaCl 0,9%, 0,45%, 3%, cristaloide de CC/RPA) e concentrado de hemácias (o item "sangue total" não tem registros no MIMIC-IV). Só a parte `-001` do `chartevents` estava disponível.

**Janela:** [início do primeiro tratamento − 1.0 h, + 12.0 h]. **Ajuste** nas primeiras 6.0 h, **teste** nas 6.0 h seguintes, com os tratamentos que o paciente recebeu de fato.

## 2. Calibração e validação

Parâmetros ajustados por paciente: `Rsa` (resistência), `Vb0` (volume sanguíneo no início) e `HR` (FC intrínseca), com o barorreflexo ligado e a FC gerada pelo modelo, restritos a faixas fisiológicas (resistência ÷3 a ×3, volume 60–120% do normal, FC intrínseca 50–120 bpm). Erro = RMSE(PAM)² + 0.25·RMSE(FC)² + 5.0·Σ(log θ/θ₀)².

| Paciente | PAM ajuste/teste (n) | Cristaloide (mL) | Hemácias (mL) | NE máx | PAM mín | RMSE ajuste | **RMSE teste** | Palpite "fica igual" | Modelo sem calibrar | RMSE FC teste (bpm) | Rsa/Rsa₀ | Vb0/Vb0₀ | FC intrínseca | CO médio (L/min) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P01 | 28/8 | 4331 | 1750 | 0.50 | 50 | 16.4 | **16.2** | 34.5 | 22.5 | 41.0 | 0.33 | 0.82 | 93 | 10.1 |
| P02 | 11/5 | 2482 | 1000 | 0.25 | 60 | 11.6 | **12.8** | 11.5 | 30.4 | 15.1 | 0.33 | 0.60 | 56 | 8.4 |
| P03 | 9/11 | 2375 | 350 | 0.32 | 54 | 7.3 | **23.5** | 13.8 | 43.4 | 4.5 | 0.33 | 0.60 | 50 | 7.0 |
| P04 | 11/6 | 3714 | 350 | 0.00 | 65 | 9.3 | **15.2** | 8.0 | 24.3 | 21.0 | 0.92 | 0.60 | 50 | 5.3 |
| P05 | 17/7 | 6500 | 1390 | 0.00 | 73 | 12.3 | **10.9** | 11.0 | 11.9 | 12.2 | 0.84 | 1.20 | 119 | 9.8 |
| P06 | 9/7 | 340 | 700 | 0.00 | 85 | 6.6 | **13.5** | 6.6 | 14.7 | 4.3 | 2.17 | 0.60 | 96 | 3.2 |
| P07 | 8/5 | 7781 | 0 | 0.50 | 56 | 21.3 | **28.6** | 60.1 | 42.5 | 25.5 | 0.42 | 0.60 | 57 | 8.4 |
| P08 | 7/7 | 5018 | 700 | 0.12 | 48 | 26.2 | **26.3** | 8.4 | 44.3 | 23.2 | 0.33 | 0.60 | 50 | 8.1 |
| P09 | 12/6 | 4809 | 0 | 0.00 | 66 | 7.3 | **6.0** | 8.1 | 11.6 | 20.0 | 0.34 | 0.99 | 108 | 12.5 |
| P10 | 5/6 | 1551 | 275 | 0.00 | 79 | 16.7 | **26.4** | 12.9 | 11.0 | 15.2 | 2.78 | 1.20 | 118 | 4.8 |

**Medianas (mmHg):** ajuste 11.6; teste 15.2; palpite "fica igual" 11.0; modelo sem calibrar 22.5. O modelo calibrado erra menos que o palpite "fica igual" em **4 de 10** pacientes.

**Identificabilidade (P01):** sensibilidade média |S| da PAM: Rsa 0.276, Vb0 0.409, FC 0.107, EmaxL 0.014.

O "palpite fica igual" repete a última PAM medida no ajuste: um modelo só é útil para prever se errar menos que isso. Figura: `docs/figuras/07_calibracao_mimic.png` (linha tracejada = fim do ajuste).

## 3. Previsão: qual tratamento para cada paciente

Com os parâmetros calibrados, cada tratamento é simulado **sozinho**, começando 1 h após o início da janela: Cristaloide 1000 mL (em 60 min), Hemácias 600 mL (em 60 min), Noradrenalina 0,1 µg/kg/min (contínua). Valores ao fim da janela (12 h). **Critério:** maior DO₂ entre os tratamentos que mantêm PAM ≥ 65 mmHg; se nenhum mantém, maior PAM.

| Paciente | Sem tratamento PAM / DO₂ | Cristaloide 1000 mL PAM / DO₂ | Hemácias 600 mL PAM / DO₂ | Noradrenalina 0,1 PAM / DO₂ | Melhor |
|---|---|---|---|---|---|
| P01 | 73 / 867 | 77 / 905 | 82 / 1150 | 87 / 777 | Hemácias 600 mL |
| P02 | 54 / 422 | 61 / 455 | 71 / 708 | 78 / 411 | Hemácias 600 mL |
| P03 | 48 / 540 | 55 / 592 | 67 / 848 | 74 / 549 | Hemácias 600 mL |
| P04 | 79 / 418 | 82 / 439 | 86 / 561 | 89 / 350 | Hemácias 600 mL |
| P05 | 98 / 1010 | 99 / 1058 | 101 / 1228 | 107 / 952 | Hemácias 600 mL |
| P06 | 92 / 247 | 94 / 268 | 97 / 364 | 100 / 222 | Hemácias 600 mL |
| P07 | 55 / 625 | 62 / 685 | 73 / 911 | 79 / 603 | Hemácias 600 mL |
| P08 | 48 / 545 | 55 / 598 | 66 / 854 | 74 / 555 | Hemácias 600 mL |
| P09 | 81 / 1341 | 83 / 1391 | 87 / 1618 | 92 / 1195 | Hemácias 600 mL |
| P10 | 112 / 463 | 115 / 500 | 122 / 626 | 154 / 489 | Hemácias 600 mL |

**Contagem do melhor tratamento:** Sem tratamento 0; Cristaloide 1000 mL 0; Hemácias 600 mL 10; Noradrenalina 0,1 0.

Figura: `docs/figuras/08_previsao_tratamentos.png`. **Atenção:** a previsão do tratamento que o paciente não recebeu não pode ser conferida com dados (nunca vemos o mesmo paciente com outro tratamento). O que sustenta a confiança nela é a validação da seção 2.
