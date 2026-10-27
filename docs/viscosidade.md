# Relação hematócrito–viscosidade μ(Hct) (entregável 1.5)

Equações e fontes: PDF da Frente B (atividade B.2) — Pries, Neuhaus & Gaehtgens (1992, *Am J
Physiol* 263:H1770) para tubos de vidro e Pries et al. (1994, *Circ Res* 75:904) para microvasos
*in vivo*. Implementação: `blood_viscosity_pries(hct; D, invivo)` em
[src/Cardio0D.jl](../src/Cardio0D.jl); os valores abaixo estão travados em `test/runtests.jl`.

## Relação adotada no modelo

**Opção A da B.2:** `μ = μ_plasma·exp(2,5·Hct)` (`blood_viscosity`). O modelo só usa a razão
`R(t) = R₀·μ(Hct(t))/μ(Hct₀)` (E18), então `μ_plasma` se cancela. Em vaso grande essa razão fica a
até 3 % da de Pries para Hct entre 0,15 e 0,55. A **opção B** (Pries *in vivo*, D = 50–100 µm) fica
disponível para a análise de sensibilidade (tarefa 3.4).

## μ_rel = μ/μ_plasma

| Hct | exp(2,5·Hct) | Pries in vitro, D = 1000 µm | Pries in vivo, D = 100 µm | Pries in vivo, D = 50 µm |
|---|---|---|---|---|
| 0,20 | 1,65 | 1,70 | 1,50 | 1,44 |
| 0,30 | 2,12 | 2,18 | 1,83 | 1,71 |
| 0,40 | 2,72 | 2,80 | 2,26 | 2,06 |
| 0,50 | 3,49 | 3,64 | 2,85 | 2,54 |

Razão μ(Hct 0,25)/μ(Hct 0,40) (hemodiluição): 0,687 (exponencial e Pries em vaso grande), 0,733
(D = 100 µm) e 0,742 (D = 30 µm).

## Krieger–Dougherty

`blood_viscosity_krieger` com os valores padrão (`hct_max = 0,62`, `[η] = 2,5`) dá ≈ 6 cP em Hct 0,40
e ≈ 35 cP em 0,55, contra ≈ 3,4 e ≈ 5,0 cP de Pries: **não usar** sem reajustar os parâmetros.

## Limitações (B.2, §3.4)

Pries usa o hematócrito de *descarga*; o `labevents` do MIMIC-IV dá o hematócrito sistêmico, e os
dois diferem em microvasos (efeito Fåhræus). A viscosidade também depende do cisalhamento.
