# CLAUDE.md — Modelo 0D do sistema cardiovascular (PBL 2026.2 · G3)

> **Projeto:** "Modelagem 0D do sistema cardiovascular para simulação de choque hipovolêmico
> e comparação de terapias de reposição volêmica em UTI" — disciplina *Modelagem e Simulação
> de Sistemas Biomédicos* (Prof. Federico Aletti, Einstein, 6º sem.).
>
> Este documento é o **guia matemático** do modelo. Ele explica *o que* modelar, *por que*
> cada equação existe e *como* levá-la para Julia. Ambiente/instalação estão no `README.md` e
> ficam fora deste escopo.

---

## Status da implementação

As etapas 1–8 abaixo **estão implementadas** em `src/` e verificadas por `test/runtests.jl`
(85 testes); detalhes, resultados e limitações em [docs/RELATORIO.md](docs/RELATORIO.md).
Divergências em relação ao texto original deste guia:

| Item | Guia original | Implementado |
|---|---|---|
| Parâmetros basais (§12) | tabela fixa | **projetados por alvos** (`baseline_params`): `R = ΔP/Q`, `C_D`/`C_S` de EDV/ESV |
| Solver do pulsátil (§4.5) | `Rodas5P` obrigatório | `Rodas5P` (funciona; ~3 s/40 s); `Tsit5` também converge |
| Fase cardíaca (§5) | `mod(t,T)` | estado `φ` (`dφ/dt = FC/60`, usa `φ − floor(φ)`): `mod` dá NaN na derivada |
| Colapso venoso | não previsto | softplus em `Psv` (hemorragia grande) |
| Cristaloide (§7.2) | níveis 1 e 2 | só o **nível 1** (plasma filtrado; refil saturável) |
| Pós-carga do modelo lento | `P_a` | `P_a + dPes` (pressão de fim de sístole > PAM) |
| Calibração (§10) | LBFGS/adjunto | Nelder–Mead em `log θ` + sensibilidade por diferenças finitas |
| ModelingToolkit (§13.3) | alternativa | **não implementado** (próximo passo: migrar usando este código como gabarito) |
| `EmaxL` na calibração | candidato | **não identificável só com a PAM** (ver relatório §4.5) |

---

## 0. Como ler este documento

O modelo é construído **de forma incremental**, exatamente como no abstract do grupo e na
sequência do *Módulo Cardiovascular* da disciplina:

| Etapa | O que se adiciona | Conceito da aula | Entregável |
|---|---|---|---|
| 1 | Windkessel arterial de 2 elementos (RC) | Caso #2–#3 | A.2 |
| 2 | Árvore venosa + volumes + pressão média de enchimento | Caso #4 | A.2 |
| 3 | Coração como capacitor de capacitância variável (PV loop) | Caso #5 | A.2 |
| 4 | Circulação pulmonar + malha fechada (coração direito/esquerdo) | Casos #6–#8 | A.2 |
| 5 | Choque hipovolêmico (hemorragia) | — | B.3 |
| 6 | Terapias: Ringer Lactato, sangue total, vasopressor; Hct → μ → R | — | B.3 |
| 7 | Barorreflexo (compensação simpática) | Aula 3–4 (controle) | B.3 / A.3 |
| 8 | Calibração com MIMIC-IV, sensibilidade e validação | Aula 1 + Hester 2020 | C.3 |

Regra de ouro: **cada etapa precisa reproduzir analiticamente a anterior** antes de se adicionar a
próxima (ver §9, Verificação). Isso é o que dá credibilidade ao modelo (Pruett, Clemmer & Hester,
*Annu. Rev. Biomed. Eng.* 2020).

### Convenções

- **Unidades:** pressão em mmHg, volume em mL, tempo em s, fluxo em mL/s.
  Resistência em mmHg·s/mL, complacência em mL/mmHg, elastância em mmHg/mL.
  (1 L/min = 16,67 mL/s. Para comparar com o MIMIC, converter na saída e nunca dentro das EDOs.)
- **Analogia eletromecânica** (Aula 1 — variáveis *across* e *through*):

  | Elétrico | Hidráulico | Símbolo |
  |---|---|---|
  | Tensão V | Pressão P | *across* |
  | Corrente I | Fluxo Q | *through* |
  | Carga q | Volume V | ∫Q dt |
  | Resistor R (Ohm: ΔV = R·I) | Resistência ao fluxo (ΔP = R·Q) | R |
  | Capacitor C (q = C·V) | Complacência (ΔV = C·ΔP) | C = 1/E |
  | Diodo | Válvula cardíaca | D |
  | Indutor L | Inércia do sangue | L (opcional) |

- **Siglas do curso** (mantidas no código): ABP (pressão arterial), CVP (pressão venosa central),
  RAP (pressão do átrio direito), CO (débito cardíaco), VR (retorno venoso), TPR (resistência
  periférica total), AC/VC (complacência arterial/venosa), Pms (pressão média de enchimento
  sistêmico), HR (frequência cardíaca), SV (volume sistólico).

---

## 1. Por que um modelo 0D (parâmetros concentrados)?

O curso distingue:

- **Modelos 0D (*lumped*)** assumem que P, Q e V são **uniformes dentro de cada compartimento** num
  dado instante. Resultado: um sistema de **EDOs**, em geral duas por compartimento (conservação de
  massa e de quantidade de movimento) mais uma relação algébrica P–V.
- **Modelos 1D/2D/3D** são EDPs (Navier–Stokes). Capturam propagação e reflexão de onda, mas exigem
  geometria e custo computacional muito maiores.
- A forma 1D de Navier–Stokes **discretizada no espaço** (e sem o termo convectivo) leva de volta a um
  modelo 0D. Por isso o 0D é a escala certa para a pergunta do grupo: "como volume, resistência e
  viscosidade alteram PAM e débito ao longo de minutos a horas de UTI?"

**Limitação a declarar no relatório:** o 0D descreve **valores médios por compartimento**. Ele não
reproduz o aumento da pressão de pulso na periferia nem o pico secundário/entalhe dicrótico
(reflexão de onda). Compare sempre a **PAM** e não a morfologia da onda radial do MIMIC.

É um modelo **grey-box** (Aula 1): a estrutura vem da física (hypothesis-driven) e alguns
parâmetros são estimados dos dados (data-driven). Isso define a estratégia de calibração da §10.

---

## 2. Etapa 1 — Windkessel de 2 elementos (árvore arterial)

### 2.1 Circuito e equação de estado

Fonte de corrente Q_in(t) (ejeção do VE) alimentando um capacitor AC em paralelo com um resistor
TPR. A variável de estado é a pressão arterial P_a(t) = ABP.

Pela conservação de massa no nó arterial (Kirchhoff das correntes):

```
fluxo que entra = fluxo que é armazenado + fluxo que sai
Q_in(t)        =  AC · dP_a/dt         + (P_a − P_v)/TPR
```

logo

```
dP_a/dt = [ Q_in(t) − (P_a − P_v)/TPR ] / AC                      (E1)
```

`P_v` é a pressão a jusante (CVP, ~4 mmHg). Na versão do curso, P_v ≈ 0.

### 2.2 Solução na diástole e a constante de tempo τ

Na diástole a válvula aórtica está fechada, então Q_in = 0. Com P_v ≈ 0:

```
P_a(t) = P_a0 · exp(−t / τ),   τ = TPR · AC                        (E2)
```

- **τ é o parâmetro fundamental do Windkessel** (≈ 1–2 s). Ele aparece separado de TPR e AC só
  quando se conhece o fluxo.
- **Identificação prática:** ajuste log(P_a) × t na fase diastólica de uma onda arterial. A
  inclinação vale −1/τ. Com CO conhecido, TPR = (PAM − CVP)/CO e AC = τ/TPR.
- **Débito pelo Windkessel (Caso #2.3):** em regime médio, `CO = (ABP − CVP)/TPR` (lei de Ohm
  hidráulica). É o mesmo princípio dos métodos de *pulse contour* de monitores de UTI.

### 2.3 Validade (Caso #3)

- O que reproduz bem: **o decaimento exponencial diastólico** e a relação PAM–CO–TPR.
- Onde falha: forma sistólica, entalhe dicrótico e diferenças entre pressão central e periférica,
  porque o modelo não tem propagação de onda.
- Extensões possíveis: o **Windkessel de 3 elementos** (resistência característica Z_c em série
  antes do RC, que melhora a sístole) e o de **4 elementos** (+ inertância L). Para o projeto, o de 2
  elementos é suficiente, porque o alvo é a PAM.

### 2.4 Em Julia

Já está implementado em [src/Cardio0D.jl](src/Cardio0D.jl) (`windkessel2_ode!`) e em
[scripts/demo_windkessel2.jl](scripts/demo_windkessel2.jl), com Q_in = meia-senoide de área SV.
**Teste de verificação obrigatório:** o τ ajustado na diástole simulada tem de bater com
`R*C` (erro < 1%).

---

## 3. Etapa 2 — Árvore venosa, volumes e Pms (Caso #4)

### 3.1 Circuito

Adiciona-se, a jusante do Windkessel arterial, um segundo RC: a complacência venosa VC (grande) e a
resistência venosa R_v (pequena) até o átrio direito.

```
 Q=CO → [AC ∥ nó ABP] → TPR → [VC ∥ nó CVP] → R_v → RAP (átrio direito) → VR
```

Comparação artéria × veia, e a consequência para os parâmetros:

| | Artérias (+ arteríolas/microcirculação) | Veias |
|---|---|---|
| Papel dominante | **Resistência** (TPR) | **Complacência / reservatório** (VC) |
| Complacência | baixa (AC ≈ 1–2 mL/mmHg) | alta (VC ≈ 20–30 × AC) |
| Volume armazenado | ~15–20% | ~65–70% |

### 3.2 Volume não estressado (V₀) e volume estressado

A relação P–V de cada capacitor **não passa pela origem**. Um vaso guarda um volume V₀ mesmo com
pressão transmural nula (volume não estressado). Só o excedente distende a parede e gera pressão:

```
V_a = AC · ABP + V_a0                                               (E3)
V_v = VC · CVP + V_v0                                               (E4)
V_t = V_a + V_v = AC·ABP + VC·CVP + V_0,   V_0 = V_a0 + V_v0         (E5)
```

Valores de referência do curso: V_t ≈ 4000 mL (≈ 85% da volemia, circulação periférica),
V₀ ≈ 3000 mL, **volume de distensão** V_t − V₀ ≈ 1000 mL.

**Este é o conceito-chave do choque hipovolêmico.** Uma hemorragia retira volume, mas a pressão só
cai na proporção do volume **estressado**. Por sua vez, a venoconstrição simpática ou o vasopressor
"convertem" volume não estressado em estressado (reduzem V₀) e sustentam a pressão sem repor volume.

### 3.3 Quedas de pressão em regime (Lei de Ohm)

```
ABP − CVP = CO · TPR                                                (E6)
CVP − RAP = VR · R_v                                                (E7)
```

### 3.4 Pressão média de enchimento sistêmico (Pms)

Definição: é a pressão que se mediria em todo o sistema **com o coração parado e o sangue
redistribuído** (CO = 0, então ABP = CVP = RAP = Pms). Substituindo em (E5):

```
Pms = (V_t − V₀) / (AC + VC) = (V_t − V₀) / C_t                     (E8)
```

Pms ≈ 7 mmHg. É o "nível de enchimento" do sistema: depende **só de volume e complacência**, não do
coração.

### 3.5 Curva do retorno venoso (Guyton)

Com o coração batendo, a média ponderada das pressões pelas complacências continua igual a Pms
(conservação de volume):

```
Pms = (AC·ABP + VC·CVP) / (AC + VC)                                 (E9)
```

Usando (E6)+(E7) com CO = VR em regime, ABP = RAP + VR(TPR + R_v) e CVP = RAP + VR·R_v. Substituindo
em (E9) e isolando VR:

```
             Pms − RAP
VR = ───────────────────────────                                    (E10)
      R_v + TPR · AC/(AC + VC)
```

O denominador é a **resistência ao retorno venoso** R_VR. Leitura gráfica no plano (RAP, VR):

- **Intercepto em RAP:** VR = 0 quando RAP = Pms.
- **Intercepto em VR:** VR(0) = Pms / R_VR.
- **Inclinação:** −1/R_VR.
- Na prática, para RAP < ~0 mmHg as veias colapsam e VR satura num platô.

Interpretação clínica (Casos #4.8–4.9), que vira **hipótese testável** no modelo:

| Intervenção | Parâmetro | Efeito na curva de VR |
|---|---|---|
| Hemorragia | ↓ V_t | ↓ Pms → curva desloca **paralelamente para baixo/esquerda** |
| Volume (RL, sangue) | ↑ V_t | ↑ Pms → desloca para cima/direita |
| Venoconstrição / noradrenalina | ↓ V₀ (e ↓ VC) | ↑ Pms |
| Vasoconstrição arteriolar | ↑ TPR | ↑ R_VR → **menor inclinação** (rotação) |
| ↑ viscosidade (↑ Hct) | ↑ TPR, ↑ R_v | ↑ R_VR → rotação |
| Hemodiluição (cristaloide) | ↓ μ → ↓ R | ↓ R_VR → maior inclinação |

Repare que o termo `TPR·AC/(AC+VC)` é pequeno porque AC ≪ VC. Por isso **R_v domina R_VR**, e esse é
um resultado contraintuitivo que vale discutir.

---

## 4. Etapa 3 — O coração como capacitor variável (Caso #5)

### 4.1 Lei de Frank–Starling em forma de circuito

O ventrículo é um capacitor cuja **elastância E = 1/C muda com o tempo**: é baixa na diástole
(relaxado, enche fácil) e alta na sístole (contraído). As duas retas limite do plano P–V são:

```
Diástole:  V_D = C_D · P + V_d      (ED = 1/C_D, pequena)           (E11)
Sístole:   V_S = C_S · P + V_d      (ES = 1/C_S, grande)            (E12)
```

V_d é o volume morto (hipótese do curso: o mesmo na sístole e na diástole).

### 4.2 PV loop do ventrículo isolado (pré-carga P_f e pós-carga P_a constantes)

| Trecho | Fase | Válvulas | O que acontece |
|---|---|---|---|
| 1 → 2 | Contração isovolumétrica | mitral e aórtica fechadas | V = EDV, P sobe de P_f até P_a |
| 2 → 3 | Ejeção | aórtica aberta | P = P_a, V cai até ESV |
| 3 → 4 | Relaxamento isovolumétrico | ambas fechadas | V = ESV, P cai até P_f |
| 4 → 1 | Enchimento | mitral aberta | P = P_f, V sobe até EDV |

```
EDV = C_D·P_f + V_d,   ESV = C_S·P_a + V_d
SV  = EDV − ESV = C_D·P_f − C_S·P_a                                 (E13)
VO  = HR · (C_D·P_f − C_S·P_a)                                      (E14)
```

### 4.3 Faixa de funcionamento (curva de débito ventricular)

```
       ⎧ 0                              se P_f ≤ (C_S/C_D)·P_a      (não abre a aórtica)
VO  =  ⎨ HR·(C_D·P_f − C_S·P_a)         se (C_S/C_D)·P_a < P_f < EDV_max/C_D
       ⎩ HR·(EDV_max − C_S·P_a)         se P_f ≥ EDV_max/C_D        (limite pericárdico)
```

A **variável independente é a pré-carga** P_f (= RAP para o coração direito). Os parâmetros são HR,
C_D (complacência diastólica, relaxamento), C_S (contratilidade: ↓C_S ⇔ ↑E_S ⇔ inotropismo
positivo) e a pós-carga P_a.

### 4.4 Do modelo de dois estados para o modelo contínuo (o que vamos simular)

Para gerar **curvas de pressão, fluxo e PV loop no tempo** (entregável A.2), troca-se a comutação
abrupta C_D ↔ C_S por uma **elastância variante no tempo**:

```
P_v(t) = E(t) · (V(t) − V_d)                                        (E15)
E(t)   = E_min + (E_max − E_min) · e(t mod T)                        (E16)
```

com T = 60/HR e uma função de ativação normalizada e ∈ [0,1], por exemplo a cosseno por partes:

```
         ⎧ ½·[1 − cos(π·τ/T_s1)]                0 ≤ τ < T_s1          (contração)
e(τ) =   ⎨ ½·[1 + cos(π·(τ − T_s1)/T_s2)]       T_s1 ≤ τ < T_s1+T_s2  (relaxamento)
         ⎩ 0                                    resto do ciclo        (diástole)
```

e T_s1 ≈ 0,3·√T, T_s2 ≈ 0,5·T_s1 (a sístole encurta pouco com a taquicardia, então a diástole é
quem perde tempo de enchimento, como no paciente do Cenário A, com FC de 124 bpm).

**Ligação com a aula:** E_min = 1/C_D e E_max = 1/C_S. Ao simular o ventrículo com pré e pós-carga
fixas, o SV calculado deve tender a (E13). Esse é o **teste de verificação** desta etapa.

### 4.5 Válvulas (diodos)

```
Q_valv = max(P_montante − P_jusante, 0) / R_valv                     (E17)
```

O `max(·,0)` não é diferenciável e atrapalha a otimização e o SciMLSensitivity. Use uma versão
suave:

```
Q_valv ≈ (ΔP/R_valv) · σ(k·ΔP),   σ(x) = 1/(1+e^{−x}),  k ≈ 5–20 mmHg⁻¹
```

ou `softplus`. Com R_valv ≈ 0,005–0,01 mmHg·s/mL o sistema fica **rígido** (stiff): use um solver
implícito (`Rodas5P`, `QNDF`, `FBDF`) ou `AutoTsit5(Rodas5P())`.

---

## 5. Etapa 4 — Malha fechada: coração direito + pulmão + coração esquerdo (Casos #6–#8)

### 5.1 Topologia (6 compartimentos de volume)

```
      ┌──────────────── circulação sistêmica ───────────────┐
      │                                                     │
  [VE] ─D_ao→ [SA: C_sa] ─R_sa(=TPR)→ [SV: C_sv] ─R_sv→ D_tri → [VD]
   ↑                                                           │
 D_mi                                                        D_pul
   │                                                           ↓
  [PV: C_pv] ←──R_pv── [PA: C_pa] ←──────────────────────── R_pa
      └──────────────── circulação pulmonar ────────────────┘
```

Compartimentos: VE (ventrículo esquerdo), SA (artérias sistêmicas), SV (veias sistêmicas), VD
(ventrículo direito), PA (artérias pulmonares), PV (veias pulmonares + átrio esquerdo).
Átrios são incorporados às veias na primeira versão (refino opcional: átrios com elastância própria).

### 5.2 Escolha das variáveis de estado: volumes

Use **volumes** como estados (não pressões). Os motivos:

1. A conservação de massa fica explícita: `Σ dV_i/dt = Q_infusão − Q_hemorragia`.
2. Hemorragia e infusão entram diretamente como termos de fonte/sumidouro em mL/s.
3. V₀ (volume não estressado) aparece naturalmente, e é ele que o vasopressor altera.

Relações constitutivas (algébricas):

```
Vasos:      P_i = (V_i − V_0,i) / C_i                (+ P_th para vasos intratorácicos)
Ventrículos: P_VE = E_VE(t)·(V_VE − V_d,VE),  P_VD = E_VD(t)·(V_VD − V_d,VD)
```

Fluxos:

```
Q_ao  = valv(P_VE − P_SA)        Q_sys = (P_SA − P_SV)/R_sa
Q_tri = valv(P_SV − P_VD)        Q_pul = valv(P_VD − P_PA)
Q_lung= (P_PA − P_PV)/R_pa       Q_mi  = valv(P_PV − P_VE)
```

EDOs, uma por compartimento, cada uma igual a "entra − sai":

```
dV_VE/dt = Q_mi  − Q_ao
dV_SA/dt = Q_ao  − Q_sys
dV_SV/dt = Q_sys − Q_tri + Q_inf(t) − Q_hem(t)       ← terapias e hemorragia entram aqui
dV_VD/dt = Q_tri − Q_pul
dV_PA/dt = Q_pul − Q_lung
dV_PV/dt = Q_lung − Q_mi
```

> Hemorragia de trauma sai fisicamente do leito arterial. Colocá-la em SA ou SV muda pouco o regime
> médio, mas documente a escolha.

### 5.3 Pressão intratorácica (Caso #7)

Câmaras e vasos dentro do tórax veem a pressão **transmural**. O curso usa um *offset*
P_th ≈ −5 mmHg (tórax fechado). Em ventilação mecânica com PEEP (Cenário A: PEEP 5 cmH₂O ≈
3,7 mmHg) P_th fica menos negativa ou positiva, e isso **reduz o retorno venoso**:

```
VR = [(Pms − P_th) − (RAP − P_th)] / R_VR   (forma do curso, Caso #8)
```

Implementação: `P_i = (V_i − V_0,i)/C_i + P_th` para VE, VD, PA e PV. É um parâmetro útil para
discutir o efeito hemodinâmico da intubação vista no Cenário A (queda da PAM de 64 para 59 mmHg).

### 5.4 Ponto de funcionamento (Caso #8)

Em regime, **CO = VR**. O ponto de operação é a interseção da curva de função cardíaca (E14,
com P_f = RAP) com a curva de retorno venoso (E10):

```
HR·(C_D·RAP − C_S·P_a) = (Pms − RAP)/R_VR     →  resolver para RAP*, CO*
```

**Duas formas de obter, e as duas devem concordar (verificação cruzada):**

1. **Analítica / algébrica:** resolver o sistema não linear com `NonlinearSolve.jl` (ou à mão no caso
   linear).
2. **Numérica:** simular a malha fechada até o regime periódico e tirar as médias por ciclo de RAP e
   CO.

Gráfico-chave do relatório: as curvas de Guyton (VR e CO × RAP) **antes da hemorragia, depois da
hemorragia e após cada terapia**, com o deslocamento do ponto de funcionamento.

### 5.5 Modelo "lento" (médio por batimento) para escalas de horas

O modelo pulsátil exige passos de ~ms e o MIMIC tem amostras de minutos a horas. Para calibração e
cenários longos, derive um **modelo de valores médios**, em que o coração vira a fonte algébrica da
aula:

```
CO_E = HR·(C_D,E·P_PV − C_S,E·P_SA)       (coração esquerdo, E14)
CO_D = HR·(C_D,D·P_SV − C_S,D·P_PA)       (coração direito)
```

Os vasos continuam como RCs, e as EDOs passam a ter só 4 estados de volume (SA, SV, PA, PV), sem
oscilação no batimento: o sistema é muito mais rápido e diferenciável.
**Verificação:** as médias por ciclo do modelo pulsátil devem coincidir com o modelo lento (erro na
PAM < 2–3 mmHg). Use o pulsátil para PV loops e formas de onda e o lento para horas de UTI e
calibração.

---

## 6. Etapa 5 — Choque hipovolêmico (hemorragia)

### 6.1 Entrada de perda

```
Q_hem(t) = taxa constante em [t_ini, t_fim]   ou   perfil por degraus
```

Cenários de referência (classes ATLS, volemia ≈ 70 mL/kg; 5 L para 70 kg):

| Classe | Perda | Volume (70 kg) | Esperado no modelo |
|---|---|---|---|
| I | < 15% | < 750 mL | PAM ~ normal (compensação) |
| II | 15–30% | 750–1500 mL | ↑FC, PAM normal/limítrofe |
| III | 30–40% | 1500–2000 mL | ↓PAM, ↓CO |
| IV | > 40% | > 2000 mL | colapso |

**Sem barorreflexo (§8) o modelo superestima a queda de pressão nas classes I–II.** Mostrar isso
(barorreflexo ON/OFF) é uma ótima figura para explicar a "compensação simpática" do abstract.

### 6.2 O que acontece matematicamente

↓V_t → ↓volume estressado → ↓Pms (E8) → curva de VR desce → ↓RAP (pré-carga) → ↓SV (E13) →
↓CO → ↓PAM (E6). A sequência causal inteira está nas equações das aulas.

---

## 7. Etapa 6 — Terapias, hematócrito e viscosidade

### 7.1 Balanço de hemácias e plasma (hematócrito dinâmico)

Separe o volume intravascular total em hemácias e plasma, com dois estados extras:

```
V_blood = Σ V_i (todos os compartimentos)
dV_rbc/dt = Hct_inf·Q_inf − Hct(t)·Q_hem
Hct(t)    = V_rbc / V_blood
```

- Hemorragia: perde sangue **com o Hct atual**, então o Hct não muda instantaneamente.
- Ringer Lactato / cristaloide: Hct_inf = 0, o que leva à **hemodiluição**.
- Sangue total: Hct_inf ≈ 0,40 (concentrado de hemácias: 0,55–0,65).

Dados do Cenário A para teste: o Hct cai de 42% → 40% → 39% com hidratação venosa.

### 7.2 Cristaloide: extravasamento para o interstício

Só ~20–30% do cristaloide permanece intravascular após ~1 h. Há dois níveis de modelagem:

**Nível 1 (recomendado para começar)**, com compartimento intersticial e troca de 1ª ordem:

```
dV_is/dt   = J_cap
J_cap      = k_f · [ (V_pl − V_pl,ref) − α·(V_is − V_is,ref) ]     (filtração líquida)
```

J_cap é subtraído de dV_SV/dt (ou distribuído entre os compartimentos). k_f (1/s) define a
constante de tempo de equilíbrio (~20–30 min).

**Nível 2 (fisiológico), Starling:**

```
J_cap = K_f · [ (P_c − P_is) − σ·(π_c − π_is) ],   π_c = π_c0 · (V_pl0 / V_pl)
```

O cristaloide dilui as proteínas plasmáticas, π_c cai e a filtração aumenta. Assim o modelo
**explica** por que o RL "vaza", e o sangue total (que traz proteína e hemácias) não. Este é o
argumento central da comparação de terapias.

### 7.3 Viscosidade e resistência (entregável B.2)

Pela lei de Poiseuille, em um tubo:

```
R = 8·μ·L / (π·r⁴)    ⇒    R(t) = R_ref · μ(Hct(t)) / μ(Hct_ref)    (E18)
```

Aplique o fator principalmente às resistências da microcirculação (TPR = R_sa e R_pa).

Relação μ(Hct), a definir com a Frente B:

| Forma | Expressão | Observação |
|---|---|---|
| Placeholder atual | μ = μ_pl·exp(2,5·Hct) | [src/Cardio0D.jl](src/Cardio0D.jl) |
| Tipo Krieger/Dougherty | μ_rel = (1 − Hct/Hct_max)^(−[η]·Hct_max) | ajustável |
| Pries et al. 1992 | μ_rel(Hct, D) in vitro / in vivo | inclui Fåhræus–Lindqvist (depende do diâmetro) |

Efeito previsto: o cristaloide **↓Hct → ↓μ → ↓TPR**, o que melhora o fluxo mas reduz a PAM para o
mesmo CO e a capacidade de transporte de O₂. O sangue total mantém μ. Quantificar esse *trade-off*
é a contribuição original do grupo.

### 7.4 Vasopressor (noradrenalina)

Use um modelo **farmacodinâmico Emax** com compartimento de efeito (atraso de ~1–2 min):

```
dC_e/dt = (Dose(t) − C_e) / τ_e
f(C_e)  = C_e^n / (EC50^n + C_e^n)
TPR  = TPR_0 · (1 + G_R · f)          (vasoconstrição arteriolar)
V_0,sv = V_0,sv0 · (1 − G_V · f)      (venoconstrição: recruta volume estressado)
E_max = E_max0 · (1 + G_E · f)        (efeito β1 inotrópico, menor)
```

Dose(t) vem do `inputevents` do MIMIC (µg/kg/min). G_R, G_V, G_E e EC50 são candidatos à
calibração. Note que o vasopressor sobe a PAM via ↑TPR **e** ↑Pms (↓V₀), sem repor volume, que é a
discussão de "vasoconstrição excessiva" do abstract.

---

## 8. Etapa 7 — Barorreflexo (compensação simpática)

É controle fisiológico em malha fechada (Aulas 3–4: teoria de sistemas e controle). Um modelo
mínimo, inspirado em Ursino (1998):

1. **Sensor (barorreceptor):** PAM filtrada, `dP̄/dt = (P_SA − P̄)/τ_p`, com τ_p ≈ 2–5 s.
2. **Sinal simpático** (sigmoide centrada no *set-point* P_n ≈ 90–93 mmHg):
   `s = 1 / (1 + exp((P̄ − P_n)/k_s))` ∈ (0,1). Com PAM baixa, s → 1.
3. **Efetores** (resposta de 1ª ordem, cada um com τ e ganho próprios):

```
dx_j/dt = ( G_j·(s − s_0) − x_j ) / τ_j,   j ∈ {HR, TPR, V_0,sv, E_max}
HR    = HR_0 + x_HR          (τ ~ 2 s, vagal/simpático)
TPR   = TPR_0 · (1 + x_R)    (τ ~ 6 s)
V_0,sv= V_0,sv0 − x_V        (τ ~ 20 s)
E_max = E_max0 · (1 + x_E)   (τ ~ 8 s)
```

s_0 é o valor de s no equilíbrio basal, de modo que os efetores começam em zero.

Use uma *flag* `baroreflex_on` e compare as curvas. **Atenção na calibração:** se HR vier medido
do MIMIC, você pode **impor HR(t)** como entrada (forçante) em vez de gerá-lo pelo reflexo. Isso
reduz parâmetros e é o que o abstract sugere ("débito cardíaco ou frequência cardíaca como input").

---

## 9. Verificação (antes de qualquer comparação com dados)

Seguindo o framework **ASME V&V 40** citado por Pruett, Clemmer & Hester (2020): **plan →
risk → verification → validation**, sempre em relação a um **contexto de uso (COU)**. Aqui, o COU
é comparar *qualitativa e semiquantitativamente* terapias de reposição em choque hipovolêmico; não é
decisão clínica individual.

Testes de **verificação de código** (colocar em [test/runtests.jl](test/runtests.jl)):

| # | Teste | Critério |
|---|---|---|
| V1 | Windkessel: τ ajustado na diástole = R·C | erro < 1% |
| V2 | Conservação de massa: Σ V_i(t) − V_total(0) − ∫(Q_inf − Q_hem) | ≈ 0 (tolerância do solver) |
| V3 | Coração parado (E constante, HR → 0): todas as pressões → Pms de (E8) | erro < 0,1 mmHg |
| V4 | Ventrículo isolado com P_f, P_a fixos: SV simulado = C_D·P_f − C_S·P_a (E13) | erro < 2% |
| V5 | Ponto de funcionamento: médias simuladas de RAP e CO = interseção de (E10) × (E14) | erro < 5% |
| V6 | Modelo lento × médias do pulsátil | PAM < 3 mmHg |
| V7 | Convergência numérica: reltol 1e-6 vs 1e-8 | diferença desprezível |
| V8 | Sanidade basal: PAM ≈ 90–95, CO ≈ 5 L/min, SV ≈ 70 mL, EF ≈ 55–65%, CVP ≈ 2–8 | faixa fisiológica |

---

## 10. Etapa 8 — Calibração, identificabilidade e sensibilidade

### 10.1 Problema de otimização

Para cada paciente *k*, sobre a janela temporal alinhada:

```
θ* = argmin_θ  J(θ),   J(θ) = (1/N)·Σ_i [ PAM_sim(t_i; θ) − PAM_med(t_i) ]²  + λ·‖(θ − θ_lit)/θ_lit‖²
```

- **Entradas impostas (forçantes)** vindas do MIMIC: HR(t), Q_inf(t) por tipo de fluido, dose de
  vasopressor(t) e Hct(t) (interpolado dos `labevents`, ou gerado pelo modelo e comparado).
- **Saídas comparadas:** PAM (principal) e, se possível, PAD.
- O termo de regularização λ evita parâmetros fisiologicamente absurdos. É uma forma "pobre" de um
  *prior* (que futuramente vira Turing.jl).
- **Otimize em log(θ)**, porque R e C são positivos e variam em ordens de grandeza. Combine uma busca
  global grosseira (ou multistart) com um refinamento local por gradiente (`LBFGS` via
  `OptimizationOptimJL`, com gradiente por `AutoForwardDiff` ou adjunto do `SciMLSensitivity`).
- Calibre sobre o **modelo lento** (§5.5). O modelo pulsátil fica só para a visualização final.

### 10.2 Identificabilidade (o ponto mais importante desta frente)

Com **só a PAM medida não é possível identificar todos os R e C**. Exemplo clássico, já visto no
próprio Windkessel: a PAM média em regime determina apenas TPR·CO, e a dinâmica diastólica só τ =
R·C. Nucci & Cobelli mostram o mesmo na mecânica respiratória: o modelo de Mead, embora mais
realista, **não é estruturalmente identificável**, e por isso modelos de 1ª ordem (RC) são usados
na prática.

Estratégia recomendada:

1. **Fixar da literatura** (tabela da Frente B): complacências, V_d, E_min, parâmetros pulmonares,
   resistências de válvulas.
2. **Estimar 2–4 parâmetros por paciente:** TPR_0, V_t (volemia efetiva) ou V₀,sv, E_max,VE e o
   ganho do vasopressor G_R.
3. Checar a **identificabilidade prática**: o perfil de J ao variar cada parâmetro (*profile
   likelihood*) deve ter um mínimo bem definido; a matriz de sensibilidade não pode ter colunas
   quase colineares (calcule `cond(SᵀS)`).

### 10.3 Análise de sensibilidade

Sensibilidade local normalizada (adimensional, comparável entre parâmetros):

```
S_ij(t) = (∂y_i/∂θ_j) · (θ_j / y_i)
```

com y = PAM, CO, Hct e θ = {μ (via Hct), taxa de infusão, dose de vasopressor, TPR, V_t, ...}. Em
Julia, use `ForwardDiff` sobre a solução (poucos parâmetros) ou `SciMLSensitivity` (adjunto, se
forem muitos). Para a discussão, complemente com uma análise **global** (variar ±20–50% em grade
ou Morris/Sobol) nos 3–4 parâmetros mais relevantes.

### 10.4 Validação

- Separar pacientes (ou janelas temporais) para **calibração** e para **validação** (Aula 1: o modelo
  identificado é validado com um *dataset* independente).
- Reportar RMSE da PAM, viés (Bland–Altman) e se o modelo prevê a **direção e a magnitude** da
  resposta a uma terapia que não entrou na calibração.

---

## 11. Métricas de comparação das terapias (saídas do modelo)

Para cada cenário (hemorragia + RL, + sangue total, + vasopressor, combinações):

| Métrica | Definição no modelo |
|---|---|
| Tempo até PAM ≥ 65 mmHg | primeiro t com média móvel da PAM ≥ 65 |
| Volume necessário | ∫Q_inf dt até atingir a meta |
| CO e SV no novo regime | médias por ciclo |
| Fração retida intravascular | ΔV_blood / V_infundido (após 30/60 min) |
| Hct final e μ final | §7.1 e §7.3 |
| Oferta de O₂ (proxy) | DO₂ ≈ CO · 1,34 · Hb · SaO₂ (Hb ≈ Hct/3) |
| Sobrecarga | V_is acumulado (edema), P_PV (risco de edema pulmonar) |

A métrica DO₂ é importante: o cristaloide pode **subir a PAM e o CO, mas baixar a DO₂** pela
hemodiluição. Essa é a nuance que justifica o projeto.

---

## 12. Parâmetros iniciais (adulto 70 kg, a substituir pela tabela da Frente B)

Valores de ordem de grandeza, coerentes com o curso (V_t ≈ 4000, V₀ ≈ 3000, Pms ≈ 7) e com modelos
clássicos (Ursino 1998; Smith et al. 2004). **Cada valor final precisa de fonte citada.**

| Parâmetro | Símbolo | Valor inicial | Unidade |
|---|---|---|---|
| Frequência cardíaca | HR | 75 | bpm |
| Volemia total | V_blood | 5000 | mL |
| Hematócrito | Hct | 0,42 | — |
| Resistência periférica total | R_sa (TPR) | 1,0 | mmHg·s/mL |
| Complacência arterial sistêmica | C_sa (AC) | 1,5 | mL/mmHg |
| Volume não estressado arterial | V_0,sa | 600 | mL |
| Resistência venosa sistêmica | R_sv (R_v) | 0,05 | mmHg·s/mL |
| Complacência venosa sistêmica | C_sv (VC) | 100 | mL/mmHg |
| Volume não estressado venoso | V_0,sv | 2500 | mL |
| Resistência arterial pulmonar | R_pa | 0,08 | mmHg·s/mL |
| Complacência arterial pulmonar | C_pa | 4 | mL/mmHg |
| Complacência venosa pulmonar | C_pv | 10 | mL/mmHg |
| V₀ pulmonares (PA + PV) | V_0,p | 150 + 350 | mL |
| Elastância máx. VE / VD | E_max | 2,5 / 0,6 | mmHg/mL |
| Elastância mín. VE / VD | E_min | 0,08 / 0,05 | mmHg/mL |
| Volume morto VE / VD | V_d | 10 / 10 | mL |
| Resistência valvar | R_valv | 0,005 | mmHg·s/mL |
| Pressão intratorácica | P_th | −4 | mmHg |

Ao final, ajuste V_0,sv para que o basal dê PAM ≈ 93 mmHg e CO ≈ 5 L/min. Esse procedimento de
"ajuste do basal" deve ser documentado.

---

## 13. Implementação em Julia: arquitetura sugerida

### 13.1 Organização do código (`src/`)

```
src/
  Cardio0D.jl          # módulo principal, exports
  params.jl            # structs de parâmetros (Base.@kwdef), conversões de unidade
  elastance.jl         # e(τ), E(t)
  blood.jl             # μ(Hct), Hct dinâmico, R(μ)
  valves.jl            # Q_valv suave
  closed_loop.jl       # EDO pulsátil (6+ estados)
  mean_model.jl        # modelo lento (§5.5)
  scenarios.jl         # Q_hem(t), Q_inf(t), Dose(t): funções/callbacks
  baroreflex.jl
  analysis.jl          # PAM, SV, EF, PV loop, médias por ciclo, Guyton
  calibration.jl       # custo, Optimization.jl
```

Mantenha o estilo atual de [src/Cardio0D.jl](src/Cardio0D.jl): `Base.@kwdef struct` para
parâmetros, funções `*_ode!(du, u, p, t)` in-place e docstrings em português com a equação.

### 13.2 Esqueleto da malha fechada (estilo do repositório)

```julia
Base.@kwdef struct CircParams{T<:Real}
    HR::T = 75.0
    Rsa::T = 1.0;   Csa::T = 1.5;  V0sa::T = 600.0
    Rsv::T = 0.05;  Csv::T = 100.0; V0sv::T = 2500.0
    Rpa::T = 0.08;  Cpa::T = 4.0;  V0pa::T = 150.0
    Cpv::T = 10.0;  V0pv::T = 350.0
    EmaxL::T = 2.5; EminL::T = 0.08; VdL::T = 10.0
    EmaxR::T = 0.6; EminR::T = 0.05; VdR::T = 10.0
    Rvalv::T = 0.005; Pth::T = -4.0; kvalv::T = 10.0
end

act(τ, T) = (Ts1 = 0.3*sqrt(T); Ts2 = 0.5*Ts1;
    τ < Ts1        ? 0.5*(1 - cos(π*τ/Ts1)) :
    τ < Ts1 + Ts2  ? 0.5*(1 + cos(π*(τ - Ts1)/Ts2)) : 0.0)

valve(ΔP, R, k) = ΔP / R / (1 + exp(-k*ΔP))      # diodo suave

function closed_loop!(du, u, p, t)
    prm, Qinf, Qhem = p                            # Qinf(t), Qhem(t) são funções
    VL, Vsa, Vsv, VR, Vpa, Vpv = u
    T = 60 / prm.HR;  e = act(mod(t, T), T)
    EL = prm.EminL + (prm.EmaxL - prm.EminL)*e
    ER = prm.EminR + (prm.EmaxR - prm.EminR)*e
    # relações constitutivas (P = (V - V0)/C [+ Pth se intratorácico])
    PL  = EL*(VL - prm.VdL) + prm.Pth
    PR  = ER*(VR - prm.VdR) + prm.Pth
    Psa = (Vsa - prm.V0sa)/prm.Csa
    Psv = (Vsv - prm.V0sv)/prm.Csv
    Ppa = (Vpa - prm.V0pa)/prm.Cpa + prm.Pth
    Ppv = (Vpv - prm.V0pv)/prm.Cpv + prm.Pth
    # fluxos
    Qao  = valve(PL - Psa,  prm.Rvalv, prm.kvalv)
    Qsys = (Psa - Psv)/prm.Rsa
    Qtri = valve(Psv - PR,  prm.Rvalv + prm.Rsv, prm.kvalv)   # R_v em série com a tricúspide
    Qpul = valve(PR - Ppa,  prm.Rvalv, prm.kvalv)
    Qlng = (Ppa - Ppv)/prm.Rpa
    Qmi  = valve(Ppv - PL,  prm.Rvalv, prm.kvalv)
    # conservação de massa
    du[1] = Qmi  - Qao
    du[2] = Qao  - Qsys
    du[3] = Qsys - Qtri + Qinf(t) - Qhem(t)
    du[4] = Qtri - Qpul
    du[5] = Qpul - Qlng
    du[6] = Qlng - Qmi
end

# prob = ODEProblem(closed_loop!, u0, (0.0, 60.0), (CircParams(), t->0.0, t->0.0))
# sol  = solve(prob, Rodas5P(); reltol=1e-6, abstol=1e-6, saveat=0.005)
```

Pontos de atenção:

- **u0:** distribua V_blood respeitando o V₀ de cada compartimento (por exemplo, comece com
  pressões "razoáveis" e converta para volumes) e simule ~20–50 ciclos até o regime periódico.
  Teste: a média por ciclo de cada V não muda mais que 0,1% entre ciclos.
- **Descontinuidades** (início/fim da hemorragia, bolus): use `tstops` ou `PresetTimeCallback` /
  `DiscreteCallback` para o solver não "pular" o degrau.
- **Parâmetros que mudam** (Hct → R, vasopressor, barorreflexo) devem virar **estados extras** do
  vetor `u` (V_rbc, V_is, C_e, P̄, x_j) e não variáveis globais, para manter a diferenciabilidade.
- Para gradientes via SciMLSensitivity/ForwardDiff, não use `Float64` fixo nos structs: o
  `T<:Real` paramétrico já permite `Dual`.

### 13.3 Alternativa com ModelingToolkit (abordagem por componentes)

O abstract propõe a biblioteca de componentes do MTK. A ideia é a mesma do circuito: cada
compartimento é um componente com um *pin* (P, Q). As conexões geram automaticamente as leis de
Kirchhoff (Σ Q = 0 nos nós e igualdade de P), e o `mtkcompile` elimina as equações algébricas.

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@connector Pin begin
    P(t)                      # pressão (across)
    Q(t), [connect = Flow]    # fluxo (through)
end

@mtkmodel Compliance begin    # vaso: V = C·P + V0  ⇔  C·dP/dt = Q
    @parameters begin C; V0; Pth = 0.0 end
    @components begin p = Pin() end
    @variables begin V(t) end
    @equations begin
        D(V) ~ p.Q
        p.P ~ (V - V0)/C + Pth
    end
end

@mtkmodel Resistor begin      # ΔP = R·Q
    @parameters begin R end
    @components begin a = Pin(); b = Pin() end
    @equations begin
        a.P - b.P ~ R*a.Q
        a.Q + b.Q ~ 0
    end
end
# ... Valve (diodo suave), Ventricle (E(t)) e montagem com connect(...)
# sys = mtkcompile(System(eqs, t; systems = [...]))
```

> A API do MTK mudou entre versões (`ODESystem` → `System`, `structural_simplify` → `mtkcompile`).
> Confira a versão no `Manifest.toml` antes de copiar código de tutoriais. A vantagem do MTK é a
> modularidade (adicionar um Windkessel de 3 elementos é trocar um componente). A da função
> "à mão" é a transparência e a facilidade de depuração. **Sugestão:** desenvolva primeiro à mão
> (§13.2), valide com os testes da §9 e só então migre para MTK, usando o modelo à mão como gabarito
> (resultados idênticos).

---

## 14. Figuras esperadas (entregáveis A.2 e A.3)

1. Windkessel: P_a(t) e o ajuste exponencial na diástole (τ).
2. Malha fechada basal: P_VE, P_SA e P_PA × t; Q_ao × t; **PV loops VE e VD**.
3. Curvas de Guyton (VR e CO × RAP) com pontos de funcionamento: basal × hemorragia × terapias.
4. Hemorragia classes I–IV: PAM, CO, HR × t (barorreflexo ON/OFF).
5. Terapias: PAM, CO, Hct, μ, DO₂ e V_is × t para RL × sangue total × vasopressor.
6. Mapas de sensibilidade (heatmap S_ij) e ajuste PAM_sim × PAM_med (MIMIC) com resíduos.

---

## 15. Limitações a declarar

- 0D: sem propagação/reflexão de onda e sem distribuição regional de fluxo (órgãos).
- Coração linear (E_max/E_min) e sem átrios ativos na primeira versão; sem acoplamento ventricular.
- Fluido intersticial e viscosidade com relações empíricas (dependentes da tabela da Frente B).
- Barorreflexo simplificado (sem quimiorreflexo, sem hormônios, sem regulação renal de longo
  prazo). Para horizontes de horas a dias, ver o escopo de modelos integrativos como o
  Guyton/HumMod, citados por Hester.
- Calibração com dados esparsos do MIMIC (PAM intermitente/horária, horários de infusão
  imprecisos): identificabilidade limitada e incerteza nos parâmetros.
- O COU é educacional/exploratório. Um "gêmeo digital" clínico exigiria o V&V 40 completo com
  quantificação de incerteza.

---

## 16. Referências do material

- Aletti F. *Módulo Cardiovascular* (slides MSB 2026: Casos #1–#8) e Aulas 1, 4, 5.
- Pré-prova P1 (15/09/2026): tópicos que o modelo deve conseguir explicar.
- Pruett WA, Clemmer JS, Hester RL. *Physiological Modeling and Simulation — Validation,
  Credibility, and Application.* Annu Rev Biomed Eng 2020;22:185–206.
- Nucci G, Cobelli C. *Mathematical Models of Respiratory Mechanics* (modelos de 1ª/2ª ordem,
  identificação por mínimos quadrados, identificabilidade).
- Cobelli C, Carson E. *Introduction to Modeling in Physiology and Medicine.*
- Programa de Simulação Realística (CSR, 6º sem.): Cenário A (choque séptico/hipovolemia
  relativa peri-intubação, PEEP) e Cenário C (noradrenalina + Ringer Lactato) como casos de teste
  qualitativos.
- Abstract PBL 2026 – G3 (escopo, ferramentas e cronograma).
