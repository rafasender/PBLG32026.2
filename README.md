# Modelagem 0D do Sistema Cardiovascular no Choque Hipovolêmico
### Comparação de Terapias de Reposição Volêmica em UTI com Calibração no MIMIC-IV

[![Julia Version](https://img.shields.io/badge/julia-v1.10%2B-blue.svg)](https://julialang.org)
[![SciML](https://img.shields.io/badge/Powered%20by-SciML-9558B2.svg)](https://sciml.ai)
[![Project](https://img.shields.io/badge/PBL-2026.2--G3-green.svg)](#)

---

## 📌 Contexto Clínico & Objetivo

O **choque hipovolêmico** é uma das causas mais prevalentes de colapso hemodinâmico e mortalidade em Unidades de Terapia Intensiva (UTI). A perda abrupta de volume intravascular reduz o retorno venoso e a pré-carga, levando à queda do volume sistólico e do débito cardíaco ($DC = FC \times VS$), comprometendo a oxigenação tecidual ($DO_2$).

A conduta terapêutica envolve opções com impactos hemodinâmicos e reológicos distintos:
- **Cristaloides (ex: Ringer Lactato)**: Expansão volêmica rápida e transitória, porém com migração substancial para o interstício e indução de hemodiluição acentuada, diminuindo o hematócrito ($\text{Hct}$) e a viscosidade sanguínea ($\mu$), além do risco de coagulopatia dilucional e acidose.
- **Sangue Total / Hemocomponentes**: Reposição conjunta de volume, hemácias, plasma e plaquetas, preservando o hematócrito, a viscosidade sanguínea e a capacidade de transporte de $O_2$.
- **Vasopressores (ex: Noradrenalina)**: Elevam a resistência vascular sistêmica ($R$) e a pressão arterial, porém sem restabelecer a volemia efetiva.

### 🎯 Objetivo do Projeto
Desenvolver um **modelo matemático-computacional 0D (parâmetros concentrados tipo Windkessel)** do sistema cardiovascular acoplado em malha fechada, calibrado com dados clínicos reais da base **MIMIC-IV**, capaz de quantificar o impacto fisiológico e comparar a resposta hemodinâmica às diferentes terapias de expansão volêmica.

---

## 👥 Divisão das Frentes de Trabalho (PBL 2026.2 - Grupo 3)

| Frente | Descrição | Principais Entregáveis | Integrantes Responsáveis |
| :--- | :--- | :--- | :--- |
| **Frente A – Modelo** | Modelagem cardiovascular em Julia com `ModelingToolkit.jl` e `DifferentialEquations.jl` | • Ambiente Julia configurado<br>• Windkessel arterial, venoso, pulmonar e cardíaco em malha fechada | André Miyake |
| **Frente B – Fisiologia e Cenários** | Parâmetros da literatura, reologia $\mu(\text{Hct})$, choque e protocolos de reposição | • Resumo teórico e tabela de parâmetros<br>• Equações de viscosidade $\mu(\text{Hct})$<br>• Critérios clínicos e discussão médica | Eduarda Klapper, Luana Metz |
| **Frente C – Dados e Calibração** | Credenciamento PhysioNet, extração MIMIC-IV, otimização e sensibilidade | • Pipeline de dados (ICU / Hosp)<br>• Ajuste de $R$ e $C$ via `Optimization.jl`<br>• Análise de sensibilidade via `SciMLSensitivity.jl` | Maria Eduarda Costa, Henrique Martins, Rafael Sender |

---

## 📂 Estrutura do Repositório

```text
├── .gitignore               # Ignora temporários, dados tabulares e arquivos sensíveis
├── Project.toml             # Manifesto de dependências do ambiente Julia
├── Manifest.toml            # Árvore exata de versões resolvidas (versionado para reprodutibilidade)
├── setup_environment.jl     # Script automatizado de instanciação do ambiente
├── README.md                # Este documento
├── data/
│   ├── raw/                 # Extrações brutas do MIMIC-IV (chartevents, inputevents, etc.)
│   └── processed/           # Coortes filtradas e séries temporais alinhadas
├── output/                  # Gráficos gerados, curvas de pressão e relatórios
├── scripts/
│   ├── check_env.jl         # Smoke test: verifica carregamento de todos os pacotes
│   └── demo_windkessel2.jl  # Simulação demonstrativa do Windkessel de 2 elementos
├── src/
│   └── Cardio0D.jl          # Módulo principal (EDOs, reologia μ(Hct), Windkessel)
└── test/
    └── runtests.jl          # Suíte de testes unitários automatizados
```

---

## 📖 Tutorial Passo a Passo: Como Rodar o Código

Este tutorial foi feito para que qualquer pessoa (mesmo sem experiência prévia com Julia) consiga clonar, configurar e rodar o projeto do zero sem erros.

---

### ⚡ Modo Rápido (TL;DR)
Se você já tem o Julia instalado e só quer rodar tudo:

```bash
# 1. Clonar e entrar na pasta
git clone https://github.com/rafasender/PBLG32026.2.git
cd PBLG32026.2

# 2. Configurar o ambiente com as versões exatas
julia setup_environment.jl

# 3. Rodar os testes unitários
julia --project=. test/runtests.jl

# 4. Rodar a simulação e gerar o gráfico
julia scripts/demo_windkessel2.jl
```

---

### 🛠️ Guia Detalhado e Explicado

#### Passo 1: Instalação do Julia
Recomenda-se a versão **Julia v1.10 LTS** ou superior (o projeto é 100% compatível com v1.10, v1.11 e v1.12+).

* **No Windows (recomendado via terminal PowerShell):**
  ```powershell
  winget install Julialang.Juliaup
  juliaup add release
  ```
  *(Ou baixe o instalador oficial diretamente em [julialang.org/downloads](https://julialang.org/downloads/))*.
* **No Linux / macOS:**
  ```bash
  curl -fsSL https://install.julialang.org | sh
  ```
* **Verifique a instalação abrindo um terminal e digitando:**
  ```bash
  julia --version
  # Deve retornar algo como: julia version 1.10.x ou superior
  ```

---

#### Passo 2: Clonar o Repositório e Navegar até a Pasta
Abra o seu terminal (PowerShell, Git Bash ou terminal do VS Code) e execute:

```bash
git clone https://github.com/rafasender/PBLG32026.2.git
cd PBLG32026.2
```

> **Atenção:** Certifique-se de que o terminal está apontado para a raiz do repositório (`PBLG32026.2`).

---

#### Passo 3: Configurar o Ambiente de Pacotes (`setup_environment.jl`)
Para que o código funcione de forma idêntica na máquina de todos os membros do grupo, as dependências são travadas pelo arquivo `Manifest.toml`.

Execute o script de configuração:
```bash
julia setup_environment.jl
```

**O que esse comando faz nos bastidores?**
1. Ativa o ambiente do projeto atual.
2. Executa `Pkg.instantiate()`, que lê o arquivo `Manifest.toml` e baixa as versões exatas das bibliotecas (DifferentialEquations, ModelingToolkit, SciMLSensitivity, Plots, etc.).
3. Executa `Pkg.precompile()`, compilando antecipadamente as bibliotecas.

> ⏳ **Nota sobre a primeira execução:** A primeira vez que você rodar esse script pode demorar entre 2 a 5 minutos, pois pacotes de diferenciação automática e compilação do SciML (como Enzyme e LLVM) estão sendo preparados. Nas próximas execuções, o carregamento será quase instantâneo.

---

#### Passo 4: Teste de Sanidade do Ambiente (`check_env.jl`)
Para garantir que todos os pacotes científicos, motores de dados e biblioteca de gráficos foram instalados sem falhas de compilação, execute o *smoke test*:

```bash
julia scripts/check_env.jl
```

**Resultado esperado no terminal:**
```text
============================================================
 TESTE DE SANIDADE DO AMBIENTE JULIA - PBL 2026.2 (G3)
============================================================

1. Testando importação dos pacotes científicos e SciML...
   [OK] DifferentialEquations carregado com sucesso.
   [OK] ModelingToolkit carregado com sucesso.
   [OK] Optimization + OptimizationOptimJL carregados com sucesso.
   [OK] SciMLSensitivity carregado com sucesso.

2. Testando importação dos pacotes de análise de dados (MIMIC-IV)...
   [OK] DataFrames, CSV, Chain e CategoricalArrays carregados com sucesso.

3. Testando biblioteca de gráficos...
   [OK] Plots (GR backend) carregado com sucesso.

============================================================
 [CONCLUÍDO] Verificação do ambiente finalizada!
============================================================
```

---

#### Passo 5: Executar os Testes Unitários (`runtests.jl`)
Para checar a integridade matemática do código (reologia de viscosidade, derivada da EDO do Windkessel e cálculo numérico da Pressão Arterial Média):

```bash
julia --project=. test/runtests.jl
```

> 💡 **Por que usar `--project=.`?**
> A flag `--project=.` avisa ao Julia para usar as dependências desta pasta. Como este repositório é uma aplicação clínica/pesquisa (e não uma biblioteca registrada no registro oficial do Julia), os testes são executados diretamente através do arquivo `test/runtests.jl` (não utilize `] test` no REPL).

**Resultado esperado:**
```text
Test Summary:                                     | Pass  Total  Time
Cardio0D.jl - Testes do Módulo de Hemodinâmica 0D |    9      9  0.2s
```

---

#### Passo 6: Executar a Simulação Demonstrativa (`demo_windkessel2.jl`)
Este script simula a hemodinâmica arterial do sistema cardiovascular em regime estacionário (30 ciclos cardíacos) comparando dois estados clínicos:
1. **Linha de Base (Normovolemia)**: $FC = 75\text{ bpm}$, $VS = 70\text{ mL}$, $DC = 5{,}25\text{ L/min}$.
2. **Choque Hipovolêmico Compensado**: $FC = 110\text{ bpm}$, $VS = 35\text{ mL}$, $DC = 3{,}85\text{ L/min}$ (queda de $26{,}7\%$ no débito cardíaco, porém com vasoconstrição reflexa elevando $R$ de $1{,}0$ para $1{,}5\text{ mmHg}\cdot\text{s/mL}$).

Execute:
```bash
julia scripts/demo_windkessel2.jl
```

**O que acontece ao rodar?**
1. O terminal exibirá as métricas hemodinâmicas calculadas no estado estacionário:
   ```text
   --- Cenário 1: Linha de Base (Normovolemia) ---
     Frequência Cardíaca (FC): 75.0 bpm
     Volume Sistólico (VS)   : 70.0 mL
     Débito Cardíaco (DC)    : 5.25 L/min
     Pressão Sistólica (PAS) : 111.2 mmHg
     Pressão Diastólica (PAD): 73.3 mmHg
     PAM (via integral)      : 91.5 mmHg
     PAM (fórmula clínica)   : 86.0 mmHg

   --- Cenário 2: Choque Hipovolêmico Compensado (Resposta Simpática) ---
     Frequência Cardíaca (FC): 110.0 bpm
     Volume Sistólico (VS)   : 35.0 mL
     Débito Cardíaco (DC)    : 3.85 L/min (redução de 26.7%)
     Pressão Sistólica (PAS) : 108.6 mmHg
     Pressão Diastólica (PAD): 88.2 mmHg
     PAM calculada           : 98.3 mmHg
     [Nota Fisiológica]: A PAM é mantida acima do basal (~98.3 vs ~91.5 mmHg)
                         devido ao aumento reflexo da RVS (1.0 -> 1.5 mmHg·s/mL), caracterizando o choque
                         compensado. A pressão aparentemente normal/elevada mascara o colapso do DC.
   ```
2. Um gráfico de alta resolução será salvo automaticamente em:
   📁 **`output/demo_windkessel2.png`**

---

### 💻 Usando pelo Visual Studio Code (VS Code)

Se preferir usar o editor visualmente:
1. Abra o VS Code e instale a extensão oficial **"Julia"** (desenvolvida pela equipe do Julia).
2. Abra a pasta do projeto no VS Code: `File -> Open Folder -> PBLG32026.2`.
3. Na barra inferior azul do VS Code, clique no ambiente do Julia e confirme que ele está selecionado como `(PBL 2026.2)` ou `./`.
4. Abra qualquer arquivo (ex: `scripts/demo_windkessel2.jl`) e pressione **Shift + Enter** para rodar linhas interativamente ou **Alt + J Alt + O** para abrir o REPL integrado.

---

### ❓ Perguntas Frequentes & Resolução de Problemas (Troubleshooting)

* **Erro: `julia: command not found` ou `julia não é reconhecido como um comando`**
  * *Causa:* O Julia não foi adicionado ao `PATH` do sistema.
  * *Solução:* Feche e reabra o terminal após a instalação. No Windows, use `winget install Julialang.Juliaup` que configura o `PATH` automaticamente.
* **Por que não rodar `] test` no REPL do Julia?**
  * *Causa:* O comando `] test` no REPL exige que o projeto tenha um `name` e `uuid` registrados no `Project.toml` como um pacote público/distribuível.
  * *Solução:* Como o repositório é um projeto de aplicação científica com módulo local, o comando correto é `julia --project=. test/runtests.jl`.
* **Como limpar saídas ou resetar o ambiente?**
  * Os gráficos gerados ficam em `output/`. Você pode excluí-los com segurança a qualquer momento; eles serão recriados ao executar o script.
  * Para forçar uma reinstalação limpa dos pacotes: `julia --project=. -e "using Pkg; Pkg.instantiate()"`

---

## 📚 Bibliotecas Principais Utilizadas

- **Modelagem e Resolução Numérica**:
  - [`DifferentialEquations.jl`](https://docs.sciml.ai/DiffEqDocs/stable/): Solvers de alta performance para equações diferenciais rígidas e não-rígidas.
  - [`ModelingToolkit.jl`](https://docs.sciml.ai/ModelingToolkit/stable/): Modelagem baseada em componentes acausais (R, C, bombas hidráulicas).
- **Otimização e Sensibilidade**:
  - [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/) & [`OptimizationOptimJL.jl`](https://github.com/SciML/Optimization.jl): Ajuste paramétrico de $R$ e $C$ por minimização de erro quadrático em relação à PAM do MIMIC-IV.
  - [`SciMLSensitivity.jl`](https://docs.sciml.ai/SciMLSensitivity/stable/): Análise de sensibilidade via gradientes adjuntos.
- **Tratamento e Engenharia de Dados**:
  - [`DataFrames.jl`](https://dataframes.juliadata.org/stable/), [`CSV.jl`](https://csv.juliadata.org/stable/), [`Chain.jl`](https://github.com/jkrumbiegel/Chain.jl).
- **Visualização**:
  - [`Plots.jl`](https://docs.juliaplots.org/latest/) com backend GR.

---

## 📖 Referências Bibliográficas

1. **Cannon, J. W.** Hemorrhagic shock. *New England Journal of Medicine*, v. 378, n. 4, p. 370–379, 2018.
2. **Hall, J. E.; Hall, M. E.** *Guyton & Hall – Tratado de Fisiologia Médica*. 14. ed. Elsevier, 2021.
3. **Pries, A. R. et al.** Blood viscosity in tube flow: dependence on diameter and hematocrit. *Am J Physiol*, 1992.
4. **Westerhof, N. et al.** The arterial Windkessel. *Med Biol Eng Comput*, v. 47, n. 2, p. 131–141, 2009.
5. **Johnson, A. E. W. et al.** MIMIC-IV, a freely accessible electronic health record dataset. *Sci Data*, 2023.


---

## Modelo implementado (resumo)

O módulo `Cardio0D` (`src/`) contém o modelo 0D em malha fechada (pulsátil e de valores médios),
hemorragia, terapias (Ringer Lactato, sangue total, vasopressor), barorreflexo, leitura do
esquema MIMIC-IV e calibração. Veja o [relatório](docs/RELATORIO.md) e o guia matemático
[CLAUDE.md](CLAUDE.md).

```bash
julia --project=. test/runtests.jl                 # 85 testes (~25 s)
julia --project=. scripts/run_model.jl             # figuras e tabelas em docs/
julia --project=. scripts/calibrate_synthetic.jl   # calibração com paciente sintético (~6 min)
```
