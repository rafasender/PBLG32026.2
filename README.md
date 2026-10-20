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
├── .gitignore               # Ignora temporários, Manifest e dados brutos pesados
├── Project.toml             # Manifesto de dependências do ambiente Julia
├── Manifest.toml            # Árvore exata de versões resolvidas (gerado localmente)
├── setup_environment.jl     # Script automatizado de instalação das dependências
├── README.md                # Este documento
├── data/
│   ├── raw/                 # Extrações brutas do MIMIC-IV (chartevents, inputevents, etc.)
│   └── processed/           # Coortes filtradas e séries temporais alinhadas
├── docs/                    # Documentação técnica, artigos de referência e notas clínicas
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

## 🚀 Instalação e Configuração do Ambiente

### 1. Pré-requisitos
- **Julia** (versão 1.10 ou superior):
  - No Windows, recomenda-se instalar via terminal utilizando o **Juliaup**:
    ```powershell
    winget install Julialang.Juliaup
    juliaup add release
    ```
  - No Linux/macOS:
    ```bash
    curl -fsSL https://install.julialang.org | sh
    ```
- **Git** instalado na máquina.

### 2. Clonar o Repositório
```bash
git clone https://github.com/rafasender/PBLG32026.2.git
cd PBLG32026.2
```

### 3. Instanciar o Ambiente de Dependências
Para garantir que todos os membros do grupo utilizem exatamente as mesmas versões de bibliotecas, execute no terminal:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
```

Ou execute o script interativo de configuração:
```bash
julia setup_environment.jl
```

---

## 🧪 Verificação e Demonstração

### Teste de Sanidade do Ambiente
Para verificar se todos os pacotes do ecossistema SciML, tratamento de dados e gráficos estão operacionais, execute:
```bash
julia scripts/check_env.jl
```

### Executar a Suíte de Testes Unitários
```bash
julia --project=. test/runtests.jl
```
*(ou via REPL do Julia: digite `]` para entrar no modo Pkg e execute `test`)*

### Executar a Simulação Demonstrativa (Windkessel 2 Elementos)
Simula a curva de pressão aórtica em normovolemia vs. choque hipovolêmico com taquicardia e vasoconstrição reflexa, gerando a figura em `output/demo_windkessel2.png`:
```bash
julia scripts/demo_windkessel2.jl
```

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
