# scripts/check_env.jl
# Script de verificação de sanidade (Smoke Test) do ambiente Julia

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

println("="^60)
println(" TESTE DE SANIDADE DO AMBIENTE JULIA - PBL 2026.2 (G3)")
println("="^60)

println("\n1. Testando importação dos pacotes científicos e SciML...")
try
    @time using DifferentialEquations
    println("   [OK] DifferentialEquations carregado com sucesso.")
catch e
    println("   [FALHA] DifferentialEquations: $e")
end

try
    @time using ModelingToolkit
    println("   [OK] ModelingToolkit carregado com sucesso.")
catch e
    println("   [FALHA] ModelingToolkit: $e")
end

try
    @time using Optimization
    @time using OptimizationOptimJL
    println("   [OK] Optimization + OptimizationOptimJL carregados com sucesso.")
catch e
    println("   [FALHA] Optimization: $e")
end

try
    @time using SciMLSensitivity
    println("   [OK] SciMLSensitivity carregado com sucesso.")
catch e
    println("   [FALHA] SciMLSensitivity: $e")
end

println("\n2. Testando importação dos pacotes de análise de dados (MIMIC-IV)...")
try
    @time using DataFrames
    @time using CSV
    @time using Chain
    @time using CategoricalArrays
    println("   [OK] DataFrames, CSV, Chain e CategoricalArrays carregados com sucesso.")
catch e
    println("   [FALHA] Pacotes de dados: $e")
end

println("\n3. Testando biblioteca de gráficos...")
try
    @time using Plots
    println("   [OK] Plots (GR backend) carregado com sucesso.")
catch e
    println("   [FALHA] Plots: $e")
end

println("\n" * "="^60)
println(" [CONCLUÍDO] Verificação do ambiente finalizada!")
println("="^60)
