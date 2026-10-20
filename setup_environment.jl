# setup_environment.jl
# Script automatizado para configurar e instanciar o ambiente Julia do PBL 2026.2

using Pkg

println("="^70)
println(" PBL 2026.2 - Grupo 3: Configuração do Ambiente Julia (Entregável A.1)")
println("="^70)

# Ativa o diretório atual do projeto
Pkg.activate(@__DIR__)

println("\n--> Verificando e instalando pacotes fundamentais do projeto...")

packages = [
    "DifferentialEquations",
    "ModelingToolkit",
    "Optimization",
    "OptimizationOptimJL",
    "CSV",
    "DataFrames",
    "Chain",
    "CategoricalArrays",
    "Plots"
]

for pkg in packages
    println("--> Garantindo dependência: $pkg")
    try
        Pkg.add(pkg)
    catch e
        @warn "Aviso ao adicionar $pkg: $e"
    end
end

println("\n--> Instanciando e pré-compilando o ambiente...")
Pkg.instantiate()
Pkg.precompile()

println("\n" * "="^70)
println(" Status do Ambiente (Project.toml):")
println("="^70)
Pkg.status()

println("\n[SUCESSO] Ambiente configurado e pronto para a Frente A, B e C!")
