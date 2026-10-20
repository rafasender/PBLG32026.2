# setup_environment.jl
# Script automatizado para instanciar e pré-compilar o ambiente Julia do PBL 2026.2

using Pkg

println("="^70)
println(" PBL 2026.2 - Grupo 3: Configuração do Ambiente Julia (Entregável A.1)")
println("="^70)

# Ativa o diretório atual do projeto
Pkg.activate(@__DIR__)

println("\n--> Instanciando o ambiente a partir de Project.toml e Manifest.toml...")
Pkg.instantiate()

println("\n--> Pré-compilando pacotes instalados...")
Pkg.precompile()

println("\n" * "="^70)
println(" Status do Ambiente (Project.toml):")
println("="^70)
Pkg.status()

println("\n[SUCESSO] Ambiente configurado e sincronizado com sucesso para as Frentes A, B e C!")
