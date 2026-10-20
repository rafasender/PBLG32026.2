# test/runtests.jl
using Test

# Carregar módulo local
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))
using Cardio0D

@testset "Cardio0D.jl - Testes do Módulo de Hemodinâmica 0D" begin

    @testset "Reologia e Viscosidade Sanguínea μ(Hct)" begin
        # Teste com hematócrito fisiológico basal (~40%)
        mu_norm = blood_viscosity(0.40; mu_plasma=1.20)
        @test mu_norm > 1.20
        @test isapprox(mu_norm, 1.20 * exp(2.5 * 0.40), atol=1e-4)

        # Teste com hemodiluição por cristaloide (queda no Hct para 25%)
        mu_diluted = blood_viscosity(0.25; mu_plasma=1.20)
        @test mu_diluted < mu_norm

        # Teste com struct BloodProperties
        bp = BloodProperties(hct=0.45, mu_plasma=1.25)
        @test blood_viscosity(bp) ≈ 1.25 * exp(2.5 * 0.45)
    end

    @testset "Parâmetros do Windkessel de 2 Elementos" begin
        params = Windkessel2Params(R=1.0, C=1.2, Pv=4.0)
        @test params.R == 1.0
        @test params.C == 1.2
        @test params.Pv == 4.0

        # Teste da derivada da EDO no repouso (sem fluxo Qin = 0)
        dP = [0.0]
        P = [104.0] # Pa = 104 mmHg
        Qin_zero = (t) -> 0.0
        windkessel2_ode!(dP, P, (params, Qin_zero), 0.0)
        
        # dP/dt deve ser (0 - (104 - 4)/1.0) / 1.2 = -100 / 1.2 ≈ -83.33 mmHg/s
        @test isapprox(dP[1], -100.0 / 1.2, atol=1e-4)
    end

    @testset "Cálculo de PAM (Pressão Arterial Média)" begin
        times = [0.0, 0.5, 1.0]
        pressures = [80.0, 100.0, 80.0]
        # Pela regra trapezoidal: 0.5*(80+100)*0.5 + 0.5*(100+80)*0.5 = 45 + 45 = 90
        pam = calculate_map(times, pressures)
        @test isapprox(pam, 90.0, atol=1e-4)
    end
end
