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

        # Pries et al. (1992): Hct 0,45 em vaso grande dá μ_rel = μ_0.45 ≈ 3,19 (→ 3,2 quando D → ∞)
        @test isapprox(blood_viscosity_pries(0.45; D=1000.0, invivo=false) / 1.20, 3.19, atol=0.01)
        @test blood_viscosity_pries(0.0; D=100.0, invivo=false) ≈ 1.20   # sem hemácias: plasma
        # valores da tabela B.2 (Frente B), μ_rel com Hct = 0,40
        @test isapprox(blood_viscosity_pries(0.40; D=1000.0, invivo=false) / 1.2, 2.80; atol=0.01)
        @test isapprox(blood_viscosity_pries(0.40; D=100.0) / 1.2, 2.26; atol=0.01)
        @test isapprox(blood_viscosity_pries(0.40; D=50.0) / 1.2, 2.06; atol=0.01)
        # hemodiluição 40 % -> 25 %: razão 0,687 (vaso grande) e 0,733 (D = 100 µm in vivo)
        @test isapprox(blood_viscosity(0.25) / blood_viscosity(0.40), 0.687; atol=0.001)
        @test isapprox(blood_viscosity_pries(0.25; D=100.0) / blood_viscosity_pries(0.40; D=100.0), 0.733; atol=0.001)
        # a razão μ(Hct)/μ(0,42) do modelo ≈ Pries em vaso grande (3%) para Hct 0,15–0,55
        for h in 0.15:0.05:0.55
            r_exp = blood_viscosity(h) / blood_viscosity(0.42)
            r_pr = blood_viscosity_pries(h; D=1000.0, invivo=false) / blood_viscosity_pries(0.42; D=1000.0, invivo=false)
            @test isapprox(r_exp, r_pr; rtol=0.03)
        end
        # Fåhræus–Lindqvist: em microvasos o Hct pesa menos
        @test blood_viscosity_pries(0.42; D=30.0) < blood_viscosity_pries(0.42; D=1000.0, invivo=false)
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

# =====================================================================================
# Verificação do modelo em malha fechada (CLAUDE.md §9: V1–V8) e da camada de dados
# =====================================================================================
using DifferentialEquations, DataFrames, CSV
using Dates: DateTime, Second

const P0 = baseline_params()

@testset "V1 - Windkessel: τ ajustado na diástole = R·C" begin
    prm = Windkessel2Params(R=1.1, C=1.4, Pv=0.0)
    sol = solve(ODEProblem(windkessel2_ode!, [100.0], (0.0, 5.0), (prm, t -> 0.0)), Tsit5();
                reltol=1e-10, abstol=1e-10)
    ts = range(0.2, 4.0; length=50)
    y = [log(sol(t)[1]) for t in ts]
    slope = (sum((ts .- sum(ts) / 50) .* (y .- sum(y) / 50))) / sum((ts .- sum(ts) / 50) .^ 2)
    @test isapprox(-1 / slope, prm.R * prm.C; rtol=0.01)
end

@testset "Peças elementares" begin
    # ativação: contínua, em [0,1], zero na diástole e 1 no pico
    T = 0.8; Ts1 = 0.3 * sqrt(T)
    @test activation(0.0, T) == 0.0
    @test isapprox(activation(Ts1, T), 1.0; atol=1e-12)
    @test activation(0.7, T) == 0.0
    @test all(0 .<= activation.(range(0, T; length=200), T) .<= 1)
    # válvula: fecha em ΔP<0, é a lei de Ohm em ΔP grande
    @test abs(valve_flow(-5.0, 0.005, 10.0)) < 1e-12
    @test isapprox(valve_flow(5.0, 0.005, 10.0), 5.0 / 0.005; rtol=1e-6)
    @test isfinite(valve_flow(-1e4, 0.005, 10.0))
    # pressão venosa e seu inverso
    for P in (0.5, 2.0, 5.0, 12.0)
        s = Cardio0D.stressed_volume_from_pressure(P, 100.0, 0.25)
        @test isapprox(Cardio0D.venous_pressure(s, 100.0, 0.25), P; atol=1e-9)
    end
    @test Cardio0D.venous_pressure(-500.0, 100.0, 0.25) < 1e-6   # veia colapsada
    # with_params troca só o que foi pedido
    q = with_params(P0; Rsa=2.0)
    @test q.Rsa == 2.0 && q.Csa == P0.Csa
    # viscosidade: hemodiluição reduz μ
    @test blood_viscosity(0.28) < blood_viscosity(0.42)
    @test blood_viscosity_krieger(0.28) < blood_viscosity_krieger(0.42)
end

@testset "V8 - Estado basal reproduz os alvos (modelo lento)" begin
    op = operating_point(P0)
    @test isapprox(op.Psa, 93.0; atol=1e-6)
    @test isapprox(op.CO, 5.25; atol=1e-6)
    @test isapprox(op.Psv, 5.0; atol=1e-6)
    @test isapprox(op.RAP, 2.5; atol=1e-6)
    @test isapprox(op.Ppa, 15.0; atol=1e-6)
    @test isapprox(op.Ppv, 8.0; atol=1e-6)
    @test 6.0 < op.Pms < 8.0                    # Pms ≈ 7 mmHg na aula
    # o basal é um equilíbrio da EDO
    u0 = initial_state_mean(P0); du = zeros(length(u0))
    mean_rhs!(du, u0, (P0, Scenario()), 0.0)
    @test maximum(abs, du[1:4]) < 1e-6
end

@testset "V3 - Coração parado: pressões sistêmicas → Pms (E8)" begin
    sc = Scenario(hr_forced = t -> 0.0, baro = false)
    u0 = initial_state_mean(P0)
    Pms = ((u0[1] + u0[2]) - (P0.V0sa + P0.V0sv)) / (P0.Csa + P0.Csv)
    sol = simulate_mean(P0, sc, (0.0, 600.0); u0)
    h = mean_hemodynamics(P0, sol.u[end], sc, 600.0)
    @test isapprox(h.Psa, Pms; atol=0.05)
    @test isapprox(h.Psv, Pms; atol=0.05)
    @test isapprox(h.Ppa, h.Ppv; atol=0.05)
    @test abs(h.QL) < 1e-6
end

@testset "V2 - Conservação de volume" begin
    # mean: Σ V_vasos + Vis − (Σ V_vasos(0)) = ∫(Qinf − Qhem) dt, com troca capilar ativa
    sc = Scenario(Qhem = window_rate(900.0, 100.0, 300.0), Qinf = window_rate(500.0, 600.0, 200.0),
                  hct_inf = 0.0, tstops = [100.0, 400.0, 600.0, 800.0])
    u0 = initial_state_mean(P0)
    sol = simulate_mean(P0, sc, (0.0, 2000.0); u0)
    net = 500.0 - 900.0
    @test isapprox(sum(sol.u[end][1:4]) + sol.u[end][6] - sum(u0[1:4]), net; atol=1e-3)
    # pulsátil: mesma conservação incluindo os ventrículos
    scp = Scenario(Qhem = window_rate(300.0, 2.0, 5.0), baro = false, tstops = [2.0, 7.0])
    up = initial_state_pulsatile(P0)
    solp = simulate_pulsatile(P0, scp, (0.0, 20.0); u0=up)
    @test isapprox(sum(solp.u[end][1:6]) + solp.u[end][9] - sum(up[1:6]), -300.0; atol=0.5)
end

@testset "V4 - Ventrículo isolado: SV = C_D·P_f − C_S·P_a (E13)" begin
    Pf, Pa, Vd, HR = 8.0, 93.0, 10.0, 75.0
    CD, CS = 1 / P0.EminL, 1 / P0.EmaxL
    T = 60 / HR
    f!(du, u, p, t) = (P = elastance(P0.EminL, P0.EmaxL, activation(t - T * floor(t / T), T)) * (u[1] - Vd);
                       du[1] = valve_flow(Pf - P, P0.Rvalv, P0.kvalv) - valve_flow(P - Pa, P0.Rvalv, P0.kvalv))
    sol = solve(ODEProblem(f!, [Vd + CD * Pf], (0.0, 25 * T)), Rodas5P(); reltol=1e-8, abstol=1e-8)
    ts = range(24 * T, 25 * T; length=800)
    V = [sol(t)[1] for t in ts]
    SV = maximum(V) - minimum(V)
    @test isapprox(SV, CD * Pf - CS * Pa; rtol=0.05)
end

@testset "V5 - Ponto de funcionamento: algébrico × Guyton × EDO" begin
    op0 = operating_point(P0)
    g0 = guyton_intersection(P0, op0)
    @test isapprox(g0.RAP, op0.RAP; atol=0.02)
    @test isapprox(g0.CO, op0.CO; rtol=0.005)

    # hemorragia de 700 mL sem barorreflexo nem interstício: a EDO deve ir ao equilíbrio algébrico
    pk = with_params(P0; kf=0.0)
    sc = Scenario(Qhem = window_rate(700.0, 60.0, 300.0), baro = false, tstops = [60.0, 360.0])
    sol = simulate_mean(pk, sc, (0.0, 4000.0))
    h = mean_hemodynamics(pk, sol.u[end], sc, 4000.0)
    op1 = operating_point(pk; Vb = pk.Vb0 - 700.0)
    @test isapprox(h.Psa, op1.Psa; atol=0.1)
    @test isapprox(mls_to_lpm(h.QL), op1.CO; rtol=0.01)
    @test isapprox(h.RAP, op1.RAP; atol=0.1)
    # Guyton (E10) vale no regime linear: hemorragia de 400 mL, antes do colapso venoso
    op4 = operating_point(pk; Vb = pk.Vb0 - 400.0)
    g4 = guyton_intersection(pk, op4)
    @test isapprox(g4.RAP, op4.RAP; atol=0.05)
    @test isapprox(g4.CO, op4.CO; rtol=0.01)
    # a hemorragia desloca a curva de VR (↓Pms) e baixa RAP e CO
    @test op1.Pms < op0.Pms && op1.RAP < op0.RAP && op1.CO < op0.CO
    # resistência alta gira a curva (↓VR com o mesmo Pms)
    @test venous_return_curve(P0, op0, 2.5; Rsa = 2 * P0.Rsa) < venous_return_curve(P0, op0, 2.5)
end

@testset "V6/V8 - Modelo pulsátil: regime periódico e coerência com o modelo lento" begin
    sc = Scenario(baro = false)
    sol = simulate_pulsatile(P0, sc, (0.0, 40.0))
    @test sol.retcode == SciMLBase.ReturnCode.Success
    m = cycle_metrics(pulsatile_observables(sol, P0, sc, last_beats(sol, P0; n=3)))
    # periodicidade: duas janelas consecutivas dão o mesmo resultado
    sc_ = Scenario(baro = false)
    ts2 = range(sol.t[end] - 6 * 0.8, sol.t[end] - 3 * 0.8; step=0.005)
    m2 = cycle_metrics(pulsatile_observables(sol, P0, sc_, ts2))
    @test isapprox(m.PAM, m2.PAM; atol=0.2)
    # faixas fisiológicas (V8)
    @test 80 < m.PAM < 100
    # Csa = 2,5 mL/mmHg (tabela B.2) dá pressão de pulso ~20 mmHg (< ~40 clínico): limitação do
    # 0D de 2 elementos, o alvo do projeto é a PAM (CLAUDE.md §1). Por isso PAS > 90, não 100.
    @test 90 < m.PAS < 135 && 60 < m.PAD < 85
    @test 4.3 < m.CO < 6.0
    @test 55 < m.SV < 85
    @test 0.5 < m.EF < 0.75
    @test 2 < m.CVP < 8
    # V6: coerência com o modelo lento
    op = operating_point(P0)
    @test abs(m.PAM - op.Psa) < 6.0
    @test isapprox(m.CO, op.CO; rtol=0.10)
    # PV loop fechado: o último batimento volta ao ponto de partida
    V, P = pv_loop(pulsatile_observables(sol, P0, sc, range(sol.t[end] - 1.6, sol.t[end]; step=0.005)), P0)
    @test abs(V[1] - V[end]) < 3.0
end

@testset "V7 - Convergência numérica" begin
    sc = Scenario(Qhem = window_rate(1000.0, 60.0, 300.0), tstops = [60.0, 360.0])
    a = simulate_mean(P0, sc, (0.0, 1800.0); reltol=1e-6, abstol=1e-6)
    b = simulate_mean(P0, sc, (0.0, 1800.0); reltol=1e-9, abstol=1e-9)
    ha = mean_hemodynamics(P0, a.u[end], sc, 1800.0); hb = mean_hemodynamics(P0, b.u[end], sc, 1800.0)
    @test isapprox(ha.Psa, hb.Psa; atol=0.01)
    @test isapprox(ha.QL, hb.QL; rtol=1e-3)
end

@testset "Fisiologia: hemorragia, barorreflexo e terapias" begin
    hem = window_rate(1500.0, 60.0, 600.0)
    run(; baro=true, qinf=0.0, hct=0.0, dose=0.0) = begin
        sc = Scenario(Qhem=hem, Qinf=window_rate(qinf, 1200.0, 900.0), hct_inf=hct,
                      dose=t -> (1200.0 <= t < 5000.0 ? dose : 0.0), baro=baro,
                      tstops=[60.0, 660.0, 1200.0, 2100.0, 5000.0])
        sol = simulate_mean(P0, sc, (0.0, 5400.0))
        return mean_observables(sol, P0, sc, [1100.0, 5000.0])
    end
    sem = run(baro=false); com = run(baro=true)
    @test sem.PAM[1] < 50                      # sem compensação: colapso
    @test com.PAM[1] > sem.PAM[1] + 15         # o barorreflexo sustenta a pressão
    @test com.HR[1] > P0.HR + 20               # taquicardia compensatória
    rl = run(qinf=1500.0, hct=0.0); sg = run(qinf=1000.0, hct=0.40); nt = run()
    @test rl.PAM[2] > nt.PAM[2] && sg.PAM[2] > nt.PAM[2]    # ambos elevam a PAM
    @test rl.Hct[2] < sg.Hct[2] < 0.43                        # RL dilui; sangue preserva o Hct
    @test sg.DO2[2] > rl.DO2[2]                               # DO₂ maior com sangue total
    @test rl.Vis[2] > sg.Vis[2]                               # RL vaza mais para o interstício
    @test run(dose=0.1).PAM[2] > nt.PAM[2]                    # vasopressor eleva a PAM
end

@testset "Dados: leitura do MIMIC (esquema), CSV e paciente sintético" begin
    t0 = DateTime(2150, 1, 1, 8, 0, 0)
    at(s) = t0 + Second(s)
    ids = MIMIC_ITEMIDS
    chart = DataFrame(stay_id = fill(7, 8),
        charttime = [at(0), at(60), at(120), at(180), at(0), at(120), at(240), at(300)],
        itemid = [fill(ids.chart.MAP_invasive, 4); fill(ids.chart.HR, 4)],
        valuenum = [60.0, 62.0, 70.0, 75.0, 110.0, 105.0, 100.0, 98.0])
    inputs = DataFrame(stay_id = [7, 7, 7],
        starttime = [at(60), at(60), at(120)], endtime = [at(180), at(180), at(240)],
        itemid = [ids.input.ringer, ids.input.norepinephrine, ids.input.prbc],
        amount = [1000.0, 0.0, 300.0], rate = [0.0, 0.08, 0.0])
    labs = DataFrame(charttime = [at(0), at(150)], itemid = [ids.lab.hct, ids.lab.hct],
                     valuenum = [40.0, 33.0])
    pd = build_patient(7, chart, inputs, labs; dt=60.0)
    @test pd.t == [0.0, 60.0, 120.0, 180.0]
    @test pd.PAM == [60.0, 62.0, 70.0, 75.0]
    @test pd.HR[1:3] == [110.0, 110.0, 105.0]              # forward-fill
    @test isapprox(sum(pd.rl_mlh) * 60 / 3600, 1000.0; atol=1e-9)    # volume conservado
    @test isapprox(sum(pd.blood_mlh) * 60 / 3600, 300.0; atol=1e-9)      # PRBC [120,240) cabe na grade
    @test pd.dose == [0.0, 0.08, 0.08, 0.0]
    @test pd.Hct[1] == 0.40 && pd.Hct[4] == 0.33
    sc = scenario_from_patient(pd)
    @test isapprox(sc.Qinf(90.0), 1000.0 / 120, atol=1e-9)       # mL/s
    @test sc.hr_forced(0.0) == 110.0

    # ida e volta pelo CSV processado
    tmp = joinpath(mktempdir(), "p.csv")
    write_patient_csv(tmp, pd)
    pd2 = load_patient(tmp)
    @test pd2.t == pd.t && pd2.PAM == pd.PAM && pd2.dose == pd.dose

    # paciente sintético é reprodutível e a calibração recupera a verdade
    ptrue = with_params(P0; Rsa = 1.25 * P0.Rsa, Vb0 = 0.88 * P0.Vb_ref)
    a, _ = synthetic_patient(ptrue; hours=1.5, dt=120.0, noise=1.5, seed=5,
                             rl = t -> (900 <= t < 1800 ? 1000.0 : 0.0), blood = t -> 0.0,
                             dose = t -> 0.0)
    b, _ = synthetic_patient(ptrue; hours=1.5, dt=120.0, noise=1.5, seed=5,
                             rl = t -> (900 <= t < 1800 ? 1000.0 : 0.0), blood = t -> 0.0,
                             dose = t -> 0.0)
    @test a.PAM == b.PAM
    r0 = rmse(simulate_patient(P0, a), a.PAM)
    fit = calibrate(P0, a, [:Rsa, :Vb0]; maxiters=120)
    @test fit.rmse < r0
    @test fit.rmse < 2.5
    # janela curta só com cristaloide: Rsa e Vb0 se compensam (identificabilidade parcial);
    # o volume inicial é razoavelmente determinado, a resistência só grosseiramente.
    # Com os parâmetros da B.2 (Csv = 160, mais volume estressado) a PAM é menos sensível ao
    # volume e a identificabilidade piora: antes 5 %/25 %, agora 15 %/50 % (ver docs/RELATORIO.md).
    @test isapprox(fit.θ[2], ptrue.Vb0; rtol=0.15)
    @test isapprox(fit.θ[1], ptrue.Rsa; rtol=0.50)
end
