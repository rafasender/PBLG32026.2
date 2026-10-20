# scripts/demo_windkessel2.jl
# Simulação rápida do modelo Windkessel 2 elementos com DifferentialEquations.jl

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))
using Cardio0D
using DifferentialEquations
using Plots

println("=================================================================")
println(" Simulação Windkessel 2 Elementos - Modelo Cardiovascular 0D")
println("=================================================================")

# Parâmetros fisiológicos baseline
# FC = 75 bpm -> Período cardíaco T = 60 / 75 = 0.8 s
# Duração da sístole Tsys ≈ 0.3 s
HR = 75.0
T = 60.0 / HR
T_sys = 0.30
SV = 70.0 # Volume sistólico ~70 mL

# Função de fluxo de ejeção aórtica Qin(t) periódica
function Qin(t)
    t_cycle = mod(t, T)
    if t_cycle < T_sys
        # Pulso de ejeção ventricular aproximado por meia senoide
        return (SV * pi / (2 * T_sys)) * sin(pi * t_cycle / T_sys)
    else
        return 0.0
    end
end

# 1. Cenário Baseline (Normovolemia)
params_base = Windkessel2Params(R=1.0, C=1.2, Pv=4.0)
tspan = (0.0, 5.0 * T) # 5 ciclos cardíacos
P0 = [80.0] # Pressão arterial inicial em mmHg

prob_base = ODEProblem(windkessel2_ode!, P0, tspan, (params_base, Qin))
sol_base = solve(prob_base, Tsit5(), reltol=1e-6, abstol=1e-6)

# Último ciclo cardíaco para métricas de estado estacionário
t_eval = range(4.0 * T, 5.0 * T, length=200)
p_eval = [sol_base(t)[1] for t in t_eval]

P_sys = maximum(p_eval)
P_dia = minimum(p_eval)
PAM_integral = calculate_map(collect(t_eval), p_eval)
PAM_approx = P_dia + (P_sys - P_dia) / 3.0

println("\n--- Cenário 1: Linha de Base (Normovolemia) ---")
println("  Pressão Sistólica (PAS) : $(round(P_sys, digits=1)) mmHg")
println("  Pressão Diastólica (PAD): $(round(P_dia, digits=1)) mmHg")
println("  PAM (via integral)      : $(round(PAM_integral, digits=1)) mmHg")
println("  PAM (fórmula clínica)   : $(round(PAM_approx, digits=1)) mmHg")

# 2. Cenário Choque Hemorrágico (Queda no débito cardíaco / SV reduzido para 35 mL e reflexo simpático elevando R e HR)
HR_shock = 110.0
T_shock = 60.0 / HR_shock
T_sys_shock = 0.25
SV_shock = 35.0

function Qin_shock(t)
    t_cycle = mod(t, T_shock)
    if t_cycle < T_sys_shock
        return (SV_shock * pi / (2 * T_sys_shock)) * sin(pi * t_cycle / T_sys_shock)
    else
        return 0.0
    end
end

# Compensação simpática inicial eleva RVS (R)
params_shock = Windkessel2Params(R=1.5, C=1.0, Pv=2.0)
prob_shock = ODEProblem(windkessel2_ode!, [60.0], (0.0, 5.0 * T_shock), (params_shock, Qin_shock))
sol_shock = solve(prob_shock, Tsit5(), reltol=1e-6, abstol=1e-6)

t_eval_shock = range(4.0 * T_shock, 5.0 * T_shock, length=200)
p_eval_shock = [sol_shock(t)[1] for t in t_eval_shock]

P_sys_shock = maximum(p_eval_shock)
P_dia_shock = minimum(p_eval_shock)
PAM_shock = calculate_map(collect(t_eval_shock), p_eval_shock)

println("\n--- Cenário 2: Choque Hemorrágico com Resposta Simpática ---")
println("  Pressão Sistólica (PAS) : $(round(P_sys_shock, digits=1)) mmHg")
println("  Pressão Diastólica (PAD): $(round(P_dia_shock, digits=1)) mmHg")
println("  PAM calculada           : $(round(PAM_shock, digits=1)) mmHg")

# Gerar gráfico comparativo
mkpath(normpath(joinpath(@__DIR__, "..", "output")))
plot_path = normpath(joinpath(@__DIR__, "..", "output", "demo_windkessel2.png"))

p = plot(sol_base, vars=(0, 1), label="Normovolemia (Baseline)", lw=2, color=:blue,
         xlabel="Tempo (s)", ylabel="Pressão Arterial (mmHg)",
         title="Modelo 0D Windkessel: Baseline vs Choque Hipovolêmico")
plot!(p, sol_shock, vars=(0, 1), label="Choque Hipovolêmico compensado", lw=2, color=:red, linestyle=:dash)
hline!(p, [65.0], label="Limiar crítico PAM (65 mmHg)", color=:gray, linestyle=:dot, lw=1.5)
savefig(p, plot_path)

println("\n[SUCESSO] Gráfico da simulação salvo em: output/demo_windkessel2.png")
