# scripts/calibrate_synthetic.jl
# Ensaio completo de calibração com um paciente SINTÉTICO (verdade conhecida):
#   dado -> CSV processado -> leitura -> calibração -> sensibilidade -> figura.
# Valida o pipeline antes de usar dados reais do MIMIC-IV (Frente C).
#
#   julia --project=. scripts/calibrate_synthetic.jl

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))

using Cardio0D
using Plots
using Printf

gr()
const ROOT = normpath(joinpath(@__DIR__, ".."))
mkpath(joinpath(ROOT, "docs", "figuras"))
mkpath(joinpath(ROOT, "data", "processed"))
const OUT = IOBuffer()
out(args...) = (println(OUT, args...); println(args...))

p0 = baseline_params()
# paciente "real": TPR 25 % maior, contratilidade 15 % menor, déficit de volume de 14 %
ptrue = with_params(p0; Rsa=1.25 * p0.Rsa, Vb0=0.86 * p0.Vb_ref, EmaxL=0.85 * p0.EmaxL)

# 1) gera e grava no formato processado; 2) relê do disco (o mesmo caminho do MIMIC)
pd_gen, _ = synthetic_patient(ptrue; id="sintetico_01", hours=3.0, dt=60.0, noise=2.0, seed=3)
csv_path = joinpath(ROOT, "data", "processed", "sintetico_01.csv")   # ignorado pelo git (*.csv)
write_patient_csv(csv_path, pd_gen)
pd = load_patient(csv_path)

out("## Paciente sintético (3 h, PAM a cada 60 s, ruído gaussiano de 2 mmHg)\n")
out("Verdade: `Rsa = $(round(ptrue.Rsa, digits=3))`, `Vb0 = $(round(ptrue.Vb0, digits=0))` mL, ",
    "`EmaxL = $(round(ptrue.EmaxL, digits=3))`. Infusões: RL 1500 mL (30–50 min), noradrenalina ",
    "0,1 µg/kg/min (60–120 min), sangue 800 mL (90–110 min).\n")
r0 = rmse(simulate_patient(p0, pd), pd.PAM)
out(@sprintf("RMSE do modelo basal sem calibrar: **%.2f mmHg**.\n", r0))

# sensibilidade local (identificabilidade prática)
names3 = [:Rsa, :Vb0, :EmaxL]
S, c = sensitivity_matrix(p0, pd, names3)
out("### Sensibilidade local normalizada (média de |S| na janela)\n")
out("| Parâmetro | média \\|S\\| |\n|---|---|")
for (j, n) in enumerate(names3)
    out(@sprintf("| %s | %.3f |", n, sum(abs, S[:, j]) / size(S, 1)))
end
out(@sprintf("\ncond(SᵀS) = %.2e  (a condição sozinha não denuncia o problema: o que aponta o parâmetro não identificável é a sensibilidade média muito menor, aqui `EmaxL`).\n", c))

# calibrações
out("### Resultados da calibração\n")
out("| Parâmetros ajustados | λ | RMSE (mmHg) | estimativa | verdade | erro relativo |\n|---|---|---|---|---|---|")
fits = Dict()
for (label, ns, λ) in (("Rsa, Vb0", [:Rsa, :Vb0], 0.0),
                       ("Rsa, Vb0, EmaxL", names3, 0.0),
                       ("Rsa, Vb0, EmaxL", names3, 0.5))
    r = calibrate(p0, pd, ns; λ, maxiters=250)
    truth = [getfield(ptrue, n) for n in ns]
    est = join([@sprintf("%.3g", x) for x in r.θ], " / ")
    tru = join([@sprintf("%.3g", x) for x in truth], " / ")
    err = join([@sprintf("%+.1f%%", 100 * (a - b) / b) for (a, b) in zip(r.θ, truth)], " / ")
    out(@sprintf("| %s | %.1f | %.2f | %s | %s | %s |", label, λ, r.rmse, est, tru, err))
    fits[(label, λ)] = r
end
out("")

best = fits[("Rsa, Vb0", 0.0)]
y_sim = simulate_patient(best.p, pd)
y_unc = simulate_patient(p0, pd)
fig = plot(pd.t ./ 60, pd.PAM, seriestype=:scatter, ms=2.5, label="PAM 'medida' (sintética)",
           xlabel="Tempo (min)", ylabel="PAM (mmHg)", title="Calibração em paciente sintético",
           color=:gray, size=(950, 480))
plot!(fig, pd.t ./ 60, y_unc, label="Modelo basal (sem calibrar)", color=:orange, linestyle=:dash)
plot!(fig, pd.t ./ 60, y_sim, label="Modelo calibrado (Rsa, Vb0)", color=:blue, lw=2.5)
savefig(fig, joinpath(ROOT, "docs", "figuras", "06_calibracao_sintetica.png"))

write(joinpath(ROOT, "docs", "resultados_calibracao.md"), String(take!(OUT)))
println("\n[OK] docs/resultados_calibracao.md e docs/figuras/06_calibracao_sintetica.png")
