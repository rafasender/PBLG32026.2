# scripts/plot_style.jl
# Estilo comum das figuras (incluído por run_model.jl, figuras_mimic.jl e prever_paciente.jl):
# cores fixas por tratamento, cabeçalho explicativo, faixas de fase com nome e rótulos diretos.

using Plots
gr()

"Versão para o relatório (`FIG_RELATORIO=1`): sem o cabeçalho de texto, que vira a legenda no documento."
const RELATORIO = get(ENV, "FIG_RELATORIO", "0") == "1"
"Pasta de saída das figuras (`docs/figuras/` ou `docs/figuras/relatorio/`)."
pasta_figuras(root) = (d = joinpath(root, "docs", "figuras", RELATORIO ? "relatorio" : ""); mkpath(d); d)

# tinta do texto, cores das séries (ordem fixa, validada para daltonismo) e fundos das fases
const TINTA = (forte="#0b0b0b", media="#52514e", fraca="#898781", grade="#e1e0d9")
const SERIE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
const FUNDO = (perda="#f6dcd8", tratamento="#d9e7f8", ajuste="#ecebe6", teste="#ffffff")
# o mesmo tratamento tem sempre a mesma cor, em todas as figuras
const COR_TRAT = Dict("Sem tratamento" => TINTA.fraca, "Sem terapia" => TINTA.fraca,
                      "Cristaloide 1000 mL" => SERIE[1], "Ringer Lactato 1500 mL" => SERIE[1],
                      "Hemácias 600 mL" => SERIE[2], "Sangue total 1000 mL" => SERIE[2],
                      "Noradrenalina 0,1" => SERIE[3], "Noradrenalina 0,1 µg/kg/min" => SERIE[3],
                      "RL 1000 mL + noradrenalina" => SERIE[4])

default(fontfamily="sans-serif", linewidth=2, framestyle=:axes, grid=true, gridalpha=0.35,
        gridcolor=TINTA.grade, foreground_color_axis=TINTA.fraca, foreground_color_border=TINTA.fraca,
        tickfontcolor=TINTA.media, guidefontcolor=TINTA.media, titlefontcolor=TINTA.forte,
        titlefontsize=11, guidefontsize=9, tickfontsize=8, legendfontsize=8,
        legend_foreground_color=TINTA.grade, background_color=:white, dpi=110)

"Quebra `s` em linhas de no máximo `cols` caracteres, entre palavras."
function quebra(s, cols)
    linhas = String[]; atual = ""
    for w in split(s, ' ')
        if !isempty(atual) && length(atual) + 1 + length(w) > cols
            push!(linhas, atual); atual = String(w)
        else
            atual = isempty(atual) ? String(w) : atual * " " * w
        end
    end
    isempty(atual) || push!(linhas, atual)
    return linhas
end

"""
    cabecalho(titulo, paragrafos; cols=150)

Painel só de texto para o topo da figura: título e parágrafos de explicação ("o que foi
usado", "como ler"), quebrados em linhas de até `cols` caracteres. Vai no layout com uma altura
relativa pequena.
"""
function cabecalho(titulo, paragrafos; cols=150)
    linhas = reduce(vcat, [quebra(s, cols) for s in paragrafos])
    p = plot(framestyle=:none, grid=false, legend=false, xlims=(0, 1), ylims=(0, 1))
    n = length(linhas) + 1
    annotate!(p, 0.0, 1 - 0.5 / n, text(titulo, :left, 13, TINTA.forte))
    for (i, l) in enumerate(linhas)
        annotate!(p, 0.0, 1 - (i + 0.5) / n, text(l, :left, 9, TINTA.media))
    end
    return p
end

"Altura em pixels do cabeçalho com estes parágrafos (≈ 24 px por linha + título)."
altura_cabecalho(paragrafos; cols=150) = 24 * (sum(length(quebra(s, cols)) for s in paragrafos) + 2)

"Faixa de fundo [a, b] (desenhar ANTES das séries)."
fase!(p, a, b, cor) = vspan!(p, [a, b], color=cor, linealpha=0, label="")

"""
    nomes_fases!(p, fases)

Escreve o nome de cada fase `(a, b, nome)` no topo do painel, abrindo espaço acima das
séries (chamar DEPOIS das séries).
"""
function nomes_fases!(p, fases; tam=8)
    y1, y2 = ylims(p)
    topo = y2 + 0.16 * (y2 - y1)
    ylims!(p, (y1, topo))
    for (a, b, nome) in fases
        annotate!(p, (a + b) / 2, topo - 0.05 * (topo - y1), text(nome, :center, tam, TINTA.media))
    end
    return p
end

"Rótulo direto no fim de uma linha (texto em tinta neutra; a cor fica só na linha)."
rotulo!(p, x, y, txt; tam=8) = annotate!(p, x, y, text("  " * txt, :left, tam, TINTA.forte))
