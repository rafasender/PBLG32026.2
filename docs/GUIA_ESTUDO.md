# Guia de estudo do projeto

Versão curta. A matemática completa está no [CLAUDE.md](../CLAUDE.md), para quem quiser se aprofundar.

---

## 1. O projeto em 5 frases

1. O programa trata a circulação como água num circuito: o coração é a bomba, e artérias, veias e pulmão são os canos e reservatórios.
2. Ele simula uma **hemorragia** (o circuito perde líquido) e mostra a pressão (PAM) e o débito caindo.
3. Ele compara **tratamentos**: cristaloide (Ringer, soro 0,9%), concentrado de hemácias e noradrenalina.
4. O modelo é **ajustado a pacientes reais** do MIMIC-IV com choque hipovolêmico, para ficar parecido com eles.
5. Com o modelo ajustado, ele **prevê** qual tratamento faria cada paciente responder melhor.

## 2. Como o programa funciona: 3 etapas

```
 1. MODELO                  2. CALIBRAR (MIMIC)                     3. PREVER
 equações da circulação  →  ajusta 3 números a cada paciente     →  simula cada tratamento
 com valores de livro       com as primeiras 6 h de dados,          sozinho e escolhe o melhor
 (Frente B)                 e testa nas 6 h seguintes               (maior O₂ com PAM ≥ 65)
```

- **Os 3 números ajustados por paciente:** resistência dos vasos, volume de sangue no início e frequência cardíaca "de base". Todo o resto é fixo, com valores de literatura iguais para todos. Isso evita *overfitting*.
- **O teste:** depois de ajustado nas 6 primeiras horas, o modelo tem que acertar as 6 horas seguintes sozinho. Comparamos com um palpite bobo ("a PAM vai ficar igual à última medida"). Se o modelo não ganhar do palpite, ele não serve para prever.
- **O critério de "melhor tratamento":** entre os tratamentos que mantêm a PAM ≥ 65 mmHg, o que entrega mais oxigênio aos tecidos (DO₂). Pressão alta sozinha não basta.

## 3. Como rodar

```bash
julia --project=. test/runtests.jl              # confere se o código está certo (~3 min)
julia --project=. scripts/run_model.jl          # figuras 1–5 (simulações, ~1 min)
julia --project=. scripts/pipeline_mimic.jl     # MIMIC: calibra, testa e prevê; figuras 7–8 (~1h30)
julia --project=. scripts/prever_paciente.jl    # pergunta PAM, FC e Hct e prevê os 4 tratamentos
```

O último é o que dá para mostrar ao vivo na apresentação. Ele pergunta só **PAM, FC e hematócrito**, porque são os únicos dados que o modelo usa. Idade, sexo e peso não mudam nenhum parâmetro do modelo (ele é de um adulto de 70 kg), então perguntar seria enganoso.

## 4. As figuras

Cada figura tem, no topo, um texto dizendo **de onde vieram os dados** e **como ler**. As faixas coloridas de fundo marcam as etapas: vermelho = perdendo sangue, azul = tratamento, cinza = período de ajuste. As cores dos tratamentos são as mesmas em todas as figuras.

As figuras 1–5 usam um **paciente hipotético padrão** (adulto de 70 kg, valores normais de literatura) em cenários escolhidos por nós. Elas não provam que o modelo prevê: mostram que ele **se comporta como a fisiologia conhecida** (verificação). Quem mostra a previsão são as figuras 7 e 8.

| Figura | O que mostra | Dados usados | Mensagem principal | Onde usar |
|---|---|---|---|---|
| 1 | Coração batendo: pressões e o ciclo pressão–volume (PV loop) | Só o modelo (valores de livro) | O coração do modelo funciona: PA 97/76, débito 4,9 L/min, volume ejetado 65 mL | Relatório (verificação) |
| 2 | Curvas de Guyton | Só o modelo | Perder 400 mL derruba o débito de 5,25 para 4,08 L/min; repor volume corrige | Relatório (verificação) |
| 3 | Hemorragia de 1500 mL (32% do sangue, início da classe III) sem tratamento | Só o modelo | O corpo compensa: a FC sobe a 113 e a PAM para de cair em 74 | Relatório (verificação) |
| 4 | Hemorragias de 600 a 2000 mL | Só o modelo | Quanto mais perda, menor a PAM. **Limitação:** 2000 mL não colapsa | Relatório (verificação, como limitação) |
| 5 | Os tratamentos depois de uma hemorragia simulada | Só o modelo | Noradrenalina sobe a PAM rápido, mas entrega menos O₂; sangue entrega mais O₂; Ringer dilui | Relatório; na apresentação, só como slide de reserva |
| 7 | O modelo ajustado a pacientes reais | MIMIC-IV | O modelo acompanha a PAM real? Ganha do palpite "fica igual"? | **Apresentação** (2–3 pacientes) e relatório |
| 8 | Previsão do melhor tratamento por paciente | MIMIC-IV + modelo | Qual tratamento o modelo recomenda para cada um | **Apresentação** e relatório |

A figura 6 não existe mais (era a calibração com dados sintéticos, removida).

**Na apresentação:** figuras 7 e 8 e a demonstração ao vivo do `prever_paciente.jl`. Guardem a figura 5 como slide de reserva para a pergunta "como sabem que o modelo é fisiologicamente plausível?".

**Como o código gera a figura 7:** (1) seleciona os pacientes; (2) monta uma tabela minuto a minuto com PAM, FC, Hct e os tratamentos que entraram; (3) o otimizador testa ~100 combinações dos 3 números, simulando as primeiras 6 h com os tratamentos reais, e fica com a de menor erro contra a PAM e a FC medidas; (4) com esses números simula as 12 h (linha azul) e mede o erro nas 6 h finais, que o modelo não viu.

**Como gera a figura 8:** com os mesmos 3 números, simula o paciente 4 vezes a partir do mesmo estado inicial, cada vez com um único tratamento no lugar dos reais, e escolhe o melhor pelo critério.

## 5. O barorreflexo, em português claro

- **O que é:** o reflexo que segura a pressão. Quando a PAM cai, o corpo acelera o coração e contrai os vasos.
- **Por que está no modelo:** para **prever**. Quando o modelo simula "e se esse paciente recebesse noradrenalina?", não existe FC medida desse cenário, porque ele nunca aconteceu. O modelo precisa calcular sozinho como o corpo reagiria.
- **E os dados do MIMIC já não têm o reflexo?** Têm, e por isso a FC medida **não é copiada**: ela é usada para conferir se o reflexo do modelo se comporta como o do paciente.
- **O "barorreflexo OFF" que aparecia antes** era só uma demonstração didática e foi removido. O reflexo fica sempre ligado.

## 6. Resultados com pacientes reais

Detalhes em `docs/resultados_mimic.md`. Em resumo:

- **Coorte:** 956 internações de UTI com CID de choque hipovolêmico ou traumático (sem sepse); 936 tratadas; 220 com dados suficientes; 10 calibradas.
- **Calibrar ajuda:** o erro no teste caiu de 22,5 mmHg (modelo genérico) para 15,2 mmHg (mediana), e melhorou em 9 de 10 pacientes.
- **Mas não ganha do palpite ingênuo:** o modelo errou menos que "a PAM fica igual" em só **4 de 10** pacientes (mediana 15,2 contra 11,0 mmHg). Com dados de hora em hora, repetir a última medida é uma previsão difícil de bater.
- **Parâmetros nas bordas:** com limites fisiológicos, o volume ficou no mínimo (60%) em 6 pacientes e a resistência no mínimo em 5. O modelo só explica PAMs tão baixas forçando os extremos, provavelmente porque o reflexo dele é de pessoa saudável e não há sangramento contínuo.
- **Sem os limites** (primeira tentativa), o modelo ganhava em 5 de 10, mas com pacientes impossíveis (0,25 L de sangue, FC de 22 bpm): sobreajuste. Por isso os limites ficaram.
- **Previsão:** hemácias foram o melhor tratamento nos 10 pacientes. A noradrenalina dá a maior PAM em todos, mas reduz a oferta de O₂ em 7. O resultado uniforme vem do critério (O₂ depende da hemoglobina), não de diferenças entre pacientes.
- **Conclusão honesta:** o modelo serve para explorar mecanismos e comparar tratamentos qualitativamente, mas ainda **não** para prever um paciente individual.

## 7. Limitações (o que dizer se perguntarem)

- **Modelo 0D:** representa médias por compartimento. Não mostra a onda de pulso detalhada, e a pressão de pulso fica baixa.
- **Reflexo de pessoa saudável:** ele compensa mais do que o de um paciente sedado na UTI. Por isso 2000 mL de perda não colapsa (figura 4).
- **MIMIC:** sinais vitais a cada ~1 h (poucos pontos por paciente); o diagnóstico (CID) é da alta, então diz *quem* teve choque, não *quando*; "sangue total" não existe na base (usamos concentrado de hemácias); só a parte 1 do `chartevents` estava disponível.
- **A previsão de um tratamento que o paciente não recebeu não pode ser conferida com dados.** A confiança vem do teste da figura 7.
- **Uso educacional:** não serve para decisão clínica.

## 8. Perguntas prováveis (resposta curta)

- **Por que um modelo 0D?** Queremos pressão e débito ao longo de horas, não a forma da onda. 0D é rápido o bastante para ajustar a cada paciente.
- **O que é calibrar?** Ajustar poucos números do modelo para ele imitar um paciente, mantendo o resto fixo.
- **Como sabem que não é *overfitting*?** Ajustamos nas 6 primeiras horas e testamos nas 6 seguintes, que o modelo não viu, comparando com um palpite ingênuo.
- **Por que a noradrenalina "perde" se ela sobe a pressão?** Ela contrai os vasos: a PAM sobe, mas o débito cai, e chega menos oxigênio aos tecidos.
- **Por que não usaram sangue total?** O MIMIC-IV não tem registros de sangue total, só de concentrado de hemácias.
- **Por que não perguntam idade e sexo?** O modelo não tem parâmetros que dependam deles.

## 9. Onde está cada coisa no código

| Arquivo | Para que serve |
|---|---|
| `src/params.jl` | Os valores de livro (Frente B) e o paciente "padrão" |
| `src/physiology.jl` | Coração, válvulas, reflexo, noradrenalina, hematócrito, troca de líquido |
| `src/mean_model.jl` | O modelo "lento" (médias por batimento), usado para horas de UTI |
| `src/closed_loop.jl` | O modelo batimento a batimento (figura 1) |
| `src/data_io.jl` | Leitura do MIMIC, seleção por CID, grupos de tratamento |
| `src/calibration.jl` | Ajuste dos 3 números por paciente |
| `src/treatments.jl` | Os 4 tratamentos comparados e o critério de "melhor" |
| `scripts/` | Os programas que você roda (seção 3) e o estilo dos gráficos |
| `test/runtests.jl` | Testes automáticos que conferem as equações |

**Dados do MIMIC nunca vão para o GitHub** (termo de uso da PhysioNet). Só resultados agregados e anônimos (P01, P02…).
