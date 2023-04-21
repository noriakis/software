

# Perform each process separately {#tidy}

The same process as functions like `wcGeneSummary` can be performed separately. This is useful for piped processing and for individual customization in each process.


```r
library(biotextgraph)
btg <- obtain_refseq(c("DDX41","PNKP","IRF3")) |>
  set_filter_words() |>
  make_corpus() |>
  make_TDM() |>
  make_graph() |>
  process_network_gene(gene_plot=TRUE, gene_path_plot="reactome") |>
  plot_biotextgraph(edge_link=FALSE) |>
  plot_wordcloud()
#> Input genes: 3
#>   Converted input genes: 3
#> Filter based on GeneSummary
#> Filtered 81 words (frequency and/or tfidf)
#> Found 21 enriched term
btg
#> Type: refseq
#> Number of words: 30
#> Query: DDX41/PNKP/IRF3
#> Graph: V(33), E(237)
#> Degree: response(32)/immune(23)/innate(23)/addition(16)/alteration(16)
#> 291 Kb

## Text of enrichment analysis results
btg2 <- obtain_enrich(c("DDX41","PNKP","IRF3"), enrich="reactome") |>
  set_filter_words() |>
  make_corpus() |>
  make_TDM() |>
  make_graph() |>
  process_network_manual() |>
  plot_biotextgraph(edge_link=FALSE)
#> Input genes: 3
#>   Converted input genes: 3
#> Performing enrichment analysis
#> Filter based on GeneSummary
#> Filtered 81 words (frequency and/or tfidf)
btg2
#> Type: enrich
#> Number of words: 30
#> Query: DDX41/PNKP/IRF3
#> Graph: V(30), E(46)
#> Degree: immune(7)/innate(7)/responses(7)/activatesmodulates(5)/adaptive(5)
#> 291.3 Kb
```