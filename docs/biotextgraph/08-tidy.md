

# Perform each process separately {#tidy}

The same process as functions like `refseq` can be performed separately. This is useful for piped processing and for individual customization in each process.


```r
library(biotextgraph)
btg <- obtain_refseq(c("DDX41","PNKP","IRF3")) |>
  set_filter_words() |>
  make_corpus() |>
  make_TDM() |>
  make_graph() |>
  graph_cluster() |>
  process_network_gene(gene_plot=TRUE, gene_path_plot="reactome") |>
  plot_biotextgraph(edge_link=FALSE) |>
  plot_wordcloud()
#> Input genes: 3
#>   Converted input genes: 3
#> Filter based on GeneSummary
#> Filtered 76 words (frequency and/or tfidf)
#> Found 25 enriched term
btg
#> Type: refseq
#> Number of words: 30
#> Query: DDX41/PNKP/IRF3
#> Graph: V(33), E(237)
#> Degree: response(32)/immune(23)/innate(23)/addition(16)/alteration(16)
#> 301.4 Kb

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
#> Filtered 76 words (frequency and/or tfidf)
btg2
#> Type: enrich
#> Number of words: 30
#> Query: DDX41/PNKP/IRF3
#> Graph: V(29), E(44)
#> Degree: innate(8)/immune(7)/responses(7)/activatesmodulates(6)/NA(6)
#> 292 Kb
```


## Use of gene descriptions in Alliance of Genome Resources

We can use gene descriptions of alliance of genome resources, that can be obtained from [here](https://www.alliancegenome.org/downloads). Download the file and call `obtain_alliance` function. The default path is set to the working directory and human gene description.


```r
btg_agr <- obtain_alliance(c("DDX41","PNKP","IRF3")) |>
  make_corpus() |> ## No filter this time
  make_TDM() |>
  make_graph() |>
  graph_cluster() |>
  process_network_gene(gene_plot=TRUE, gene_path_plot="reactome") |>
  plot_biotextgraph(edge_link=FALSE)
#> Input genes: 3
#> Found 25 enriched term
```
