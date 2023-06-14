--- 
title: "ggkegg"
author: "Noriaki Sato"
date: "2023-06-14"
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib]
# url: your book url like https://bookdown.org/yihui/bookdown
# cover-image: path to the social sharing image like images/cover.jpg
biblio-style: apalike
csl: chicago-fullnote-bibliography.csl
---



# About

`ggkegg` fetches information from KEGG and parse, analyze and visualize them using `ggplot2` and `ggraph`, combined with the other packages investigating biological functions using KEGG. This package aims to visualize the complex components of KEGG using the grammar of graphics. For Python, please use [`pykegg`](https://pykegg.readthedocs.io) using `plotnine`, which offers the almost same functionality as `ggkegg` used in conjuction with the package such as `gseapy` and `PyDESeq2` and single-cell transcriptomics analysis library `scanpy`.


```r
# devtools::install_github("noriakis/ggkegg")
library(ggkegg)
```

One of the main aims of `ggkegg` or `tidykegg` is manupilating KEGG information in tidy ways using `tidygraph`, and offers the customized visualization of KEGG information including KEGG PATHWAY, MODULE, and NETWORK.


```r
library(dplyr)
library(tidygraph)
pathway("hsa04110") |> ## Obtain and parse the KEGG pathway
  activate(nodes) |> ## node manipulation
  mutate(convert_hsa=convert_id("hsa"),
         convert_map=convert_id("pathway")) |> ## convert IDs for organism hsa and pathway
  ggraph(x=x, y=y)+ ## ggraph plot
  geom_edge_parallel(arrow=arrow(length=unit(1,"mm")),
                     aes(linetype=subtype_name),
                     end_cap=circle(7.5,"mm"))+ ## Parallel edges
  geom_node_rect(aes(filter=type=="gene",
                     fill=I(bgcolor)),
                 color="black")+ ## rectangular nodes
  geom_node_text(aes(label=convert_hsa),
                 size=2, family="serif")+ ## text
  geom_node_text(aes(label=convert_hsa,
                     filter=!is.na(convert_hsa) & convert_hsa=="TP53"),
                 size=2, color="red", family="serif")+ ## highlight text
  theme_void()
```

<img src="index_files/figure-html/ggkegg-1.png" width="100%" style="display: block; margin: auto;" />
