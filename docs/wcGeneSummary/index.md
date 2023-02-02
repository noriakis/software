--- 
title: "wcGeneSummary (OSplot)"
author: "Noriaki Sato"
date: "2023-02-02"
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib]
url: https://noriakis.github.io/software/wcGeneSummary
# cover-image: path to the social sharing image like images/cover.jpg
biblio-style: apalike
csl: chicago-fullnote-bibliography.csl
---

# About

This package aims to visualize the word and text information contained in the gene or the other omics identifiers such as microbiome, and identify important words among the clusters, and compare the clusters based on those information. It contributes to understanding the gene clusters and aid in easy interpretation and visualization. The documentation using bookdown is available [here](https://noriakis.github.io/software/wcGeneSummary), and the web server using `shinyapps.io` is [here](https://nsato.shinyapps.io/osplotweb/).



```r
library(wcGeneSummary)
#> 
#> Registered S3 method overwritten by 'pvclust':
#>   method       from      
#>   text.pvclust dendextend
knitr::include_url("https://noriakis.github.io/cyjs_test/")
```

<iframe src="https://noriakis.github.io/cyjs_test/" width="672" height="400px" data-external="1"></iframe>