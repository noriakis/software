--- 
title: "scstruc: causal assessment of gene regulatory network using Bayesian network"
author: "Noriaki Sato"
date: "2025-04-15"
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib, packages.bib, article.bib]
# url: your book url like https://bookdown.org/yihui/bookdown
# cover-image: path to the social sharing image like images/cover.jpg
description: |
  This is the documentation of the package scstruc.
link-citations: yes
github-repo: rstudio/bookdown-demo
---




# Introduction

`scstruc` is a package designed to estimate gene regulatory networks (GRNs) from single-cell transcriptomics (SCT) data using Bayesian networks (BN). Notably, the package focuses on following points:

- Implementation of algorithms tailored to SCT data characteristics, such as dropouts (or zero-inflation).
- Evaluating the inferred networks based on their causal validity.
- Probabilistic reasoning and the extraction of differences in regulatory relationships between groups
- Comprehensive visualization for understanding the resulting GRNs. 

BN estimate directed relationships between genes and are useful for GRN inference. The package also supports DAG estimation from the other popular GRN inference software tools, significantly expanding the applicability of GRNs in SCT data analysis.

## Installation and prerequiresties


``` r
devtools::install_github("noriakis/scstruc")
```

We need some packages that needs to be installed before using the full functions of scstruc.

CCDr algorithm from the following repository:


``` r
devtools::install_github("itsrainingdata/sparsebnUtils")
devtools::install_github("noriakis/ccdrAlgorithm")
```

HurdleNormal package from the following repository:


``` r
devtools::install_github("amcdavid/HurdleNormal")

```

For CCDr algorithm and `HurdleNormal` package, please consult the original papers [@mcdavid_graphical_2019; @aragam_concave_2015].

