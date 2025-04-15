# Hurdle model



For handling the zero-inflated nature of the single-cell transcriptomics data, the use of the hurdle model is presented. The hurdle model is a two-part model that models whether the observation is zero or not, and non-zero part separately. For this function, `Hurdle` algorithm should be specified in `algorithm` argument in `scstruc`.

Based on the inferred network in `HurdleNormal`, which proposed the multivariate Hurdle model and grouped lasso to learn the undirected network [@mcdavid_graphical_2019], the directed acyclic graph is inferred based on score maximization using the constraints. By default, score will be `BIC` in `bnlearn`. The customized score is also introduced in this section. The first undirected network is selected from multiple lambdas based on the BIC criteria implemented in `HurdleNormal`.



``` r
library(HurdleNormal)
library(scstruc)

sce <- mockSCE()
sce <- logNormCounts(sce)
included_genes <- sample(row.names(sce), 20)
gs <- scstruc(sce, included_genes,
	changeSymbol=FALSE, algorithm="Hurdle")
gs$net
#> 
#>   Bayesian network learned via Score-based methods
#> 
#>   model:
#>    [Gene_0182][Gene_0419][Gene_0475][Gene_0708][Gene_0742]
#>    [Gene_0780][Gene_0795][Gene_0848][Gene_0921][Gene_0940]
#>    [Gene_1003][Gene_1110][Gene_1194][Gene_1291][Gene_1555]
#>    [Gene_1557][Gene_1609][Gene_1663][Gene_1808][Gene_1844]
#>   nodes:                                 20 
#>   arcs:                                  0 
#>     undirected arcs:                     0 
#>     directed arcs:                       0 
#>   average markov blanket size:           0.00 
#>   average neighbourhood size:            0.00 
#>   average branching factor:              0.00 
#> 
#>   learning algorithm:                    Hill-Climbing 
#>   score:                                 BIC (Gauss.) 
#>   penalization coefficient:              2.649159 
#>   tests used in the learning procedure:  0 
#>   optimized:                             TRUE
```

The score maximization function can be set arbitrarily (`maximizeFun`), set to `hc` by default. Greedy Equivalence Search can be performed via setting `maximize` argument to `ges`.


``` r
gs.tabu <- scstruc(sce, included_genes,
    changeSymbol=FALSE, algorithm="Hurdle",
    algorithm.args=list(maximizeFun=bnlearn::tabu))
gs.tabu$net
#> 
#>   Bayesian network learned via Score-based methods
#> 
#>   model:
#>    [Gene_0182][Gene_0419][Gene_0475][Gene_0708][Gene_0742]
#>    [Gene_0780][Gene_0795][Gene_0848][Gene_0921][Gene_0940]
#>    [Gene_1003][Gene_1110][Gene_1194][Gene_1291][Gene_1555]
#>    [Gene_1557][Gene_1609][Gene_1663][Gene_1808][Gene_1844]
#>   nodes:                                 20 
#>   arcs:                                  0 
#>     undirected arcs:                     0 
#>     directed arcs:                       0 
#>   average markov blanket size:           0.00 
#>   average neighbourhood size:            0.00 
#>   average branching factor:              0.00 
#> 
#>   learning algorithm:                    Tabu Search 
#>   score:                                 BIC (Gauss.) 
#>   penalization coefficient:              2.649159 
#>   tests used in the learning procedure:  0 
#>   optimized:                             TRUE
```


``` r
## This performs GES as score-based learning
gs.ges <- scstruc(sce, included_genes,
    changeSymbol=FALSE, algorithm="Hurdle",
    algorithm.args=list(maximize="ges"))
gs.ges$net
#> 
#>   Random/Generated Bayesian network
#> 
#>   model:
#>    [Gene_0182][Gene_0419][Gene_0475][Gene_0708][Gene_0742]
#>    [Gene_0780][Gene_0795][Gene_0848][Gene_0921][Gene_0940]
#>    [Gene_1003][Gene_1110][Gene_1194][Gene_1291][Gene_1555]
#>    [Gene_1557][Gene_1609][Gene_1663][Gene_1808][Gene_1844]
#>   nodes:                                 20 
#>   arcs:                                  0 
#>     undirected arcs:                     0 
#>     directed arcs:                       0 
#>   average markov blanket size:           0.00 
#>   average neighbourhood size:            0.00 
#>   average branching factor:              0.00 
#> 
#>   generation algorithm:                  Empty
```


## Customized score for hurdle model

Additional score can be specified by using `hurdle.bic` function in `algorithm.args` argument. The score is defined as the sum of BIC values from continuous and logistic regression part of the hurdle model. Let $Y = [y_{ij}]$ denote the log-normalized expression value of the gene i in subset cell j and $Z = [z_{ij}]$ a 0-1 indicator value of whether the gene expression is zero. The score is defined as:

$$
\text{logit}\left(\Pr(Z_{ij} = 1)\right) \sim \beta_i G_j,
$$

$$
\Pr(Y_{ij} = y \mid Z_{ij} = 1) \sim N(\beta_i G_j, \sigma^2).
$$



``` r
gs2 <- scstruc(sce, included_genes,
	changeSymbol=FALSE, algorithm="Hurdle",
	algorithm.args=list("score"=hurdle.bic))
gs2$net
#> 
#>   Bayesian network learned via Score-based methods
#> 
#>   model:
#>    [Gene_0182][Gene_0419][Gene_0475][Gene_0708][Gene_0742]
#>    [Gene_0780][Gene_0795][Gene_0848][Gene_0921][Gene_0940]
#>    [Gene_1003][Gene_1110][Gene_1194][Gene_1291][Gene_1555]
#>    [Gene_1557][Gene_1609][Gene_1663][Gene_1808][Gene_1844]
#>   nodes:                                 20 
#>   arcs:                                  0 
#>     undirected arcs:                     0 
#>     directed arcs:                       0 
#>   average markov blanket size:           0.00 
#>   average neighbourhood size:            0.00 
#>   average branching factor:              0.00 
#> 
#>   learning algorithm:                    Hill-Climbing 
#>   score:                                 
#>                                  User-Provided Score Function 
#>   tests used in the learning procedure:  0 
#>   optimized:                             TRUE
```

As proposed in the `MAST` library, the cellular detection rate adjustment (CDR) can be performed in the scoring phase by `cdrAdjuetment` to `TRUE`. This applies inclusion of CDR term in the hurdle modeling and score maximizing phase.

## `add.dropout`

Like in the other SCT data simulator, the excessive zero in the matrix can be simulated by `add.dropout`. The function takes absolute expression values if any of the expression in the matrix is negative.


``` r
net <- readRDS("../arth150.rds")
sim <- rbn(net, 100)
table(sim == 0)
#> 
#> FALSE 
#> 10700
sim.do <- sim * add.dropout(sim)
table(sim.do == 0)
#> 
#> FALSE  TRUE 
#>  4361  6339
```



