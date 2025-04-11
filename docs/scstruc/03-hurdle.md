# Hurdle model



For handling the zero-inflated nature of the single-cell transcriptomics data, the use of the hurdle model is presented. The hurdle model is a two-part model that models whether the observation is zero or not, and non-zero part separately.

Based on the inferred network in `HurdleNormal`, which proposed the multivariate Hurdle model and grouped lasso to learn the undirected network [@mcdavid_graphical_2019], the directed acyclic graph is inferred based on score maximization using the constraints. The undirected network is selected from multiple lambdas based on the BIC criteria implemented in `HurdleNormal`. For this function, `Hurdle` algorithm should be specified in `algorithm` argument in `scstruc`. By default, score will be `BIC` in `bnlearn`.


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
#>    [Gene_0013][Gene_0072][Gene_0347][Gene_0377][Gene_0516]
#>    [Gene_0535][Gene_0633][Gene_0638][Gene_0708][Gene_0858]
#>    [Gene_0860][Gene_0914][Gene_0932][Gene_0951][Gene_1079]
#>    [Gene_1268][Gene_1586][Gene_1655][Gene_1905][Gene_1972]
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
#>    [Gene_0013][Gene_0072][Gene_0347][Gene_0377][Gene_0516]
#>    [Gene_0535][Gene_0633][Gene_0638][Gene_0708][Gene_0858]
#>    [Gene_0860][Gene_0914][Gene_0932][Gene_0951][Gene_1079]
#>    [Gene_1268][Gene_1586][Gene_1655][Gene_1905][Gene_1972]
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
## This performs GES
gs.ges <- scstruc(sce, included_genes,
    changeSymbol=FALSE, algorithm="Hurdle",
    algorithm.args=list(maximize="ges"))
gs.ges$net
#> 
#>   Random/Generated Bayesian network
#> 
#>   model:
#>    [Gene_0013][Gene_0072][Gene_0347][Gene_0377][Gene_0516]
#>    [Gene_0535][Gene_0633][Gene_0638][Gene_0708][Gene_0858]
#>    [Gene_0860][Gene_0914][Gene_0932][Gene_0951][Gene_1079]
#>    [Gene_1268][Gene_1586][Gene_1655][Gene_1905][Gene_1972]
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
#>    [Gene_0013][Gene_0072][Gene_0347][Gene_0377][Gene_0516]
#>    [Gene_0535][Gene_0633][Gene_0638][Gene_0708][Gene_0858]
#>    [Gene_0860][Gene_0914][Gene_0932][Gene_0951][Gene_1079]
#>    [Gene_1268][Gene_1586][Gene_1655][Gene_1905][Gene_1972]
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
