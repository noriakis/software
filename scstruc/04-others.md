# Other algorithms and software




## Bootstrapping

For obtaining the reliable network, we can use bootstrapping approach by randomly subsampling the data, infer the network from it, and average the bootstrapped networks. Most algorithms can use bootstrapping by specifying `boot=TRUE`. For instance, to use bootstrapping for CCDr algorithm (which is fast), the example codes are shown below. Replication number can be specified by `R` and sampling number can be specified by `m` argument. In this way, `bn.strength` object will be returned as `net`.


``` r
library(scstruc)

sce <- mockSCE()
sce <- logNormCounts(sce)
included_genes <- sample(row.names(sce), 40)
gs <- scstruc(sce, included_genes, R=30,
              changeSymbol=FALSE, algorithm="ccdr", boot=TRUE)
#> Bootstrapping specified
#> Returning the list of bn.strength

## It consists of strength from multiple lambdas
head(gs$net[[5]])
#>        from        to strength direction
#> 1 Gene_0268 Gene_0332        0         0
#> 2 Gene_0268 Gene_0399        0         0
#> 3 Gene_0268 Gene_0413        0         0
#> 4 Gene_0268 Gene_0518        0         0
#> 5 Gene_0268 Gene_0566        0         0
#> 6 Gene_0268 Gene_0596        0         0
```

## `PIDC`

PIDC algorithm is based on information theory and leverages partial information decomposition (PID) to analyze and quantify the relationships between multiple variables [@chan_gene_2017]. `scstruc` uses the undirected network produced from `PIDC`, implemented in Julia. `JuliaCall` is used to call Julia from R, and the users should first set up the Julia environment. Also, the path to the main library ([NetworkInference.jl](https://github.com/Tchanders/NetworkInference.jl)) should be specified to `NetworkInference_HOME` argument. The directed acyclic graph is then obtained based on the score maximization with the constraints produced by `PIDC`. The example code is shown as follows:



``` r
library(JuliaCall)
julia <- julia_setup(JULIA_HOME = "./Julia-1.10.5/bin")
gs <- scstruc(sce, included_genes,
              changeSymbol=FALSE,
              algorithm="pidc",
              algorithm.args=list(NetworkInference_HOME="./NetworkInference.jl"))
gs$net
```

The bootstrapped-based inference can be also performed. In this case, the averaged network is obtained per `thresholds` parameters in `algorithm.args`. By default, the network is thresholded by `seq(0.1, 0.4, 0.1)`.


``` r
gs <- scstruc(sce, included_genes,
              changeSymbol=FALSE,
              algorithm="pidc", boot=TRUE, R=30,
              algorithm.args=list(NetworkInference_HOME="./NetworkInference.jl"))
gs$net
```

The maximization can be chosen from the scoring function available in `bnlearn` or Greedy Equivalence Search (`ges`), which can be specified by the parameter `maximize` in `algorithm.args`.

Note that raw network output that is not scaled and thresholded is stored in the temporary directory.

## `GENIE3`

GENIE3 is a widely used method for inferring directed interactions between genes based on an ensemble of regression trees. By specifying `genie3` in `algorithm`, the function performs `GENIE3` algorithm and subsequently threshold the arcs. Then, the function returns the list of the original results and the networks if the thresholded network is DAG.


``` r
gs <- scstruc(sce, included_genes,
              changeSymbol=FALSE,
              algorithm="genie3")
head(gs$net$original)
#>            Gene_0268  Gene_0332  Gene_0399  Gene_0413
#> Gene_0268 0.00000000 0.04806280 0.03070923 0.02950172
#> Gene_0332 0.03812961 0.00000000 0.02597613 0.02584433
#> Gene_0399 0.02731085 0.02948143 0.00000000 0.02595806
#> Gene_0413 0.02788775 0.02794079 0.02462729 0.00000000
#> Gene_0518 0.02786071 0.04559585 0.04388638 0.02711829
#> Gene_0566 0.03400589 0.03555564 0.02683549 0.02723738
#>            Gene_0518  Gene_0566  Gene_0596  Gene_0598
#> Gene_0268 0.02667846 0.03775870 0.03543492 0.02975625
#> Gene_0332 0.04566691 0.03769726 0.06203825 0.04807715
#> Gene_0399 0.03792178 0.03876648 0.04997654 0.03264433
#> Gene_0413 0.02799680 0.02985924 0.03269626 0.03735218
#> Gene_0518 0.00000000 0.03087539 0.03659631 0.03943075
#> Gene_0566 0.03633175 0.00000000 0.02488180 0.03061393
#>            Gene_0634  Gene_0649  Gene_0671  Gene_0696
#> Gene_0268 0.03041270 0.02941533 0.02677473 0.03643082
#> Gene_0332 0.02959120 0.02812892 0.02661490 0.02698970
#> Gene_0399 0.04353496 0.03329067 0.02547989 0.03879671
#> Gene_0413 0.02922098 0.03568671 0.02853269 0.02694774
#> Gene_0518 0.03136156 0.02140210 0.03969633 0.02308859
#> Gene_0566 0.03063246 0.02983986 0.03281514 0.03084788
#>            Gene_0728  Gene_0737  Gene_0770  Gene_0795
#> Gene_0268 0.03066103 0.03773379 0.02844980 0.02407317
#> Gene_0332 0.03210181 0.03371274 0.02408368 0.02536610
#> Gene_0399 0.03298882 0.02472669 0.03035161 0.02891186
#> Gene_0413 0.03720801 0.03131826 0.03122379 0.03656113
#> Gene_0518 0.02566118 0.05191985 0.02593010 0.07423479
#> Gene_0566 0.03009736 0.02629928 0.02884896 0.03102382
#>            Gene_0805  Gene_0868  Gene_0942  Gene_0945
#> Gene_0268 0.04687927 0.03797992 0.02472566 0.03571946
#> Gene_0332 0.02658287 0.02870913 0.02361332 0.02629898
#> Gene_0399 0.03164461 0.02950199 0.04783306 0.03902401
#> Gene_0413 0.02725425 0.02521442 0.02443385 0.03326018
#> Gene_0518 0.02541573 0.02351137 0.02189042 0.03392706
#> Gene_0566 0.02916401 0.04044313 0.03458140 0.03063736
#>            Gene_0985  Gene_0995  Gene_1021  Gene_1085
#> Gene_0268 0.03051196 0.02847305 0.03275989 0.03003508
#> Gene_0332 0.03308601 0.03090671 0.02534707 0.03064466
#> Gene_0399 0.03071515 0.03178637 0.05252687 0.03137185
#> Gene_0413 0.02803579 0.05766842 0.02553859 0.03277104
#> Gene_0518 0.02504285 0.04407958 0.02352920 0.03660910
#> Gene_0566 0.03226340 0.03033738 0.02638007 0.03131746
#>            Gene_1091  Gene_1152  Gene_1212  Gene_1303
#> Gene_0268 0.02888809 0.02762018 0.07199049 0.04025052
#> Gene_0332 0.03638740 0.02670459 0.02347369 0.02341440
#> Gene_0399 0.03672325 0.02925462 0.03168339 0.02866118
#> Gene_0413 0.03474642 0.03967002 0.02910728 0.03573482
#> Gene_0518 0.02096500 0.02933942 0.02573424 0.02449710
#> Gene_0566 0.02848721 0.02956891 0.02894976 0.02607299
#>            Gene_1343  Gene_1404  Gene_1478  Gene_1493
#> Gene_0268 0.03408796 0.02783098 0.03319151 0.02119563
#> Gene_0332 0.02858623 0.02250102 0.04324024 0.02208566
#> Gene_0399 0.03249124 0.02671115 0.03251746 0.03009837
#> Gene_0413 0.04176923 0.03338838 0.02985457 0.03398381
#> Gene_0518 0.02491081 0.05369118 0.03898393 0.02688708
#> Gene_0566 0.02576682 0.02669661 0.04282265 0.02293499
#>            Gene_1504  Gene_1663  Gene_1701  Gene_1716
#> Gene_0268 0.02620566 0.03591880 0.03497526 0.04206038
#> Gene_0332 0.02059514 0.02900346 0.02181436 0.03954160
#> Gene_0399 0.02555442 0.03104142 0.02362603 0.03148938
#> Gene_0413 0.03201972 0.02857039 0.02521726 0.02862532
#> Gene_0518 0.04494758 0.02299887 0.02005112 0.02731326
#> Gene_0566 0.02265762 0.02990485 0.03274941 0.03050448
#>            Gene_1761  Gene_1771  Gene_1926  Gene_1970
#> Gene_0268 0.03795837 0.02962692 0.02862905 0.02906427
#> Gene_0332 0.02384148 0.02770915 0.02613935 0.02820980
#> Gene_0399 0.03462029 0.02537008 0.04733737 0.02784980
#> Gene_0413 0.03817140 0.03378951 0.03289098 0.02656692
#> Gene_0518 0.01999025 0.03072816 0.02638434 0.03204576
#> Gene_0566 0.05504639 0.04605090 0.03006919 0.03185736
```
