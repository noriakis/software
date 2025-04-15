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
#> 1 Gene_0156 Gene_0229        0         0
#> 2 Gene_0156 Gene_0231        0         0
#> 3 Gene_0156 Gene_0323        0         0
#> 4 Gene_0156 Gene_0336        0         0
#> 5 Gene_0156 Gene_0350        0         0
#> 6 Gene_0156 Gene_0376        0         0
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
#>            Gene_0156  Gene_0229  Gene_0231   Gene_0323
#> Gene_0156 0.00000000 0.03028585 0.03391460 0.028297127
#> Gene_0229 0.03485903 0.00000000 0.02685464 0.023199691
#> Gene_0231 0.02504192 0.02325473 0.00000000 0.018069640
#> Gene_0323 0.01330440 0.01289305 0.01165647 0.000000000
#> Gene_0336 0.01387871 0.01000074 0.01248417 0.009755696
#> Gene_0350 0.03365878 0.02816917 0.03527866 0.024367847
#>             Gene_0336   Gene_0350   Gene_0376  Gene_0444
#> Gene_0156 0.070875001 0.035628644 0.034072305 0.03610615
#> Gene_0229 0.025028168 0.030756747 0.032017148 0.03047891
#> Gene_0231 0.022989495 0.026319057 0.016737978 0.03054703
#> Gene_0323 0.009869578 0.013729302 0.021490612 0.01224105
#> Gene_0336 0.000000000 0.007330711 0.008323914 0.01051796
#> Gene_0350 0.020088861 0.000000000 0.030250537 0.02957796
#>             Gene_0487  Gene_0804  Gene_0846  Gene_0856
#> Gene_0156 0.027183744 0.02686804 0.03104218 0.02501029
#> Gene_0229 0.033594750 0.02750657 0.02897075 0.04739677
#> Gene_0231 0.019500411 0.02036525 0.03566455 0.01975604
#> Gene_0323 0.013909248 0.01418789 0.02271942 0.02982123
#> Gene_0336 0.007838251 0.01182224 0.01161484 0.01398596
#> Gene_0350 0.033563407 0.03037086 0.02856881 0.03013948
#>            Gene_0955  Gene_1024  Gene_1025  Gene_1095
#> Gene_0156 0.04415820 0.03707054 0.02987893 0.03411332
#> Gene_0229 0.02551375 0.02526397 0.02492996 0.03685898
#> Gene_0231 0.01765163 0.05755167 0.02365340 0.02120446
#> Gene_0323 0.00924323 0.01134262 0.01520570 0.03066596
#> Gene_0336 0.01456278 0.01820598 0.01109332 0.01147988
#> Gene_0350 0.03159843 0.03257884 0.02872809 0.03191452
#>             Gene_1117   Gene_1153  Gene_1186  Gene_1202
#> Gene_0156 0.033519299 0.023575979 0.03278377 0.04166837
#> Gene_0229 0.023120704 0.034732414 0.03311580 0.03056318
#> Gene_0231 0.015526007 0.020404076 0.02300868 0.02740384
#> Gene_0323 0.008892518 0.012305456 0.01849592 0.01489978
#> Gene_0336 0.023706508 0.008458542 0.01176431 0.01090599
#> Gene_0350 0.048786133 0.033617082 0.03125257 0.03529774
#>             Gene_1257  Gene_1290   Gene_1297   Gene_1321
#> Gene_0156 0.031691666 0.02550102 0.044217887 0.026045983
#> Gene_0229 0.029435506 0.02901794 0.032874545 0.029900634
#> Gene_0231 0.021768284 0.02010626 0.020335334 0.026729682
#> Gene_0323 0.017240952 0.01289683 0.011303454 0.012148152
#> Gene_0336 0.009658607 0.02648709 0.008851073 0.008587002
#> Gene_0350 0.028932234 0.02582708 0.028632917 0.023811084
#>            Gene_1372  Gene_1457  Gene_1471  Gene_1537
#> Gene_0156 0.02784701 0.02988726 0.02845569 0.02707732
#> Gene_0229 0.02780712 0.03744497 0.02609542 0.03712001
#> Gene_0231 0.02405152 0.01902755 0.02167770 0.02410179
#> Gene_0323 0.01398449 0.03161001 0.01170533 0.01381456
#> Gene_0336 0.01299827 0.02251519 0.01067377 0.01350765
#> Gene_0350 0.02822658 0.02237824 0.03072959 0.02867282
#>             Gene_1545  Gene_1593  Gene_1670   Gene_1727
#> Gene_0156 0.033948275 0.02985123 0.02447357 0.041362762
#> Gene_0229 0.030545531 0.02088252 0.03696399 0.040073972
#> Gene_0231 0.029501297 0.01641544 0.02041282 0.017931735
#> Gene_0323 0.011327183 0.01587072 0.01330209 0.011527915
#> Gene_0336 0.009613473 0.01149877 0.01696952 0.009199474
#> Gene_0350 0.030764519 0.03571193 0.03207682 0.040063411
#>            Gene_1849  Gene_1905  Gene_1921  Gene_1938
#> Gene_0156 0.03306036 0.03233509 0.03281073 0.02866882
#> Gene_0229 0.03090804 0.02514331 0.03191349 0.03733548
#> Gene_0231 0.01896607 0.03050193 0.02650402 0.01851383
#> Gene_0323 0.01234698 0.02280683 0.02842592 0.01311272
#> Gene_0336 0.02067467 0.01005650 0.02645692 0.01462112
#> Gene_0350 0.03954748 0.03086294 0.02116200 0.02998779
#>            Gene_1966   Gene_1981   Gene_1983  Gene_1997
#> Gene_0156 0.02514435 0.037165455 0.024047831 0.01921417
#> Gene_0229 0.02649093 0.035184183 0.045474490 0.03456340
#> Gene_0231 0.01920766 0.020855721 0.022368184 0.01654087
#> Gene_0323 0.01745645 0.009293087 0.013671040 0.01523694
#> Gene_0336 0.01356252 0.014113312 0.007202735 0.00902023
#> Gene_0350 0.02427172 0.031094622 0.029315544 0.03629895
```
