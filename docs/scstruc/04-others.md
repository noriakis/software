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
#> 1 Gene_0136 Gene_0302        0         0
#> 2 Gene_0136 Gene_0307        0         0
#> 3 Gene_0136 Gene_0350        0         0
#> 4 Gene_0136 Gene_0392        0         0
#> 5 Gene_0136 Gene_0401        0         0
#> 6 Gene_0136 Gene_0457        0         0
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
#>             Gene_0136  Gene_0302  Gene_0307   Gene_0350
#> Gene_0136 0.000000000 0.03815198 0.03033649 0.041082345
#> Gene_0302 0.038278390 0.00000000 0.04384009 0.027663166
#> Gene_0307 0.030111652 0.03673630 0.00000000 0.032976816
#> Gene_0350 0.016747356 0.01259568 0.01475499 0.000000000
#> Gene_0392 0.009417026 0.00905313 0.02542493 0.007675829
#> Gene_0401 0.026928721 0.03057284 0.02406323 0.023513165
#>            Gene_0392  Gene_0401   Gene_0457   Gene_0604
#> Gene_0136 0.04648226 0.02337528 0.026669667 0.030160030
#> Gene_0302 0.02192453 0.02674055 0.029739473 0.032546259
#> Gene_0307 0.06882408 0.02316867 0.028509166 0.036082195
#> Gene_0350 0.01194514 0.01299722 0.013846909 0.017583805
#> Gene_0392 0.00000000 0.01883923 0.009379975 0.009928166
#> Gene_0401 0.03699605 0.00000000 0.026166402 0.025136633
#>             Gene_0713   Gene_0750  Gene_0864  Gene_0896
#> Gene_0136 0.025855796 0.021457870 0.02811444 0.02466951
#> Gene_0302 0.024273916 0.023555467 0.03804105 0.03145132
#> Gene_0307 0.046605300 0.025364098 0.02651592 0.02677376
#> Gene_0350 0.012690999 0.012314518 0.01280650 0.01707095
#> Gene_0392 0.009130004 0.006187762 0.01058934 0.00804439
#> Gene_0401 0.030871800 0.032225826 0.02977111 0.03693539
#>            Gene_0900  Gene_0952   Gene_1019   Gene_1027
#> Gene_0136 0.02664938 0.02794521 0.050803463 0.025559635
#> Gene_0302 0.02333425 0.02743567 0.028156731 0.022775713
#> Gene_0307 0.04213985 0.04642870 0.033456730 0.030720022
#> Gene_0350 0.02097608 0.03077974 0.013197604 0.009960530
#> Gene_0392 0.01831480 0.00720375 0.009492629 0.005857273
#> Gene_0401 0.02755639 0.02627090 0.023445935 0.027167016
#>             Gene_1044   Gene_1099   Gene_1141   Gene_1147
#> Gene_0136 0.030426695 0.031965744 0.019138151 0.030471918
#> Gene_0302 0.029496119 0.046808942 0.013459680 0.028850265
#> Gene_0307 0.030023075 0.053413942 0.053407348 0.036376766
#> Gene_0350 0.014516524 0.015075026 0.007550333 0.013949234
#> Gene_0392 0.008448556 0.006533706 0.005959179 0.008713158
#> Gene_0401 0.048214400 0.027910790 0.035217399 0.046985978
#>             Gene_1212   Gene_1270   Gene_1313   Gene_1357
#> Gene_0136 0.022230224 0.028804101 0.022075456 0.021049147
#> Gene_0302 0.029052077 0.023988628 0.040953583 0.015107704
#> Gene_0307 0.046269274 0.040661734 0.025518922 0.027737605
#> Gene_0350 0.014653098 0.017480964 0.008155586 0.011108790
#> Gene_0392 0.005272744 0.007265429 0.024178523 0.009322583
#> Gene_0401 0.041266070 0.025051974 0.024298585 0.020992240
#>             Gene_1362   Gene_1448  Gene_1449   Gene_1456
#> Gene_0136 0.024161605 0.044487650 0.03141692 0.027353634
#> Gene_0302 0.036407367 0.027787208 0.03280477 0.035540784
#> Gene_0307 0.033298731 0.042767262 0.03369469 0.032973421
#> Gene_0350 0.029097022 0.012382960 0.01976695 0.016965853
#> Gene_0392 0.005977574 0.009708553 0.01339975 0.006919641
#> Gene_0401 0.031397431 0.045816305 0.04198381 0.051922271
#>             Gene_1459  Gene_1509   Gene_1521   Gene_1677
#> Gene_0136 0.023175863 0.03086472 0.028955898 0.030377296
#> Gene_0302 0.029205700 0.03027946 0.045265261 0.024402367
#> Gene_0307 0.033317710 0.04266091 0.032629080 0.038115352
#> Gene_0350 0.020894528 0.04888214 0.016485845 0.012038425
#> Gene_0392 0.008546227 0.01131372 0.009456666 0.007903083
#> Gene_0401 0.029065808 0.03615510 0.022199438 0.027399003
#>             Gene_1743   Gene_1766   Gene_1822  Gene_1841
#> Gene_0136 0.026207363 0.037843780 0.028362633 0.04885591
#> Gene_0302 0.022705805 0.031373235 0.028968008 0.02407048
#> Gene_0307 0.037131385 0.032466744 0.027821938 0.03447045
#> Gene_0350 0.031977465 0.014364427 0.016248574 0.02234902
#> Gene_0392 0.008937382 0.008547491 0.006495196 0.01031908
#> Gene_0401 0.071085403 0.028437538 0.038390071 0.03079519
#>             Gene_1885   Gene_1893  Gene_1903  Gene_1906
#> Gene_0136 0.028342254 0.030577202 0.02564471 0.03109586
#> Gene_0302 0.032120416 0.043366644 0.02661660 0.03050130
#> Gene_0307 0.032491832 0.027815304 0.05801220 0.03257176
#> Gene_0350 0.017024285 0.029956807 0.01982418 0.01370361
#> Gene_0392 0.006693124 0.009723063 0.01181725 0.00840537
#> Gene_0401 0.044197370 0.030558183 0.02218325 0.03685088
```
