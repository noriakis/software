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
#> 1 Gene_0006 Gene_0024        0         0
#> 2 Gene_0006 Gene_0035        0         0
#> 3 Gene_0006 Gene_0068        0         0
#> 4 Gene_0006 Gene_0155        0         0
#> 5 Gene_0006 Gene_0178        0         0
#> 6 Gene_0006 Gene_0268        0         0
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
#>            Gene_0006  Gene_0024  Gene_0035  Gene_0068
#> Gene_0006 0.00000000 0.03088093 0.03834214 0.02447807
#> Gene_0024 0.01957930 0.00000000 0.01893201 0.01754363
#> Gene_0035 0.02940798 0.02949616 0.00000000 0.02877513
#> Gene_0068 0.02461398 0.02816775 0.02747613 0.00000000
#> Gene_0155 0.02534968 0.01813408 0.02943684 0.03197248
#> Gene_0178 0.01897838 0.01895927 0.01880557 0.02050140
#>            Gene_0155  Gene_0178  Gene_0268  Gene_0285
#> Gene_0006 0.02904142 0.02439686 0.03957938 0.02679946
#> Gene_0024 0.01628740 0.01759800 0.01877264 0.02292221
#> Gene_0035 0.02826699 0.02700899 0.02617555 0.02605871
#> Gene_0068 0.03897850 0.03458415 0.03934045 0.03165985
#> Gene_0155 0.00000000 0.01974402 0.01963389 0.01967546
#> Gene_0178 0.01788149 0.00000000 0.02186202 0.03218175
#>            Gene_0354  Gene_0407  Gene_0461  Gene_0474
#> Gene_0006 0.02850355 0.02917745 0.03266967 0.02571844
#> Gene_0024 0.01606783 0.02395562 0.01966547 0.01893470
#> Gene_0035 0.02683567 0.03307719 0.02994596 0.02386858
#> Gene_0068 0.03237185 0.02726279 0.02734305 0.02415926
#> Gene_0155 0.01717159 0.02254696 0.02825549 0.02070662
#> Gene_0178 0.02039013 0.01681829 0.01839037 0.01885531
#>            Gene_0613  Gene_0621  Gene_0680  Gene_0681
#> Gene_0006 0.03358618 0.02437391 0.02878870 0.03107345
#> Gene_0024 0.03415715 0.01830281 0.01928706 0.01511642
#> Gene_0035 0.02688471 0.02964478 0.03129199 0.02744632
#> Gene_0068 0.02968370 0.02863495 0.02710226 0.02435848
#> Gene_0155 0.01912834 0.02435755 0.02494557 0.02167478
#> Gene_0178 0.02741718 0.02005268 0.02154768 0.01997393
#>            Gene_0690  Gene_0693  Gene_0724  Gene_0754
#> Gene_0006 0.02117753 0.02516211 0.02302914 0.02903713
#> Gene_0024 0.01918560 0.01560746 0.01935694 0.01796692
#> Gene_0035 0.07474305 0.02469208 0.02533006 0.02589232
#> Gene_0068 0.02973997 0.03769346 0.03891703 0.03198661
#> Gene_0155 0.02300340 0.02003054 0.01969889 0.02790590
#> Gene_0178 0.02342628 0.03023963 0.02036297 0.01985533
#>            Gene_0784  Gene_0836  Gene_0846  Gene_0947
#> Gene_0006 0.02654190 0.02522234 0.03080537 0.04013842
#> Gene_0024 0.02381901 0.04265656 0.03259398 0.02239880
#> Gene_0035 0.02410693 0.02934640 0.02363304 0.03234239
#> Gene_0068 0.02552641 0.02492546 0.02285037 0.03852429
#> Gene_0155 0.02996549 0.01904311 0.02677543 0.01825622
#> Gene_0178 0.01990014 0.02781590 0.01843462 0.02738928
#>            Gene_0975  Gene_1023  Gene_1063  Gene_1140
#> Gene_0006 0.03074997 0.03291853 0.02754159 0.02421725
#> Gene_0024 0.01867381 0.01925868 0.01577452 0.01421171
#> Gene_0035 0.02611413 0.02420067 0.02620930 0.02706129
#> Gene_0068 0.02724259 0.02900529 0.02391209 0.03909870
#> Gene_0155 0.02572670 0.02679155 0.02425918 0.02579700
#> Gene_0178 0.02148480 0.02041351 0.02616502 0.01902780
#>            Gene_1242  Gene_1299  Gene_1365  Gene_1374
#> Gene_0006 0.02584367 0.02797067 0.02987608 0.04080190
#> Gene_0024 0.01991268 0.01750479 0.02221940 0.01767164
#> Gene_0035 0.03177332 0.02396212 0.02861492 0.04031653
#> Gene_0068 0.03291799 0.02486527 0.03217252 0.02250598
#> Gene_0155 0.02607348 0.02427945 0.03154361 0.02210219
#> Gene_0178 0.02074743 0.02089802 0.01930644 0.02431743
#>            Gene_1446  Gene_1464  Gene_1574  Gene_1576
#> Gene_0006 0.03338292 0.02854316 0.03283427 0.03580522
#> Gene_0024 0.02074742 0.01822139 0.02073954 0.03753362
#> Gene_0035 0.03286167 0.03049468 0.02893392 0.02882686
#> Gene_0068 0.02353721 0.02729962 0.03067578 0.03270469
#> Gene_0155 0.02419758 0.02007164 0.02684538 0.04662410
#> Gene_0178 0.02050590 0.02056526 0.01820344 0.01708360
#>            Gene_1594  Gene_1631  Gene_1946  Gene_1976
#> Gene_0006 0.03441058 0.03193840 0.03100051 0.04085575
#> Gene_0024 0.01465738 0.01546324 0.01749811 0.01991336
#> Gene_0035 0.05522610 0.02543445 0.04158306 0.03529794
#> Gene_0068 0.02501158 0.02039993 0.02746443 0.02887373
#> Gene_0155 0.01917789 0.01540677 0.02308983 0.03455652
#> Gene_0178 0.02527642 0.01895051 0.04094081 0.02435783
```
