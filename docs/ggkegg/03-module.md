

# Module

Module information can be obtained and parsed. Parsing `DEFINITION` and `REACTION` is supported. For the definition, first the function breaks down the definition to steps, and make graphical representation using `ggraph` or text itself.


```r
library(ggkegg)
mod <- obtain_module("M00004")
mod
#> $name
#> [1] "Pentose phosphate pathway (Pentose phosphate cycle)"
#> 
#> $definition
#> [1] "(K13937,((K00036,K19243) (K01057,K07404))) K00033 K01783 (K01807,K01808) K00615 ((K00616 (K01810,K06859,K15916)),K13810)"
#> 
#> $reaction
#>  [1] "R02736,R10907  C01172 -> C01236"           
#>  [2] "R02035  C01236 -> C00345"                  
#>  [3] "R01528,R10221  C00345 -> C00199"           
#>  [4] "R01529  C00199 -> C00231"                  
#>  [5] "R01056  C00199 -> C00117"                  
#>  [6] "R01641  C00117 + C00231 -> C05382 + C00118"
#>  [7] "R01827  C05382 + C00118 -> C00279 + C05345"
#>  [8] "R01830  C00279 + C00231 -> C05345 + C00118"
#>  [9] "R02740  C05345 -> C00668"                  
#> [10] "R02739  C00668 -> C01172"
```
