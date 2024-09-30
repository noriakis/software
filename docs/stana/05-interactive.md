


# Interactive inspection

`stana` offers the interactive inspection of metagenotyping results for examining the intra-species diversity.


``` r
library(stana)
load(system.file("extdata", "sysdata.rda", package = "stana"))
```

For this purpose, `exportInteractive()` function is prepared. This expects gene family (KO) abundance tables set to the `kos` slot, and the tree is set to `treeList` slot by using the functions like `consensusSeq`. However, if they are not available, the gene copy number table alone can be used to export the interactive application. The function outputs the necessary data and shiny codes to the specified directory and run the app.
Specifying `notRun` to TRUE will not run the application.


``` r
## Default to export to the current directory
exportInteractive(stana, notRun=TRUE)
#> # No tree for 100003
#> # Tree number: 0 KO (or gene) number: 1
#> # Exporting ...
#> # A stana: MIDAS2
#> # Database: uhgg
#> # Loaded directory: midas2_sample_merge_uhgg
#> # Species number: 1
#> # Group info (list): Group1/Group2
#> # Group column (DF): label/group
#> # Loaded KO table: 1 ID: 100003
#> # Size: 0.882912 MB
```
