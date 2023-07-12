# Importing

`stana` is aimed to import the metagenotyping results of various pipelines.
The package is designed primarily for `MIDAS` and `MIDAS2`, which outputs the 
gene abundances by default.


```r
library(stana)
```

## MIDAS

For `MIDAS`, `loadMIDAS` function can be used to import the output of `merge` command.

## MIDAS2

For `MIDAS2`, `loadMIDAS2` function can be used to import the output of `merge` command.

## inStrain

For inStrain, we need `compare` command output, after `profile` command. First we would like to know which species are profiled. `just_species` option only list species in the table. Here, we load an example dataset profiled by `inStrain`, using the default database described in their original tutorial. For `inStrain`, loading all the SNV tables is often impossible, we must specify `candidate_species` for investigation.


```r
sp <- loadInStrain("../inStrain_out", just_species=TRUE)
sp
#>  [1] "GUT_GENOME000022" "GUT_GENOME000024"
#>  [3] "GUT_GENOME000147" "GUT_GENOME000220"
#>  [5] "GUT_GENOME000221" "GUT_GENOME000224"
#>  [7] "GUT_GENOME000225" "GUT_GENOME000231"
#>  [9] "GUT_GENOME000509" "GUT_GENOME007566"
#> [11] "GUT_GENOME009103" "GUT_GENOME031782"
#> [13] "GUT_GENOME034989" "GUT_GENOME044231"
#> [15] "GUT_GENOME067546" "GUT_GENOME068725"
#> [17] "GUT_GENOME080972" "GUT_GENOME090701"
#> [19] "GUT_GENOME094995" "GUT_GENOME096045"
#> [21] "GUT_GENOME096080" "GUT_GENOME096083"
#> [23] "GUT_GENOME096473" "GUT_GENOME096573"
#> [25] "GUT_GENOME102034" "GUT_GENOME103721"
#> [27] "GUT_GENOME104570" "GUT_GENOME109880"
#> [29] "GUT_GENOME112794" "GUT_GENOME113322"
#> [31] "GUT_GENOME114679" "GUT_GENOME115272"
#> [33] "GUT_GENOME115357" "GUT_GENOME116258"
#> [35] "GUT_GENOME116897" "GUT_GENOME117271"
#> [37] "GUT_GENOME132077" "GUT_GENOME135463"
#> [39] "GUT_GENOME140076" "GUT_GENOME142015"
#> [41] "GUT_GENOME142390" "GUT_GENOME143131"
#> [43] "GUT_GENOME143211" "GUT_GENOME143348"
#> [45] "GUT_GENOME143497" "GUT_GENOME143505"
#> [47] "GUT_GENOME149497" "GUT_GENOME156849"
#> [49] "GUT_GENOME174809" "GUT_GENOME175554"
#> [51] "GUT_GENOME189814" "GUT_GENOME195293"
#> [53] "GUT_GENOME208589" "GUT_GENOME210309"
#> [55] "GUT_GENOME210710" "GUT_GENOME217823"
#> [57] "GUT_GENOME217842" "GUT_GENOME217850"
#> [59] "GUT_GENOME234840" "GUT_GENOME252930"
#> [61] "GUT_GENOME258721" "GUT_GENOME261411"
#> [63] "GUT_GENOME272874" "GUT_GENOME274362"
#> [65] "GUT_GENOME275708" "GUT_GENOME277090"
#> [67] "GUT_GENOME284693" "GUT_GENOME286118"
```

This recalculates the minor allele frequency based on pooled SNV information (not on individual SNV information), thus takes a long time in case the number of profiled positions is large. We can set `cl` argument if needed.


```r
instr_chk <- "GUT_GENOME142015"
instr <- loadInStrain("../inStrain_out", instr_chk, skip_pool=FALSE) ## Load MAF table
#> Loading allele count table
#> Loading key table
#> Loading info table
#> Candidate species: GUT_GENOME142015
#>   Candidate key numbers: 1
#>   Dimension of pooled SNV table for species: 2359827
#> Calculating MAF
instr
#> Type: InStrain
#> Directory: ../inStrain_out
#> Species number: 68
#> Loaded SNV table: 1
#> Loaded gene table (): 0
#> 70.4 Mb
```

## metaSNV

## GT-Pro
