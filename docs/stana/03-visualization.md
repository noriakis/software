# Visualization

`stana` offers visualization functions that can interpret the profiled data.


```r
library(stana)
library(ComplexHeatmap)
```


## Visualization of gene abundances


```r
load("../hd_meta.rda")
stana <- loadMIDAS2("../merge_uhgg", cl=hd_meta, candSp="102478")
#> SNPS
#>   102478
#>     Number of snps: 77431
#>     Number of samples: 31
#>       Number of samples in HC: 9
#>       Number of samples in R: 22
#>       Passed the filter
#> Genes
#>   102478
#>     Number of genes: 150996
#>     Number of samples: 32
#>       Number of samples in HC: 10
#>       Number of samples in R: 22
#>       Passed the filter
plotGenes(stana, "102478", c("UHGG000186_00531","UHGG000186_00521"))
```

<img src="03-visualization_files/figure-html/gab-1.png" width="672" />

## Visualization of functional analysis results

See \@ref(function).


## Visualization of `inStrain` results

The imported `inStrain` `compare` profiles can be visualized. The loaded genome-wide comparison table and strain cluster table can be visualized using `genomeHeatmap` and `strainClusterHeatmap` by `ComplexHeatmap`. For `genomeHeatmap`, typically population ANI or consensus ANI are plotted, but all the columns listed in `genomeWide_compare.tsv` can be plotted. The parameters to be passed to `Heatmap` can be specified with `heatmapArgs`. If cluster information (`stana@cl`) is available or `cl` is specified, the columns will be split to present the grouping.


```r
instr_chk <- "GUT_GENOME142015"
instr <- loadInStrain("../inStrain_out", instr_chk)
genomeHeatmap(instr, instr_chk, column = "popANI", heatmapArgs = list(show_column_name=FALSE))
```

<img src="03-visualization_files/figure-html/heatmaps_instrain-1.png" width="672" />

```r
strainClusterHeatmap(instr, instr_chk, heatmapArgs = list(show_column_name=FALSE))
```

<img src="03-visualization_files/figure-html/heatmaps_instrain-2.png" width="672" />
