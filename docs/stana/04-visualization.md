# Visualization

`stana` offers visualization functions that can interpret the profiled data.


```r
library(stana)
library(ComplexHeatmap)
```

## `plotSNVSummary` and `plotSNVInfo`

These functions can be used to plot the information related to SNV, which can be used to inspect the overview of profiled SNVs.


```r
load("../hd_meta.rda")
stana <- loadMIDAS2("../merge_uhgg", cl=hd_meta, candSp="102478")
#>   102478
#>     Number of snps: 77431
#>     Number of samples: 31
#>   102478
#>     Number of genes: 150996
#>     Number of samples: 32
plotSNVSummary(stana, "102478")
```

<img src="04-visualization_files/figure-html/p-1.png" width="672" />

```r
plotSNVInfo(stana, "102478")
```

<img src="04-visualization_files/figure-html/p-2.png" width="672" />

## `plotMAF`

This functions plots the boxplot of MAF between groups.


```r
plotMAF(stana, "102478", row.names(getSlot(stana, "snpsInfo")[[1]])[1])
```

<img src="04-visualization_files/figure-html/pm-1.png" width="672" />

## `plotCoverage`

This function plots the coverage of the canddiate species across the group.


```r
plotCoverage(stana, "102478")
```

<img src="04-visualization_files/figure-html/pc-1.png" width="672" />


## `plotPCA`

You can plot principal component analysis results using `prcomp` by `plotPCA` function across grouping specified. If no group is specified, `No_Group` label is assigned to all the samples.


```r
mt <- loadmetaSNV("../metasnv_sample_out/",just_species = TRUE)
mt <- loadmetaSNV("../metasnv_sample_out/", candSp=mt[1])
#>   Loading refGenome1clus
plotPCA(mt, species=getID(mt)[1])
#> After filtering: 1022 SNVs
#> $refGenome1clus
```

<img src="04-visualization_files/figure-html/plotpca-1.png" width="672" />


## `plotCirclize`

This function can be used to circlize plot, which can link the information related to MAF or SNV with gene copy number analysis. The function needs the `snpsInfo` slot filled with the data having `genome_id`, `gene_id` and `position` column. For MIDAS2 and inStrain, the packge automatically infers the `genome_id`.


```r
library(circlize)
#> ========================================
#> circlize version 0.4.15
#> CRAN page: https://cran.r-project.org/package=circlize
#> Github page: https://github.com/jokergoo/circlize
#> Documentation: https://jokergoo.github.io/circlize_book/book/
#> 
#> If you use it in published research, please cite:
#> Gu, Z. circlize implements and enhances circular visualization
#>   in R. Bioinformatics 2014.
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(circlize))
#> ========================================
#> 
#> Attaching package: 'circlize'
#> The following object is masked from 'package:igraph':
#> 
#>     degree
plotCirclize(stana, "102478", genomeId="UHGG143505", thresh_snp_gene = 150)
#> Type is MIDAS2
#> Features not provided, default to sample_counts
#> Genome ID in SNV information:
#>   UHGG143505_1: 179 - 5439351, number of position: 77431
#> Genome ID: UHGG143505
#> Included position: 987
#> Note: 1 point is out of plotting region in sector
#> 'UHGG143505_00169', track '1'.
#> Note: 16 points are out of plotting region in sector
#> 'UHGG143505_00169', track '1'.
#> Note: 1 point is out of plotting region in sector
#> 'UHGG143505_01855', track '1'.
#> Note: 16 points are out of plotting region in sector
#> 'UHGG143505_01855', track '1'.
#> Note: 1 point is out of plotting region in sector
#> 'UHGG143505_01894', track '1'.
#> Note: 16 points are out of plotting region in sector
#> 'UHGG143505_01894', track '1'.
#> Note: 1 point is out of plotting region in sector
#> 'UHGG143505_02021', track '1'.
#> Note: 16 points are out of plotting region in sector
#> 'UHGG143505_02021', track '1'.
#> Note: 1 point is out of plotting region in sector
#> 'UHGG143505_03685', track '1'.
#> Note: 16 points are out of plotting region in sector
#> 'UHGG143505_03685', track '1'.
```

<img src="04-visualization_files/figure-html/plotcirc-1.png" width="672" />

## Visualization of gene copy numbers

The violin plot of gene copy numbers can be plotted by `plotGenes`.


```r
load("../hd_meta.rda")
stana <- loadMIDAS2("../merge_uhgg", cl=hd_meta, candSp="102478")
#>   102478
#>     Number of snps: 77431
#>     Number of samples: 31
#>   102478
#>     Number of genes: 150996
#>     Number of samples: 32
plotGenes(stana, "102478", c("UHGG000186_00531","UHGG000186_00521"))
```

<img src="04-visualization_files/figure-html/gab-1.png" width="672" />

Default color mapping can be changed by `changeColors`.


```r
stana <- changeColors(stana, c("blue","red"))
plotGenes(stana, "102478", c("UHGG000186_00531","UHGG000186_00521"))
```

<img src="04-visualization_files/figure-html/gab2-1.png" width="672" />

## Visualization of phylogenetic tree

See \@ref(cons) for the functions like `inferAndPlotTree` for visualizing the phylogenetic tree inferred by various methods.

## Visualization of functional analysis results

See \@ref(function) for the functions like `plotHeatmap` and `plotKEGGPathway`.

## Visualization of `inStrain` results

For the specific software, the imported `inStrain` `compare` profiles can be visualized. The loaded genome-wide comparison table and strain cluster table can be visualized using `genomeHeatmap` and `strainClusterHeatmap` by `ComplexHeatmap`. For `genomeHeatmap`, typically population ANI or consensus ANI are plotted, but all the columns listed in `genomeWide_compare.tsv` can be plotted. The parameters to be passed to `Heatmap` can be specified with `heatmapArgs`. If cluster information (`getCl(stana)`) is available or `cl` is specified, the columns will be split to present the grouping. Please refer to the documentation of [`inStrain`](https://instrain.readthedocs.io/en/latest/important_concepts.html) for `popANI` and `conANI`.


```r
instr_chk <- "GUT_GENOME142015"
instr <- loadInStrain("../inStrain_out", instr_chk)
genomeHeatmap(instr, instr_chk, column = "popANI", heatmapArgs = list(show_column_name=FALSE))
```

<img src="04-visualization_files/figure-html/heatmaps_instrain-1.png" width="672" />

```r
strainClusterHeatmap(instr, instr_chk, heatmapArgs = list(show_column_name=FALSE))
```

<img src="04-visualization_files/figure-html/heatmaps_instrain-2.png" width="672" />
