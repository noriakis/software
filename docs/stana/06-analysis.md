# Analysis example of gut microbiome of ESRD patients

We examine the dataset investigating the gut microbiome of ESRD patients ([Zhang et al. 2023](https://doi.org/10.1186/s13059-023-03056-y)) using `stana`. The SNV was profiled by MIDAS2 pipeline, and loaded to stana object. This time, we investigate one of the profiled species, `Faecalicatena gnavus` genotyping results.

Load the MIDAS2 merge results to stana object with setting the clinical grouping.


```r
library(dplyr)
library(stana)

meta <- read.table("../clinical/HDSubset/metadata.tsv", sep="\t", header=1)
hd <- meta$Run %>% setNames(meta$Group)
cl <- split(hd, names(hd))
cl <- lapply(cl, unname)
cand_species <- "101380"
stana <- loadMIDAS2("../clinical/HDSubset", candSp=cand_species, cl=cl, db="uhgg")
#>   101380
#>   g__Faecalicatena;s__Faecalicatena gnavus
#>     Number of snps: 26540
#>     Number of samples: 102
#>   101380
#>   g__Faecalicatena;s__Faecalicatena gnavus
#>     Number of genes: 59925
#>     Number of samples: 124
stana <- changeColors(stana, c("steelblue","gold","tomato"))
stana
#> # A stana: MIDAS2
#> # Database: uhgg
#> # Loaded directory: ../clinical/HDSubset
#> # Species number: 1
#> # Group info (list): CKD/HC/HD
#> # Loaded SNV table: 1 ID: 101380
#> # Loaded gene table: 1 ID: 101380
#> # Size: 108701376 B
```
Get a brief overview of SNVs.


```r
plotSNVInfo(stana, cand_species)
```

<img src="06-analysis_files/figure-html/app2-1.png" width="672" />

```r
plotSNVSummary(stana, cand_species) +
    scale_y_log10()
```

<img src="06-analysis_files/figure-html/app2-2.png" width="672" />

Based on the SNV and the related statistics of SNV, the consensus multiple sequence alignment is made by `consensusSeq` function. It can be accessed by `getFasta` function.


```r
stana <- consensusSeq(stana, cand_species)
#> # Beginning calling for 101380
#> # Original Site number: 26540
#>   Profiled samples: 102
#>   Included samples: 102
getFasta(stana)[[cand_species]]
#> 102 sequences with 24641 character and 24559 different site patterns.
#> The states are a c g t
```

Based on the MSA, the phylogenetic tree can be inferred by `inferAndPlotTree`. Inferring function can be specified by `treeFun` and is default to `upgma` function in phangorn, but can be stated the other functions. The tree plot shows the cladogram by default, but can be overridden by `branch.length` argument passed to `ggtree`.


```r
stana <- inferAndPlotTree(stana, cand_species, target="fasta", treeFun="FastTree")
#> File already exists!
#> Tree file already exists!
getTree(stana)[[cand_species]]
#> 
#> Phylogenetic tree with 102 tips and 100 internal nodes.
#> 
#> Tip labels:
#>   ERR11865846, ERR11865866, ERR11865898, ERR11865921, ERR11865925, ERR11865952, ...
#> 
#> Unrooted; includes branch lengths.
getTreePlot(stana)[[cand_species]]
```

<img src="06-analysis_files/figure-html/app4-1.png" width="672" />

Using cophenetic distance matrix from tree, the PERMANOVA is performed and the principal coordinate analysis plot based on the distance matrix is plotted.


```r
stana <- doAdonis(stana, cand_species, target="tree", pcoa=TRUE)
#> # Performing adonis in 101380 target is tree
#> Warning in att$heading[2] <- deparse(match.call(),
#> width.cutoff = 500L): number of items to replace is not a
#> multiple of replacement length
#> #  F: 2.70213130791564, R2: 0.0517628541251898, Pr: 0.031
```

<img src="06-analysis_files/figure-html/app5-1.png" width="672" />

```r
getAdonis(stana)[[cand_species]]
#> Permutation test for adonis under reduced model
#> Terms added sequentially (first to last)
#> Permutation: free
#> Number of permutations: 999
#> 
#> adonis2(formula = d ~ ., data = structure(list(group = c("CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", "HD", 
#>           Df SumOfSqs      R2      F Pr(>F)  
#> group      2   0.2860 0.05176 2.7021  0.031 *
#> Residual  99   5.2391 0.94824                
#> Total    101   5.5251 1.00000                
#> ---
#> Signif. codes:  
#> 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Based on the results, assuming there are multiple factors within the species, we can directly examine their estimated abundances and the functional implications by using the NMF approach in KO table. First we calculate KO abundances based on eggNOG-mapper annotation.


```r
stana <- setAnnotation(stana,
                       annotList=list("101380"="../annotations_uhgg/101380_eggnog_out.emapper.annotations"))
stana <- calcGF(stana, candSp=cand_species)
```

Subsequently, using the cross-validation approach replacing random variables to `NA`, we estimate the rank based on the KO tables. The function in `NNLM` is used to compute the loss (`mse`) with replacing 30% of the data to `NA`.


```r
library(NNLM)
cvl <- list()
for (i in seq_len(5)) {
  cvl[[i]] <- NMF(stana, cand_species, nnlm_flag=TRUE,
                  nnlm_args=list("loss"="mse"),
                  estimate=TRUE, tss=TRUE)[1,]
}
#> # NMF started 101380, target: kos, method: NNLM::nnmf
#> # Original features: 2568
#> # Original samples: 124
#> # Original matrix NA: NA
#> # Original matrix zero: 0.487
#> # Performing TSS
#> # Filtered features: 2568
#> # Filtered samples: 124
#> # NNLM flag enabled, the error matrix only will be returned.
#> # NMF started 101380, target: kos, method: NNLM::nnmf
#> # Original features: 2568
#> # Original samples: 124
#> # Original matrix NA: NA
#> # Original matrix zero: 0.487
#> # Performing TSS
#> # Filtered features: 2568
#> # Filtered samples: 124
#> # NNLM flag enabled, the error matrix only will be returned.
#> # NMF started 101380, target: kos, method: NNLM::nnmf
#> # Original features: 2568
#> # Original samples: 124
#> # Original matrix NA: NA
#> # Original matrix zero: 0.487
#> # Performing TSS
#> # Filtered features: 2568
#> # Filtered samples: 124
#> # NNLM flag enabled, the error matrix only will be returned.
#> # NMF started 101380, target: kos, method: NNLM::nnmf
#> # Original features: 2568
#> # Original samples: 124
#> # Original matrix NA: NA
#> # Original matrix zero: 0.487
#> # Performing TSS
#> # Filtered features: 2568
#> # Filtered samples: 124
#> # NNLM flag enabled, the error matrix only will be returned.
#> # NMF started 101380, target: kos, method: NNLM::nnmf
#> # Original features: 2568
#> # Original samples: 124
#> # Original matrix NA: NA
#> # Original matrix zero: 0.487
#> # Performing TSS
#> # Filtered features: 2568
#> # Filtered samples: 124
#> # NNLM flag enabled, the error matrix only will be returned.

do.call(rbind, cvl) %>% data.frame(check.names=FALSE) %>% mutate(group=1:5) %>%
  tidyr::pivot_longer(1:6) %>% ggplot(aes(x=name, y=value))+
  geom_boxplot() + cowplot::theme_cowplot() + xlab("rank")
```

<img src="06-analysis_files/figure-html/hd7-1.png" width="672" />

Based on the information, the factor number of two is selected.


```r
set.seed(1)
stana <- NMF(stana, cand_species, rank=2,
	nnlm_flag=TRUE, nnlm_args=list("loss"="mse"))
#> # NMF started 101380, target: kos, method: NNLM::nnmf
#> # Original features: 2568
#> # Original samples: 124
#> # Original matrix NA: NA
#> # Original matrix zero: 0.487
#> # Filtered features: 2568
#> # Filtered samples: 124
#> # Rank 2
#> Mean relative abundances: 0.5826988 0.4173012 
#> Present feature per factor: 2231 2394

## Plot the results
plotAbundanceWithinSpecies(stana, cand_species, by="coef")
```

<img src="06-analysis_files/figure-html/hd9-1.png" width="672" />

```r
plotStackedBarPlot(stana, cand_species, by="coef") + scale_fill_manual(values=c("tomato","gold"))
#> Using sample, group as id variables
```

<img src="06-analysis_files/figure-html/hd9-2.png" width="672" />

Using these two factors, we summarize KO abundance information to KEGG PATHWAY information, and plot the relationship between the pathway abundance within two factors by scatter plot and heatmap.


```r
library(ggrepel)
#> Warning: package 'ggrepel' was built under R version 4.3.2
pw <- data.frame(pathwayWithFactor(stana, cand_species, tss=TRUE, change_name=TRUE,
	mat = getSlot(stana, "NMF")[[cand_species]]$W))
colnames(pw) <- c("1","2")
pw[["name"]] <- row.names(pw)
pw[["size"]] <- (pw[,1] + pw[,2])/2
ggplot(pw, aes(x=pw[,1], y=pw[,2]))+
    geom_point(aes(size=size))+
    geom_text_repel(aes(label=name, size=size), bg.colour="white")+
    geom_smooth(method=lm, se=TRUE)+xlab("1")+ylab("2")+
    cowplot::theme_cowplot()
#> `geom_smooth()` using formula = 'y ~ x'
#> Warning: ggrepel: 225 unlabeled data points (too many
#> overlaps). Consider increasing max.overlaps
```

<img src="06-analysis_files/figure-html/app9-1.png" width="768" />

```r

## Sort by absolute difference
fc <- pw[,1] - pw[,2]
names(fc) <- pw[["name"]]
nms <- names(sort(abs(fc[!is.infinite(fc)]), decreasing=TRUE) %>% head(40))

library(pheatmap)
pheatmap(pw[nms, 1:2])
```

<img src="06-analysis_files/figure-html/app9-2.png" width="768" />
Of these, cysteine and methionine metabolism pathway is interesting as the pathway is reported to be related to the species. The KEGG PATHWAY scheme of the pathway is plotted by ggkegg (For group comparison, use `plotKEGGPathway`. The returned object is ggplot object and the users can modify the visualization by stacking the layers).

The colors in the nodes of left-side is abundance for factor 1 and right side is factor 2.


```r
## Built-in `plotKEGGPathway` function. The statistics to be shown is moderated t-value
# kegg <- plotKEGGPathway(stana, cand_species,
#                         pathway_id="ko00270",
#                         statMethod="mod.t")

library(ggkegg)
library(tidygraph)

pp <- ggkegg::pathway("ko00270") %N>% 
    mutate(
    	f1=ggkegg::node_numeric(getSlot(stana, "NMF")[[cand_species]]$W[,1]),
    	f2=ggkegg::node_numeric(getSlot(stana, "NMF")[[cand_species]]$W[,2])
    )

gg <- ggraph(pp, layout="manual", x=x, y=y)+
    geom_node_rect(aes(fill=f1, xmin=xmin, xmax=x, filter=type=="ortholog"))+
    geom_node_rect(aes(fill=f2, xmin=x, xmax=xmax, filter=type=="ortholog"))+
    scale_fill_gradient(low="blue",high="pink", name="abundance")+
    overlay_raw_map() +
    stamp("ko:K00789")+
    theme_void()
gg
```

<img src="06-analysis_files/figure-html/ggkegg-1.png" width="672" />

In this map, we can find interesting findings like one of the enzymes AdoMet synthetase (2.5.1.6, indicated in red rectangle), is enriched in the factor 2, and the factor 2 is elevated in HD. The corresponding enzyme is reported to be in relation to the hemodialysis ([Loehrer et al. 1998.](https://doi.org/10.1093/ndt/13.3.656)). The results suggest the library and function can link the intra-species diversity and clinical factors in the R environment.


Finally, the results can be exported to the interactive inspection by `exportInteractive` function for the sharing with the other researchers.


```r
exportInteractive(stana)
#> Warning in dir.create(paste0(out, "/data")): '.\data'
#> already exists
#> Tree number: 1, KO (or gene) number: 1
#> Exporting ...
#> # A stana: MIDAS2
#> # Database: uhgg
#> # Loaded directory: ../clinical/HDSubset
#> # Species number: 1
#> # Group info (list): CKD/HC/HD
#> # Loaded SNV table: 1 ID: 101380
#> # Loaded gene table: 1 ID: 101380
#> # Loaded KO table: 1 ID: 101380
#> # Inferred fasta: 1 ID: 101380
#> # Size: 125115648 B
```
