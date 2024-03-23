# Importing and filtering

`stana` is aimed to import the metagenotyping results of various pipelines.
The package is designed to handle some types of the metagenotyping software such as `MIDAS` and `MIDAS2`, which outputs the 
gene copy numbers by default. Note that the slot name `snps` here refers to just the variable, and not to reflect the actual meaning. The species ID in `stana` object is interpreted as **character**.


```r
library(stana)
```

Along with specifying which directory to load, the grouping variables should be set to stana object. By default, named list of groups is used (categorical). The `cl` argument in the functions and `cl` slot corresponds to this information. Additionally, the metadata in ordinaly data.frame can be set by `setMetadata`. The row names or `id` column should represent the sample names.

## MIDAS

For `MIDAS`, `loadMIDAS` function can be used to import the output of `merge` command. In `MIDAS` and `MIDAS2` examples, we load the example dataset deposited by the study investigating gut microbiome of hemodialysis patients ([Shi et al. 2022](https://doi.org/10.3389/fcimb.2022.904284)). `hd_meta` includes named list of grouping. First we would like to see how many samples were profiled in the species.


```r
load("../hd_meta.rda")
stana <- loadMIDAS("../merge_midas1", cl=hd_meta, only_stat=TRUE)
stana$snps |> head() |> DT::datatable()
```


```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-8635d10534f1a103e9e9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-8635d10534f1a103e9e9">{"x":{"filter":"none","vertical":false,"data":[["Acidaminococcus_intestini_54097","Akkermansia_muciniphila_55290","Alistipes_finegoldii_56071","Alistipes_indistinctus_62207","Alistipes_onderdonkii_55464","Alistipes_putredinis_61533"],["Acidaminococcus_intestini_54097","Akkermansia_muciniphila_55290","Alistipes_finegoldii_56071","Alistipes_indistinctus_62207","Alistipes_onderdonkii_55464","Alistipes_putredinis_61533"],["1","3","3","0","7","7"],["5","8","5","1","14","9"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>species<\/th>\n      <th>HC<\/th>\n      <th>R<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

We will load the interesting species.


```r
stana <- loadMIDAS("../merge_midas1/", cl=hd_meta, candSp="Bacteroides_uniformis_57318")
#> Bacteroides_uniformis_57318
#>   Snps
#>     HC 13
#>     R 16
#>     Bacteroides_uniformis_57318 cleared filtering threshold in SNV
#>   Genes
#>     HC 13
#>     R 16
#>     Bacteroides_uniformis_57318 cleared filtering threshold in genes
#> Overall, 1 species met criteria in SNPs
#> Overall, 1 species met criteria in genes
stana
#> # A stana: MIDAS1
#> # Loaded directory: ../merge_midas1/
#> # Species number: 1
#> # Group info (list): HC/R
#> # Loaded SNV table: 1 ID: Bacteroides_uniformis_57318
#> # Loaded gene table: 1 ID: Bacteroides_uniformis_57318
#> # Size: 14586224 B
```

## MIDAS2

For `MIDAS2`, `loadMIDAS2` function can be used to import the output of `merge` command.

::: rmdwarning
For `MIDAS2`, the function assumes the default database of UHGG or GTDB is used.
:::


```r
load("../hd_meta.rda")
hd_meta
#> $HC
#>  [1] "ERR9492498" "ERR9492502" "ERR9492507" "ERR9492508"
#>  [5] "ERR9492512" "ERR9492497" "ERR9492499" "ERR9492501"
#>  [9] "ERR9492504" "ERR9492505" "ERR9492509" "ERR9492511"
#> [13] "ERR9492500" "ERR9492503" "ERR9492506" "ERR9492510"
#> 
#> $R
#>  [1] "ERR9492489" "ERR9492490" "ERR9492492" "ERR9492494"
#>  [5] "ERR9492519" "ERR9492524" "ERR9492526" "ERR9492527"
#>  [9] "ERR9492528" "ERR9492513" "ERR9492514" "ERR9492518"
#> [13] "ERR9492520" "ERR9492522" "ERR9492523" "ERR9492525"
#> [17] "ERR9492491" "ERR9492493" "ERR9492495" "ERR9492496"
#> [21] "ERR9492515" "ERR9492516" "ERR9492517" "ERR9492521"
```

We can check stats of how many samples are profiled for each species, by `only_stat`. This returns the list of tibbles with names `snps` and `genes`.


```r
stana <- loadMIDAS2("../merge_uhgg", only_stat=TRUE, cl=hd_meta)
stana$snps |> dplyr::filter(group=="HC") |> dplyr::arrange(desc(n)) |> head()
#> # A tibble: 6 × 4
#> # Groups:   species_id [6]
#>   species_id group     n species_description           
#>   <chr>      <chr> <int> <chr>                         
#> 1 101346     HC       12 s__KS41 sp003584895           
#> 2 102438     HC       10 s__Synechococcus_C sp001632165
#> 3 101378     HC        9 s__CAG-873 sp002490635        
#> 4 102478     HC        9 s__Clostridium_A leptum       
#> 5 102492     HC        8 s__Thalassospira lucentensis_A
#> 6 100044     HC        7 s__Helicobacter pylori_BU
```
As the long output is expected, only one species is loaded here. 


```r
stana <- loadMIDAS2("../merge_uhgg", candSp="100002", cl=hd_meta)
#>   100002
#>     Number of snps: 2058
#>     Number of samples: 5
#>   100002
#>     Number of genes: 23427
#>     Number of samples: 7
```

The data is profiled against UHGG. `loadSummary` and `loadInfo` can be specified to load the SNV summary and SNV info per species, which is default to `TRUE`.


```r
getSlot(stana, "snps")[["100002"]] |> head()
#>                                 ERR9492497 ERR9492515
#> gnl|Prokka|UHGG000004_1|2901|A           1          1
#> gnl|Prokka|UHGG000004_1|4071|C           1          0
#> gnl|Prokka|UHGG000004_1|11094|T          1          0
#> gnl|Prokka|UHGG000004_1|11148|T          1          0
#> gnl|Prokka|UHGG000004_1|11940|G          0          0
#> gnl|Prokka|UHGG000004_1|11970|C          0          0
#>                                 ERR9492526 ERR9492527
#> gnl|Prokka|UHGG000004_1|2901|A           0      0.429
#> gnl|Prokka|UHGG000004_1|4071|C           1      0.000
#> gnl|Prokka|UHGG000004_1|11094|T          0      0.000
#> gnl|Prokka|UHGG000004_1|11148|T          0      0.429
#> gnl|Prokka|UHGG000004_1|11940|G          1      0.444
#> gnl|Prokka|UHGG000004_1|11970|C          1      0.556
#>                                 ERR9492528
#> gnl|Prokka|UHGG000004_1|2901|A           0
#> gnl|Prokka|UHGG000004_1|4071|C           0
#> gnl|Prokka|UHGG000004_1|11094|T          1
#> gnl|Prokka|UHGG000004_1|11148|T          1
#> gnl|Prokka|UHGG000004_1|11940|G          0
#> gnl|Prokka|UHGG000004_1|11970|C          0
getSlot(stana, "freqTableSnps") |> head()
#> data frame with 0 columns and 0 rows
```

For ID conversion, the metadata accompanied with the default database can be used.


```r
taxtbl <- read.table("../metadata_uhgg.tsv", sep="\t",
                     header=1, row.names=1, check.names = FALSE)
taxtbl |> head()
#>          representative MGnify_accession species_closest
#> 100001 GUT_GENOME000001  MGYG-HGUT-00001          100049
#> 100002 GUT_GENOME000004  MGYG-HGUT-00002          100201
#> 100003 GUT_GENOME000008  MGYG-HGUT-00003          103279
#> 100004 GUT_GENOME000010  MGYG-HGUT-00004          103876
#> 100005 GUT_GENOME000017  MGYG-HGUT-00005          101623
#> 100006 GUT_GENOME000020  MGYG-HGUT-00006          100011
#>        ani_closest gtpro_kmer_counts phyeco_marker_counts
#> 100001    83.35600                NA                   15
#> 100002    84.86955             29822                   15
#> 100003    93.41555              1115                   15
#> 100004    78.17620                NA                   15
#> 100005    85.84275                NA                   15
#> 100006    93.07075                NA                   15
#>        phyeco_pass_ratio pangene_counts genome_counts
#> 100001                 1           4893             4
#> 100002                 1         147601           358
#> 100003                 1         113409          1178
#> 100004                 1          12599            24
#> 100005                 1           5488             2
#> 100006                 1           2656             1
#>        Genome_type  Length N_contigs    N50 GC_content
#> 100001     Isolate 3219614       137  47258      28.26
#> 100002     Isolate 4433090       100 109266      42.60
#> 100003     Isolate 3229507        35 158570      58.52
#> 100004     Isolate 3698872       105  90296      54.19
#> 100005     Isolate 3930422        32 350032      28.59
#> 100006     Isolate 2822523        36 121380      32.65
#>        Completeness Contamination
#> 100001        98.59          0.70
#> 100002        99.37          0.00
#> 100003       100.00          0.00
#> 100004        98.66          0.22
#> 100005        99.30          0.00
#> 100006        99.26          1.39
#>                                                                                                                                                Lineage
#> 100001                                 d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Peptostreptococcales;f__Peptostreptococcaceae;g__GCA-900066495;s__
#> 100002                            d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia_A;s__Blautia_A sp900066165
#> 100003                                   d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;g__Alistipes;s__Alistipes shahii
#> 100004                   d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Anaerotruncus;s__Anaerotruncus colihominis
#> 100005 d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Peptostreptococcales;f__Peptostreptococcaceae;g__Terrisporobacter;s__Terrisporobacter glycolicus_A
#> 100006                       d__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus xylosus
#>        Continent
#> 100001    Europe
#> 100002    Europe
#> 100003    Europe
#> 100004    Europe
#> 100005    Europe
#> 100006    Europe
```

The taxonomy table can be loaded with providing to `taxtbl` argument.


```r
loadMIDAS2("../merge_uhgg", cl=hd_meta, candSp="100002", taxtbl=taxtbl, db="uhgg")
#>   100002
#>   d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia_A;s__Blautia_A sp900066165
#>     Number of snps: 2058
#>     Number of samples: 5
#>   100002
#>   d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia_A;s__Blautia_A sp900066165
#>     Number of genes: 23427
#>     Number of samples: 7
#> # A stana: MIDAS2
#> # Database: uhgg
#> # Loaded directory: ../merge_uhgg
#> # Species number: 1
#> # Group info (list): HC/R
#> # Loaded SNV table: 1 ID: 100002
#> # Loaded gene table: 1 ID: 100002
#> # Size: 4525192 B
#> # 
#> # SNV description
#> # A tibble: 2 × 3
#> # Groups:   group [2]
#>   group species_id                                         n
#>   <chr> <chr>                                          <int>
#> 1 HC    d__Bacteria;p__Firmicutes_A;c__Clostridia;o__…     1
#> 2 R     d__Bacteria;p__Firmicutes_A;c__Clostridia;o__…     4
```

The coverage for each species per sample is plotted by `plotCoverage`.


```r
plotCoverage(stana, "100002", pointSize=5)
```

<img src="01-importing_files/figure-html/plotcoverage-1.png" width="672" />

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
#> Loading the large table...
#> Loading key table
#> Loading info table
#> Candidate species: GUT_GENOME142015
#>   Candidate key numbers: 1
#>   Dimension of pooled SNV table for species: 2359827
#> Calculating MAF
instr
#> # A stana: InStrain
#> # Loaded directory: ../inStrain_out
#> # Species number: 1
#> # Loaded SNV table: 1 ID: GUT_GENOME142015
#> # Size: 81502904 B
```

## metaSNV

For loading the output of `metaSNV`, `metaSNV.py` and `metaSNV_Filtering.py` are typically performed beforehand. You can use `just_species` to return species ID. Note that the `loadmetaSNV` is currently supported to load only the SNV profiles.


```r
meta <- loadmetaSNV("../metasnv_sample_out")
#>   Loading refGenome1clus
#>   Loading refGenome2clus
#>   Loading refGenome3clus
```

## Manual

The loading of the manually created data.frame is possible by `snv` and `gene` function. The metagenotyping results produced by the other software can be loaded with this function, although some statistics should be inserted manually. This accepts the named list of data.frame.


```r
man <- snv(list("manual"=head(getSlot(meta, "snps")[[1]])))
man
#> # A stana: manual
#> # Loaded directory: 
#> # Species number: 1
#> # Loaded SNV table: 1 ID: manual
#> # Size: 49432 B
```

### Conversion to MAF matrix

Some metagenotype data comes with TSV with the number of allele count with long or wide format.
The package offers some functions to convert these to MAF matrix.


```r
df1 <- data.frame(rbind(
  c("SPID1","SP1_NZ_GG770218.1_131457", "G", "A", "7", "2"),
  c("SPID1","SP1_NZ_GG770218.1_131458", "C", "T", "2", "5"),
  c("SPID2","SP2_NZ.1_131459", "A", "T", "7", "13"),
  c("SPID2","SP2_NZ.1_131470", "A", "T", "12", "2")
))

head(df1)
#>      X1                       X2 X3 X4 X5 X6
#> 1 SPID1 SP1_NZ_GG770218.1_131457  G  A  7  2
#> 2 SPID1 SP1_NZ_GG770218.1_131458  C  T  2  5
#> 3 SPID2          SP2_NZ.1_131459  A  T  7 13
#> 4 SPID2          SP2_NZ.1_131470  A  T 12  2

convertWideToMaf(df1)[[1]]
#> $SPID1
#>                                maf
#> SP1_NZ_GG770218.1_131457 0.2222222
#> SP1_NZ_GG770218.1_131458 0.2857143
#> 
#> $SPID2
#>                       maf
#> SP2_NZ.1_131459 0.3500000
#> SP2_NZ.1_131470 0.1428571



df2 <- data.frame(rbind(
  c("SPID3","SP3_AZ.1_131457", "G", "A", "7", "2"),
  c("SPID3","SP3_AZ.1_131458", "C", "T", "2", "5"),
  c("SPID1","SP1_NZ_GG770218.1_131457", "A", "T", "7", "13"),
  c("SPID1","SP1_NZ_GG770218.1_131458", "A", "T", "12", "2")
))
```

They can be merged to make stana object by `combineMaf()` function.



```r
mafs <- list("sample1"=convertWideToMaf(df1), "sample2"=convertWideToMaf(df2))
stana <- combineMaf(mafs)
stana
#> # A stana: manual
#> # Loaded directory: 
#> # Species number: 0
#> # Loaded SNV table: 3 ID: SPID1
#> # Size: 23368 B
getSlot(stana, "snpsSummary")
#> data frame with 0 columns and 0 rows

## If needed, SNV-level statistics should be set to stana
# stana <- setSlot(stana, "snpsInfo", list("SPID1"=spid1_statistics))
```

This object can be used with the subsequent downstream analysis. The function can be useful in importing the metagenotyping software such as `GT-Pro`.

## Filtering by species ID

The `check` method applied to `stana` object can check the samples or species met the specific criteria.


```r
stana <- loadMIDAS2("../merge_uhgg", cl=hd_meta, candSp="100002", taxtbl=taxtbl, db="uhgg")
#>   100002
#>   d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia_A;s__Blautia_A sp900066165
#>     Number of snps: 2058
#>     Number of samples: 5
#>   100002
#>   d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia_A;s__Blautia_A sp900066165
#>     Number of genes: 23427
#>     Number of samples: 7
filt <- stana %>% check(mean_coverage > 4)
head(filt)
#> # A tibble: 6 × 3
#>   species_id     n species_description                      
#>        <int> <int> <chr>                                    
#> 1     100002     5 d__Bacteria;p__Firmicutes_A;c__Clostridi…
#> 2     100003    15 d__Bacteria;p__Bacteroidota;c__Bacteroid…
#> 3     100022     7 d__Bacteria;p__Firmicutes_A;c__Clostridi…
#> 4     100036     3 d__Bacteria;p__Firmicutes_A;c__Clostridi…
#> 5     100039     2 d__Bacteria;p__Firmicutes_A;c__Clostridi…
#> 6     100041     3 d__Bacteria;p__Firmicutes_A;c__Clostridi…
```

Subsequently, `filter` method can be used to subset the stana object for interesting species.


```r
stana <- stana %>% filter(unique(filt$species_id))
```

Most of the functions in stana is designed to perform the analysis in all the species within stana object unless specified, so the filtering beforehand with grouping variables is important (such as number of samples profiled per group).
