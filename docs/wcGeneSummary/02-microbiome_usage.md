

# The basic usage for microbiome-related data

The general usage for mining the information related to microbiome-related signatures.


```r
library(wcGeneSummary)
library(ggplot2)
library(ggraph)
library(RColorBrewer)
```

We use BugSigDB, and its R port bugsigdbr to obtain the curated dataset of the relationship with bacterial taxonomy and human diseases ([Geistlinger et al. 2022](https://bugsigdb.org/Main_Page)). Users can query microbiome names, which will be searched for MetaPhlAn taxonomic annotation. If `target="title"`, the title of the corresponding articles will be summarized.


```r
basic <- wcBSDB(c("Veillonella dispar","Neisseria flava"),tag=FALSE,
    curate=TRUE,target="title",pre=TRUE,cl=snow::makeCluster(12),
    pal=RColorBrewer::brewer.pal(10, "Set2"),numWords=80,argList=list(min.freq=1))
#> Input microbes: 2
#>   Found 17 entries for Veillonella dispar
#>   Found 1 entries for Neisseria flava
#> Including 28 entries
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
basic@freqDf |> head(n=20)
#>                    word freq
#> gut                 Gut    8
#> oral               Oral    5
#> patients       patients    5
#> study             study    3
#> arthritis     arthritis    2
#> association association    2
#> bacterial     bacterial    2
#> covid19         COVID19    2
#> diabetes       diabetes    2
#> dysbiosis     Dysbiosis    2
#> infant           infant    2
#> infection     Infection    2
#> obese             obese    2
#> respiratory respiratory    2
#> sequencing   sequencing    2
#> surgery         surgery    2
#> tract             tract    2
#>  features      features    1
#> 16s                 16S    1
#> aerosol         Aerosol    1
basic@wc
```

<img src="02-microbiome_usage_files/figure-html/bsdb_basic-1.png" width="768" />

If `target="abstract"`, the corresponding abstract will be fetched and be summarized.


```r
basic2 <- wcBSDB(c("Veillonella dispar","Neisseria flava"),tag=TRUE,
    curate=TRUE,target="abstract",pre=TRUE,cl=snow::makeCluster(12),
    pal=RColorBrewer::brewer.pal(10, "Set2"),numWords=80)
#> Input microbes: 2
#>   Found 17 entries for Veillonella dispar
#>   Found 1 entries for Neisseria flava
#> Including 28 entries
#> Target is abstract
#>   Querying PubMed for 17 pmids
#>   Querying without API key
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
#> Multiscale bootstrap... Done.
basic2@freqDf |> head()
#>                word freq
#> patients   patients   37
#> gut             Gut   30
#> oral           Oral   27
#> bacterial Bacterial   25
#> species     species   24
#> microbial Microbial   23
basic2@wc
```

<img src="02-microbiome_usage_files/figure-html/bsdb_basic2-1.png" width="768" />

For successful visualization, pre-caculated TF-IDF and frequency data frame is available and one can use them to filter the highly occurring words, or the other prefiltering option used in `wcGeneSummary`.


```r
rmwords <- wcGeneSummary:::allFreqBSDB
filter <- rmwords[rmwords$freq>quantile(rmwords$freq, 0.95),]
filter |> head(n=20)
#>             freq        word
#> microbiota   275  microbiota
#> gut          242         gut
#> microbiome   175  microbiome
#> patients     146    patients
#> study         69       study
#> cancer        64      cancer
#> composition   62 composition
#> oral          55        oral
#> human         49       human
#> intestinal    48  intestinal
#> children      45    children
#> microbial     45   microbial
#> disease       43     disease
#> fecal         41       fecal
#> alterations   36 alterations
#> analysis      34    analysis
#> association   34 association
#> dysbiosis     34   dysbiosis
#> infection     34   infection
#> risk          33        risk
```
The network visualization is possible by enabling `plotType="network"`.
The same parameters that can be passed to `wcGeneSummary` can be used.


```r
net <- wcBSDB(c("Neisseria","Veillonella"),
    curate=TRUE,target="title",pre=TRUE,plotType="network",
    additionalRemove=filter$word, corThresh=0.2, edgeLink=FALSE,
    numWords=60)
#> Input microbes: 2
#>   Found 76 entries for Neisseria
#>   Found 221 entries for Veillonella
#> Including 502 entries
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
net@net
```

<img src="02-microbiome_usage_files/figure-html/bsdb_basic_network-1.png" width="768" />

The words-to-species relationship can be plotted by `mbPlot=TRUE`.


```r
net2 <- wcBSDB(c("Veillonella dispar","Neisseria flava",
    "Veillonella parvula","Akkermansia muciniphila"), mbPlot=TRUE,
    curate=TRUE,target="title",pre=TRUE,plotType="network",
    additionalRemove=filter$word,
    numWords=40, corThresh=0.2, colorText=TRUE)
#> Input microbes: 4
#>   Found 17 entries for Veillonella dispar
#>   Found 1 entries for Neisseria flava
#>   Found 20 entries for Veillonella parvula
#>   Found 21 entries for Akkermansia muciniphila
#> Including 90 entries
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
net2@net
```

<img src="02-microbiome_usage_files/figure-html/bsdb_basic_network_mb-1.png" width="768" />

As the BugSigDB contains the relationship between bacterial taxonomy and disease, disease name can also be plotted. When `disPlot=TRUE`, the `mbPlot`
 will be set to `TRUE` by default.
 

```r
net3 <- wcBSDB(c("Veillonella dispar","Neisseria flava",
    "Veillonella parvula","Akkermansia muciniphila"), mbPlot=TRUE,
    curate=TRUE,target="title",pre=TRUE,plotType="network",
    additionalRemove=filter$word, disPlot=TRUE,
    numWords=40, corThresh=0.2, colorText=TRUE)
#> Input microbes: 4
#>   Found 17 entries for Veillonella dispar
#>   Found 1 entries for Neisseria flava
#>   Found 20 entries for Veillonella parvula
#>   Found 21 entries for Akkermansia muciniphila
#> Including 90 entries
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
net3@net
```

<img src="02-microbiome_usage_files/figure-html/bsdb_basic_network_dis-1.png" width="960" />

Other than curated databases, the PubMed query can also be performed with setting `curate=FALSE`. This way, the text information of the latest literature for the microbes and diseases can be plotted. The options for use in function obtaining PubMed information can be specified to `abstArg` in list format, like `sortOrder="pubdate"`.


```r
net4 <- wcBSDB(c("Veillonella dispar","Neisseria flava",
    "Veillonella parvula","Akkermansia muciniphila"), mbPlot=TRUE,
    curate=FALSE,target="title",pre=TRUE,plotType="network",
    additionalRemove=filter$word, disPlot=TRUE,
    numWords=40, corThresh=0.2, colorText=TRUE,
    abstArg = list(retMax=80, sortOrder="pubdate"))
#> Input microbes: 4
#>   Found 17 entries for Veillonella dispar
#>   Found 1 entries for Neisseria flava
#>   Found 20 entries for Veillonella parvula
#>   Found 21 entries for Akkermansia muciniphila
#> Including 90 entries
#> Proceeding without API key
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
net4@net
```

<img src="02-microbiome_usage_files/figure-html/no_curate-1.png" width="960" />

For microbiome analysis, it is often the case that investigating coded enzymes is important. Using `wcEC` function and `getUPtax` function, the queried species or genus can be linked to possible interaction with enzymes using following databases. The downloaded file path should be specified to the function like below to link the queried taxonomy and enzymes. Specifically, enzymes listed in `enzyme.dat` are searched, and corresponding UniProt identifiers are obtained, followed by mapping using `speclist.txt`. This way, the links to microbe - textual information - enzyme can be plotted. 

- [enzyme.dat from Expasy](https://enzyme.expasy.org/)
- [UniProt Controlled vocabulary of species](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist.txt)


```r
# Not run
vp <- wcBSDB(c("Veillonella parvula"),
                 plotType="network",
                curate=TRUE, target="title", edgeLink=FALSE,
                mbPlot = TRUE, ecPlot=TRUE, disPlot=TRUE, tag=TRUE,
                cl=snow::makeCluster(10),colorText=TRUE, pre=TRUE, numWords=30,
                nodePal=RColorBrewer::brewer.pal(10, "Set2"),
                ecFile="../enzyme.dat",
                upTaxFile = "../speclist.txt")
#> Input microbes: 1
#>   Found 20 entries for Veillonella parvula
#> Including 31 entries
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
#> Multiscale bootstrap... Done.
#> Processing EC file
#>   Linking taxonomy to EC
vp@net
```

<img src="02-microbiome_usage_files/figure-html/enzyme-1.png" width="960" />

Further, the relationship between metabolites and microbiome is of interest. Recent studies have revealed various associations in gut microbiome composition and human plasma metabolites, as well as in the other environments.

- Wishart DS, Oler E, Peters H, et al. MiMeDB: the Human Microbial Metabolome Database. Nucleic Acids Res. 2023;51(D1):D611-D620. doi:10.1093/nar/gkac868
- Muller E, Algavi YM, Borenstein E. The gut microbiome-metabolome dataset collection: a curated resource for integrative meta-analysis. npj Biofilms and Microbiomes. 2022;8(1):1-7. doi:10.1038/s41522-022-00345-5
- Dekkers KF, Sayols-Baixeras S, Baldanzi G, et al. An online atlas of human plasma metabolite signatures of gut microbiome composition. Nat Commun. 2022;13(1):5370. doi:10.1038/s41467-022-33050-0

We now use data obtained in [Dekkers et al.](https://www.nature.com/articles/s41467-022-33050-0).
First, read the downloaded file using `readxl`.


```r
metab <- readxl::read_excel(
  "../41467_2022_33050_MOESM8_ESM.xlsx",skip = 7)
```

Pass this tibble, as well as the columns to represent `taxonomy`, `metabolites`, and `quantitative values` to `metabCol`.


```r

metabEx <- wcBSDB(c("Akkermansia muciniphila"),
                edgeLink=FALSE,
                curate=TRUE,
                corThresh=0.3,
                nodePal=RColorBrewer::brewer.pal(10, "Dark2"),
                pre=TRUE, tag=TRUE,
                additionalRemove = filter$word,
                target="abstract", colorText=TRUE,
                plotType="network", numWords=50,
                mbPlot=TRUE,
                metab =metab, metabThresh=0.15,
                preserve = TRUE,
                cl=snow::makeCluster(10),
                abstArg = list(retMax=80,
                               sortOrder="relevance"))
#> Input microbes: 1
#>   Found 21 entries for Akkermansia muciniphila
#> Including 31 entries
#> Target is abstract
#>   Querying PubMed for 20 pmids
#>   Querying without API key
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
#> Multiscale bootstrap... Done.
#> Checking metabolites
metabEx@net
```

<img src="02-microbiome_usage_files/figure-html/metabex-1.png" width="960" />

In this way, we can plot links between microbes - metabolites - textual information.



```r
sessionInfo()
#> R version 4.2.2 (2022-10-31 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 22621)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=Japanese_Japan.utf8 
#> [2] LC_CTYPE=Japanese_Japan.utf8   
#> [3] LC_MONETARY=Japanese_Japan.utf8
#> [4] LC_NUMERIC=C                   
#> [5] LC_TIME=Japanese_Japan.utf8    
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets 
#> [6] methods   base     
#> 
#> other attached packages:
#> [1] RColorBrewer_1.1-3   ggraph_2.1.0        
#> [3] ggplot2_3.4.0        wcGeneSummary_0.99.0
#> 
#> loaded via a namespace (and not attached):
#>   [1] GeneSummary_0.99.4     colorspace_2.0-3      
#>   [3] rjson_0.2.21           ISOcodes_2022.09.29   
#>   [5] ellipsis_0.3.2         XVector_0.38.0        
#>   [7] GlobalOptions_0.1.2    base64enc_0.1-3       
#>   [9] ggdendro_0.1.23        fs_1.5.2              
#>  [11] rstudioapi_0.14        farver_2.1.1          
#>  [13] graphlayouts_0.8.4     ggrepel_0.9.2         
#>  [15] bit64_4.0.5            AnnotationDbi_1.60.0  
#>  [17] fansi_1.0.3            xml2_1.3.3            
#>  [19] codetools_0.2-18       downlit_0.4.2         
#>  [21] cachem_1.0.6           knitr_1.41            
#>  [23] polyclip_1.10-4        jsonlite_1.8.4        
#>  [25] dbplyr_2.3.0           png_0.1-8             
#>  [27] graph_1.76.0           ggforce_0.4.1         
#>  [29] shiny_1.7.4            bugsigdbr_1.5.3       
#>  [31] rentrez_1.2.3          compiler_4.2.2        
#>  [33] httr_1.4.4             assertthat_0.2.1      
#>  [35] fastmap_1.1.0          cli_3.6.0             
#>  [37] later_1.3.0            tweenr_2.0.2          
#>  [39] htmltools_0.5.4        tools_4.2.2           
#>  [41] igraph_1.3.5           NLP_0.2-1             
#>  [43] gtable_0.3.1           glue_1.6.2            
#>  [45] GenomeInfoDbData_1.2.9 dplyr_1.0.10          
#>  [47] rappdirs_0.3.3         Rcpp_1.0.9            
#>  [49] slam_0.1-50            Biobase_2.58.0        
#>  [51] cellranger_1.1.0       jquerylib_0.1.4       
#>  [53] vctrs_0.5.1            Biostrings_2.66.0     
#>  [55] xfun_0.36              stringr_1.5.0         
#>  [57] stopwords_2.3          mime_0.12             
#>  [59] lifecycle_1.0.3        pvclust_2.2-0         
#>  [61] XML_3.99-0.13          dendextend_1.16.0     
#>  [63] org.Hs.eg.db_3.16.0    zlibbioc_1.44.0       
#>  [65] MASS_7.3-58.1          scales_1.2.1          
#>  [67] tidygraph_1.2.2        promises_1.2.0.1      
#>  [69] parallel_4.2.2         cyjShiny_1.0.34       
#>  [71] curl_5.0.0             yaml_2.3.6            
#>  [73] memoise_2.0.1          gridExtra_2.3         
#>  [75] yulab.utils_0.0.6      sass_0.4.4            
#>  [77] stringi_1.7.12         RSQLite_2.2.20        
#>  [79] highr_0.10             S4Vectors_0.36.1      
#>  [81] filelock_1.0.2         BiocGenerics_0.44.0   
#>  [83] GenomeInfoDb_1.34.6    rlang_1.0.6           
#>  [85] pkgconfig_2.0.3        bitops_1.0-7          
#>  [87] evaluate_0.19          purrr_1.0.1           
#>  [89] labeling_0.4.2         patchwork_1.1.2       
#>  [91] htmlwidgets_1.6.1      cowplot_1.1.1         
#>  [93] bit_4.0.5              tidyselect_1.2.0      
#>  [95] magrittr_2.0.3         bookdown_0.31         
#>  [97] R6_2.5.1               snow_0.4-4            
#>  [99] IRanges_2.32.0         generics_0.1.3        
#> [101] DBI_1.1.3              pillar_1.8.1          
#> [103] withr_2.5.0            KEGGREST_1.38.0       
#> [105] RCurl_1.98-1.9         tibble_3.1.8          
#> [107] crayon_1.5.2           wordcloud_2.6         
#> [109] utf8_1.2.2             BiocFileCache_2.6.0   
#> [111] rmarkdown_2.19         viridis_0.6.2         
#> [113] GetoptLong_1.0.5       readxl_1.4.1          
#> [115] grid_4.2.2             blob_1.2.3            
#> [117] digest_0.6.31          xtable_1.8-4          
#> [119] tm_0.7-10              tidyr_1.2.1           
#> [121] httpuv_1.6.8           gridGraphics_0.5-1    
#> [123] stats4_4.2.2           munsell_0.5.0         
#> [125] ggplotify_0.1.0        viridisLite_0.4.1     
#> [127] bslib_0.4.2
```
