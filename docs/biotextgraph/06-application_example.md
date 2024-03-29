

# Application examples {#app}

## Actual application of the package

### The selection of genes

The input of genes can be selected in various ways. In most cases, the list obtained by differential expression analysis can be used as an input to the function. Additionally, the cluster of genes identified in gene clustering analysis such as weighted gene co-expression network analysis can be used as input.

## Gene cluster

### Load the necessary packages and genes


```r
library(biotextgraph)
library(org.Hs.eg.db)
library(ggplot2)
library(ReactomePA);library(clusterProfiler)
library(ggraph);library(igraph)
load(system.file("extdata", "sysdata.rda", package = "biotextgraph"))
```

For the omics analysis involving transcript or genes, we would not obtain a list containing the single gene type such as ERCC genes shown in the basic usage example. In almost all the cases, the gene comes from various biological pathways like obtained in the analysis mentioned in the above section. Here, we introduce an example using the example cluster of WGCNA analysis obtained from the transcriptomic dataset investigating bladder cancer ([Chen et al. 2019](https://doi.org/10.1038/s41556-019-0361-y)). The gene set contains 621 genes.

First, we perform over-representation analysis on the gene set (KEGG) to grasp the biological functions of these genes.


```r
converted <- clusterProfiler::bitr(exampleWGCNAcluster,
    fromType="ENSEMBL", toType="ENTREZID",
    OrgDb="org.Hs.eg.db")[,2]
ora <- clusterProfiler::enrichKEGG(converted)
ora |> dplyr::filter(p.adjust<0.05) |> data.frame() |> dplyr::pull("Description")
#> [1] "p53 signaling pathway"                
#> [2] "Pathogenic Escherichia coli infection"
#> [3] "PI3K-Akt signaling pathway"           
#> [4] "ECM-receptor interaction"             
#> [5] "Small cell lung cancer"               
#> [6] "Lipid and atherosclerosis"            
#> [7] "Hepatocellular carcinoma"             
#> [8] "Influenza A"
enrichplot::dotplot(ora, showCategory=20)
```

<img src="06-application_example_files/figure-html/aa_kegg-1.png" width="100%" style="display: block; margin: auto;" />

We use the function in `biotextgraph` to make a summarized visualization of textual information, along with associated genes. As for the input reaching to a hundred of genes, there is an option `filterByGO` term, which filters the text mining results to those used in GO terms. This is useful for limiting the visualization to biologically-relevant terms.


```r
check <- refseq(converted, genePlot=TRUE,
                filterByGO=TRUE, keyType="ENTREZID", autoThresh=FALSE)
#> Input genes: 168
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> `filterByGO` option enabled. Filtering by GO terms ...
plotNet(check, asis=TRUE)
```

<img src="06-application_example_files/figure-html/aa_btg-1.png" width="100%" style="display: block; margin: auto;" />

Some genes are not related to the significantly enriched pathway, and one would like to inspect the biological function of these genes.


```r
## Extraction of non-related genes
enr_genes <- ora@result %>% data.frame() %>% 
  filter(p.adjust<0.05) %>% dplyr::pull(geneID) %>%
  strsplit("/") %>% unlist() %>% unique()

no_enr <- converted[!(converted %in% enr_genes)]
length(no_enr)
#> [1] 135
```

This time, we enable the `filterByGO` option along with `ngram=2`, which produces 2-gram visualization.


```r
check_noenr_WGO <- refseq(no_enr,
    layout="nicely",
    ngram=2,
    keyType="ENTREZID",
    filterByGO=TRUE,
    autoThresh=FALSE,
    corThresh=0.1,
    docsum=TRUE)
#> Input genes: 135
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> `filterByGO` option enabled. Filtering by GO terms ...
plotNet(check_noenr_WGO, asis=TRUE)
```

<img src="06-application_example_files/figure-html/aa_noenr-1.png" width="100%" style="display: block; margin: auto;" />

This way, we can filter the unnecessary words from the many terms included in the textual information of identifiers, and focus on biologically relevant terms contained in the identifiers.

## DEGs

We demonstrate an use case of the package, which investigates transcriptomic changes induced by BK polyomavirus (BKPyV) infection in renal proximal tubular epithelial cells ([Assetta et al. 2016](https://pubmed.ncbi.nlm.nih.gov/27381292/)). Differentially expressed mRNAs in 3 days post-infection were obtained, and down-regulated mRNAs in BKPyV infected cells were examined.

### Load the necessary packages and genes


```r
library(biotextgraph)
library(org.Hs.eg.db)
library(ggplot2)
library(ReactomePA);library(clusterProfiler)
library(ggraph);library(igraph)
library(ggforce) ## For genePathPlot
load(system.file("extdata", "sysdata.rda", package = "biotextgraph"))
degs <- d3degDownAssetta2016
length(degs)
#> [1] 191
degs
#>   [1] "ABCB4"           "ABCB7"           "AKTIP"          
#>   [4] "ALS2"            "ANKRA2"          "ANTXR1"         
#>   [7] "APH1B"           "ARHGEF28"        "ARNTL"          
#>  [10] "ATMIN"           "BDH2"            "BIRC2"          
#>  [13] "BTG1"            "BTG2"            "C4A"            
#>  [16] "C7orf60"         "C8orf4"          "CALCOCO1"       
#>  [19] "CAPRIN2"         "CARS"            "CBLB"           
#>  [22] "CCNDBP1"         "CCNG2"           "CCZ1B"          
#>  [25] "CDC42EP3"        "CDH6"            "CFLAR"          
#>  [28] "CHMP1B"          "CHMP4C"          "CLUAP1"         
#>  [31] "COG2"            "COG3"            "CPQ"            
#>  [34] "CROT"            "CTTNBP2NL"       "CYP4V2"         
#>  [37] "DAB2"            "DDB2"            "DDX17"          
#>  [40] "DDX5"            "DGKA"            "DHX32"          
#>  [43] "DLG1"            "DYNC2LI1"        "DYRK2"          
#>  [46] "EFHC1"           "EIF4A2"          "ERMAP"          
#>  [49] "ERMARD"          "EXOC1"           "FAM134B"        
#>  [52] "FAM160B1"        "FAM21C"          "FAM84B"         
#>  [55] "FANK1"           "FAS"             "FBXO38"         
#>  [58] "FCHO2"           "FGD6"            "FLJ22447"       
#>  [61] "FMNL2"           "GADD45A"         "GJA1"           
#>  [64] "GLIDR"           "GLT8D1"          "GOPC"           
#>  [67] "GPBP1L1"         "GPR155"          "GPR75-ASB3"     
#>  [70] "GRAMD3"          "HADHB"           "HCG11"          
#>  [73] "HDAC9"           "HDHD2"           "HERPUD1"        
#>  [76] "HSDL2"           "ICA1"            "ICK"            
#>  [79] "IFNGR1"          "IFT46"           "IFT74"          
#>  [82] "IRF6"            "ITGA2"           "ITGA6"          
#>  [85] "ITGAV"           "ITGB6"           "KDM5B"          
#>  [88] "KIF3A"           "KIF5B"           "KLHL20"         
#>  [91] "KLHL24"          "KLHL9"           "KPNA5"          
#>  [94] "KRCC1"           "L3MBTL3"         "LINC00657"      
#>  [97] "LOC100131564"    "LZTFL1"          "MAMDC2"         
#> [100] "MAP4K5"          "MAT2B"           "MBNL2"          
#> [103] "MDM2"            "MECOM"           "MFSD1"          
#> [106] "MGEA5"           "MICU3"           "MSANTD4"        
#> [109] "NBPF11"          "NCBP2"           "NEAT1"          
#> [112] "NTM"             "OGT"             "PAFAH1B2"       
#> [115] "PAFAH2"          "PCMTD2"          "PDE4D"          
#> [118] "PDP1"            "PELI1"           "PEX1"           
#> [121] "PHF14"           "PHOSPHO2-KLHL23" "PIK3IP1"        
#> [124] "PLA2R1"          "PLCB4"           "PLK2"           
#> [127] "POLI"            "POSTN"           "PPAN-P2RY11"    
#> [130] "PPFIBP1"         "PPP2CB"          "PRICKLE1"       
#> [133] "PROS1"           "PSMD5-AS1"       "RAD17"          
#> [136] "RHOQ"            "RIMKLB"          "RNA18S5"        
#> [139] "RNF170"          "RNF20"           "RNU1-28P"       
#> [142] "RPL23AP53"       "RRM2B"           "RRN3"           
#> [145] "RRN3P1"          "SEMA3C"          "SERINC1"        
#> [148] "SESN1"           "SESN3"           "SGK1"           
#> [151] "SLC22A5"         "SLC37A3"         "SVIL"           
#> [154] "SYT11"           "TARSL2"          "TBC1D19"        
#> [157] "TBCK"            "TBRG1"           "TGFA"           
#> [160] "TGFB2"           "TIPARP"          "TMEM136"        
#> [163] "TNFRSF10B"       "TNFRSF10D"       "TOM1L1"         
#> [166] "TP53INP1"        "TRIM13"          "TRIM32"         
#> [169] "TRIM4"           "TSC1"            "TSPYL5"         
#> [172] "UGT2B7"          "UNC13B"          "UPRT"           
#> [175] "VPS41"           "VPS8"            "WDR11"          
#> [178] "WDR19"           "XPC"             "YPEL2"          
#> [181] "ZC2HC1A"         "ZFAND5"          "ZFP90"          
#> [184] "ZFR"             "ZMAT3"           "ZNF12"          
#> [187] "ZNF248"          "ZNF322"          "ZNF561"         
#> [190] "ZNF626"          "ZSCAN30"
set.seed(1)
```

### Enrichment analysis

First, we perform enrichment analysis using ReactomePA.
From the enrichment analysis results, the cluster is related to transcriptional regulation by TP53.


```r
## Convert to ENTREZID
entre <- AnnotationDbi::select(org.Hs.eg.db, keytype="SYMBOL",
                               keys = degs, columns = "ENTREZID")$ENTREZID
pway <- setReadable(enrichPathway(entre), org.Hs.eg.db)
sigpway <- subset(pway |> data.frame(), p.adjust<0.05)
sigpway$Description
#>  [1] "Transcriptional Regulation by TP53"                              
#>  [2] "RIPK1-mediated regulated necrosis"                               
#>  [3] "Regulation of necroptotic cell death"                            
#>  [4] "Regulated Necrosis"                                              
#>  [5] "Intraflagellar transport"                                        
#>  [6] "Regulation by c-FLIP"                                            
#>  [7] "CASP8 activity is inhibited"                                     
#>  [8] "Dimerization of procaspase-8"                                    
#>  [9] "TP53 Regulates Transcription of Death Receptors and Ligands"     
#> [10] "Caspase activation via Death Receptors in the presence of ligand"
#> [11] "FOXO-mediated transcription of cell cycle genes"                 
#> [12] "Cilium Assembly"                                                 
#> [13] "TP53 Regulates Transcription of Cell Death Genes"
cnetplot(pway)
```

<img src="06-application_example_files/figure-html/ea-1.png" width="100%" style="display: block; margin: auto;" />

```r

## Genes involved in significant pathways
excheck <- unlist(unique(sapply(sigpway$geneID,
                                function (x) strsplit(x,"/"))))
length(excheck)
#> [1] 57

## Not involved in significant pathways
no_enr <- degs[!(degs %in% excheck)]
length(no_enr)
#> [1] 159
```

### Text mining the gene summaries

Next we perform the plain function producing a correlation network, with showing the top-genes related to high-frequency words in the text in RefSeq summary. We make two networks, those related to significantly enriched pathways and those not related to the significant pathways.


```r
net1 <- refseq(excheck, colorText=TRUE, genePlot=TRUE, numWords=50, docsum=TRUE)
#> Input genes: 57
#>   Converted input genes: 32
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> Ignoring corThresh, automatically determine the value
#> threshold = 0.3
net2 <- refseq(no_enr, colorText=TRUE, genePlot=TRUE, numWords=50, docsum=TRUE)
#> Input genes: 159
#>   Converted input genes: 142
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> Ignoring corThresh, automatically determine the value
#> threshold = 0.1

plotNet(net2, asis=TRUE)
```

<img src="06-application_example_files/figure-html/basic2app-1.png" width="100%" style="display: block; margin: auto;" />

### Combine and inspect merged network

These networks can be combined to find intersections and differences, from textual perspective. We can see from the results that the genes related to non-enriched pathways could be associated with the membrane transport and adhesion, while genes in the enriched pathways have textual information resembling to those in curated biological databases (such as DNA damage response and tumor suppression).


```r
compareWordNet(list(net1, net2),
               titles=c("Enrich","No_Enrich")) |>
plotNet()
```

<img src="06-application_example_files/figure-html/combineapp-1.png" width="100%" style="display: block; margin: auto;" />



```r
sessionInfo()
#> R version 4.3.1 (2023-06-16 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 11 x64 (build 22621)
#> 
#> Matrix products: default
#> 
#> 
#> locale:
#> [1] LC_COLLATE=Japanese_Japan.utf8 
#> [2] LC_CTYPE=Japanese_Japan.utf8   
#> [3] LC_MONETARY=Japanese_Japan.utf8
#> [4] LC_NUMERIC=C                   
#> [5] LC_TIME=Japanese_Japan.utf8    
#> 
#> time zone: Asia/Tokyo
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils    
#> [6] datasets  methods   base     
#> 
#> other attached packages:
#>  [1] ggforce_0.4.1         igraph_1.5.1         
#>  [3] ggraph_2.1.0.9000     clusterProfiler_4.9.5
#>  [5] ReactomePA_1.46.0     org.Hs.eg.db_3.17.0  
#>  [7] AnnotationDbi_1.63.2  IRanges_2.35.2       
#>  [9] S4Vectors_0.38.1      Biobase_2.61.0       
#> [11] BiocGenerics_0.47.0   biotextgraph_0.99.0  
#> [13] ggplot2_3.4.4        
#> 
#> loaded via a namespace (and not attached):
#>   [1] splines_4.3.1                
#>   [2] later_1.3.1                  
#>   [3] bitops_1.0-7                 
#>   [4] ggplotify_0.1.2              
#>   [5] filelock_1.0.2               
#>   [6] tibble_3.2.1                 
#>   [7] polyclip_1.10-4              
#>   [8] graph_1.79.1                 
#>   [9] XML_3.99-0.14                
#>  [10] lifecycle_1.0.3              
#>  [11] bugsigdbr_1.8.1              
#>  [12] NLP_0.2-1                    
#>  [13] lattice_0.21-8               
#>  [14] MASS_7.3-60                  
#>  [15] dendextend_1.17.1            
#>  [16] magrittr_2.0.3               
#>  [17] sass_0.4.7                   
#>  [18] rmarkdown_2.25               
#>  [19] jquerylib_0.1.4              
#>  [20] yaml_2.3.7                   
#>  [21] httpuv_1.6.11                
#>  [22] cowplot_1.1.1                
#>  [23] DBI_1.1.3                    
#>  [24] RColorBrewer_1.1-3           
#>  [25] zlibbioc_1.47.0              
#>  [26] purrr_1.0.2                  
#>  [27] downlit_0.4.3                
#>  [28] RCurl_1.98-1.12              
#>  [29] yulab.utils_0.1.0            
#>  [30] tweenr_2.0.2                 
#>  [31] rappdirs_0.3.3               
#>  [32] pvclust_2.2-0                
#>  [33] GenomeInfoDbData_1.2.10      
#>  [34] ISOcodes_2022.09.29          
#>  [35] enrichplot_1.21.3            
#>  [36] cyjShiny_1.0.42              
#>  [37] tm_0.7-11                    
#>  [38] ggrepel_0.9.5                
#>  [39] tidytree_0.4.5               
#>  [40] rentrez_1.2.3                
#>  [41] reactome.db_1.86.0           
#>  [42] codetools_0.2-19             
#>  [43] DOSE_3.27.2                  
#>  [44] xml2_1.3.5                   
#>  [45] tidyselect_1.2.0             
#>  [46] aplot_0.2.1                  
#>  [47] farver_2.1.1                 
#>  [48] viridis_0.6.4                
#>  [49] GeneSummary_0.99.6           
#>  [50] BiocFileCache_2.9.1          
#>  [51] base64enc_0.1-3              
#>  [52] jsonlite_1.8.7               
#>  [53] GetoptLong_1.0.5             
#>  [54] ellipsis_0.3.2               
#>  [55] tidygraph_1.2.3              
#>  [56] tools_4.3.1                  
#>  [57] treeio_1.25.4                
#>  [58] HPO.db_0.99.2                
#>  [59] Rcpp_1.0.11                  
#>  [60] glue_1.6.2                   
#>  [61] gridExtra_2.3                
#>  [62] xfun_0.40                    
#>  [63] qvalue_2.33.0                
#>  [64] GenomeInfoDb_1.37.4          
#>  [65] dplyr_1.1.4                  
#>  [66] withr_2.5.0                  
#>  [67] BiocManager_1.30.22          
#>  [68] fastmap_1.1.1                
#>  [69] fansi_1.0.4                  
#>  [70] digest_0.6.33                
#>  [71] R6_2.5.1                     
#>  [72] mime_0.12                    
#>  [73] gridGraphics_0.5-1           
#>  [74] colorspace_2.1-0             
#>  [75] GO.db_3.17.0                 
#>  [76] RSQLite_2.3.1                
#>  [77] utf8_1.2.3                   
#>  [78] tidyr_1.3.0                  
#>  [79] generics_0.1.3               
#>  [80] data.table_1.14.8            
#>  [81] graphlayouts_1.0.0           
#>  [82] stopwords_2.3                
#>  [83] httr_1.4.7                   
#>  [84] htmlwidgets_1.6.2            
#>  [85] scatterpie_0.2.1             
#>  [86] graphite_1.48.0              
#>  [87] pkgconfig_2.0.3              
#>  [88] gtable_0.3.4                 
#>  [89] blob_1.2.4                   
#>  [90] XVector_0.41.1               
#>  [91] shadowtext_0.1.2             
#>  [92] htmltools_0.5.6              
#>  [93] bookdown_0.35                
#>  [94] fgsea_1.27.1                 
#>  [95] scales_1.3.0                 
#>  [96] png_0.1-8                    
#>  [97] wordcloud_2.6                
#>  [98] ggfun_0.1.3                  
#>  [99] ggdendro_0.1.23              
#> [100] knitr_1.44                   
#> [101] rstudioapi_0.15.0            
#> [102] reshape2_1.4.4               
#> [103] rjson_0.2.21                 
#> [104] nlme_3.1-163                 
#> [105] curl_5.0.2                   
#> [106] cachem_1.0.8                 
#> [107] GlobalOptions_0.1.2          
#> [108] stringr_1.5.0                
#> [109] BiocVersion_3.18.0           
#> [110] parallel_4.3.1               
#> [111] HDO.db_0.99.1                
#> [112] pillar_1.9.0                 
#> [113] grid_4.3.1                   
#> [114] vctrs_0.6.5                  
#> [115] slam_0.1-50                  
#> [116] promises_1.2.1               
#> [117] dbplyr_2.3.3                 
#> [118] xtable_1.8-4                 
#> [119] evaluate_0.21                
#> [120] cli_3.6.1                    
#> [121] compiler_4.3.1               
#> [122] rlang_1.1.1                  
#> [123] crayon_1.5.2                 
#> [124] labeling_0.4.3               
#> [125] plyr_1.8.8                   
#> [126] fs_1.6.3                     
#> [127] stringi_1.7.12               
#> [128] viridisLite_0.4.2            
#> [129] BiocParallel_1.35.4          
#> [130] MPO.db_0.99.7                
#> [131] munsell_0.5.0                
#> [132] Biostrings_2.69.2            
#> [133] lazyeval_0.2.2               
#> [134] GOSemSim_2.27.3              
#> [135] Matrix_1.6-5                 
#> [136] patchwork_1.2.0              
#> [137] bit64_4.0.5                  
#> [138] KEGGREST_1.41.0              
#> [139] shiny_1.7.5                  
#> [140] interactiveDisplayBase_1.39.0
#> [141] AnnotationHub_3.9.2          
#> [142] memoise_2.0.1                
#> [143] bslib_0.5.1                  
#> [144] ggtree_3.9.1                 
#> [145] fastmatch_1.1-4              
#> [146] bit_4.0.5                    
#> [147] ape_5.7-1                    
#> [148] gson_0.1.0
```
