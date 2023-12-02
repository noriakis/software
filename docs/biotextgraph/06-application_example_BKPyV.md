

# Application examples {#app}

We demonstrate an use case of the package, which investigates transcriptomic changes induced by BK polyomavirus (BKPyV) infection in renal proximal tubular epithelial cells ([Assetta et al. 2016](https://pubmed.ncbi.nlm.nih.gov/27381292/)). Differentially expressed mRNAs in 3 days post-infection were obtained, and down-regulated mRNAs in BKPyV infected cells were examined.

## Load the necessary packages and genes


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

## Enrichment analysis

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

<img src="06-application_example_BKPyV_files/figure-html/ea-1.png" width="100%" style="display: block; margin: auto;" />

```r

## Genes involved in significant pathways
excheck <- unlist(unique(sapply(sigpway$geneID,
                                function (x) strsplit(x,"/"))))
excheck
#>  [1] "BTG2"      "DDB2"      "DYRK2"     "FAS"      
#>  [5] "GADD45A"   "MDM2"      "PLK2"      "PPP2CB"   
#>  [9] "RAD17"     "RRM2B"     "SESN1"     "SESN3"    
#> [13] "SGK1"      "TNFRSF10B" "TNFRSF10D" "TP53INP1" 
#> [17] "TSC1"      "BIRC2"     "CFLAR"     "FAS"      
#> [21] "OGT"       "PELI1"     "TNFRSF10B" "BIRC2"    
#> [25] "CFLAR"     "CHMP4C"    "FAS"       "OGT"      
#> [29] "PELI1"     "TNFRSF10B" "CLUAP1"    "DYNC2LI1" 
#> [33] "IFT46"     "IFT74"     "KIF3A"     "WDR19"    
#> [37] "CFLAR"     "FAS"       "TNFRSF10B" "FAS"      
#> [41] "TNFRSF10B" "TNFRSF10D" "BTG1"      "CCNG2"    
#> [45] "GADD45A"   "CLUAP1"    "DYNC2LI1"  "EXOC1"    
#> [49] "IFT46"     "IFT74"     "KIF3A"     "LZTFL1"   
#> [53] "WDR19"     "FAS"       "TNFRSF10B" "TNFRSF10D"
#> [57] "TP53INP1"
```

We store the name of enriched pathways in the network for the downstream analysis.


```r
netreac <- refseq(degs,
            enrich="reactome",
            plotType="network",
            numWords=50,
            colorText=TRUE)
#> Input genes: 191
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 174
#> Performing enrichment analysis
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> Scale for size is already present.
#> Adding another scale for size, which will replace the
#> existing scale.
```

## Text mining the gene summaries
Next we perform the plain function producing a correlation network, with showing the top-genes related to high-frequency words in the text in RefSeq summary. We obtained the list of these genes from geneCount slot.


```r
net1 <- refseq(excheck,
                plotType="network",
				colorText=TRUE,
				numWords=30,
				corThresh=0.5,
				genePlot=TRUE,
				genePlotNum=5,
				edgeLink=FALSE,
				genePathPlot="reactome")
#> Input genes: 57
#>   Converted input genes: 32
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> Found 28 enriched term
plotNet(net1)
```

<img src="06-application_example_BKPyV_files/figure-html/basic2app-1.png" width="100%" style="display: block; margin: auto;" />

```r
gnc <- getSlot(net1, "geneCount")
gnc
#> 
#>     RAD17     SESN1      TSC1      DDB2 TNFRSF10B      BTG1 
#>         9         9         6         5         5         4 
#>     DYRK2     BIRC2     CCNG2       FAS   GADD45A    LZTFL1 
#>         4         3         3         3         3         3 
#>      MDM2     SESN3      SGK1     CFLAR    CHMP4C       OGT 
#>         3         3         3         2         2         2 
#> TNFRSF10D     WDR19      BTG2  DYNC2LI1     EXOC1     IFT74 
#>         2         2         1         1         1         1 
#>    PPP2CB     RRM2B 
#>         1         1

top <- names(gnc[gnc>=5])
top
#> [1] "RAD17"     "SESN1"     "TSC1"      "DDB2"     
#> [5] "TNFRSF10B"
```

## Text mining the available literature

These genes are further queried for PubMed information.
The results were obtained in April, 2023.
First, show the network for the article titles.


```r
# titlenet <- pubmed(top,
#                    sortOrder="relevance",
#                    target="title",
#                    plotType="network",
#                    colorText=TRUE,
#                    madeUpperGenes=FALSE,
#                    corThresh=0.2,
#                    pre=TRUE,
#                    retMax=40,
#                    numWords=40,
#                    edgeLink=FALSE)
load("titlenet_20230420.rda")
titlenet |> plotNet()
```

<img src="06-application_example_BKPyV_files/figure-html/absttitleapp-1.png" width="100%" style="display: block; margin: auto;" />

Obtain and show the network for the article abstract.


```r
# abstnet <- pubmed(top,
#                    target="abstract",
#                    plotType="network",
#                    colorText=TRUE,
#                    madeUpperGenes=FALSE,
#                    corThresh=0.2,
#                    pre=TRUE,
#                    retMax=40,
#                    numWords=40,
#                    edgeLink=FALSE)
load("abstnet_20230420.rda")
abstnet |> plotNet()
```

<img src="06-application_example_BKPyV_files/figure-html/abstmainapp-1.png" width="100%" style="display: block; margin: auto;" />

## Combine and inspect merged network

From the RefSeq summary and articles related to important genes, the cluster could have functionality of DNA damage response, which is also upregulated by BKPyV infection. These networks can be combined to find intersections and differences. We can see that in addition to Reactome pathway names, plenty of information could be obtained and summarized by querying other databases, which could aid in interpreting clusters of genes and hypothesis generation.


```r
compareWordNet(list(abstnet, titlenet, netreac, net1),
               titles=c("Abstract","Title","Reactome","RefSeq")) |> plotNet()
```

<img src="06-application_example_BKPyV_files/figure-html/combineapp-1.png" width="100%" style="display: block; margin: auto;" />

From the network, DNA damage repair pathway, especially nucleotide excision repair related to DDB2, METTL14,and RAD17 might be related to BKPyV infection, which cannot be prioritize based on log2FoldChange or enrichment analysis.


```r
# btg <- compareWordNet(list(net1,titlenet,abstnet,netreac))
load("btg_20230420.rda")
conet <- btg |> getSlot("igraphRaw")

ddrNms <- NULL
for (nm in names(V(conet))) {
    if (tolower(nm) %in% c("dna","damage","repair")) {
        ddrNms <- c(ddrNms, names(neighbors(conet, nm)))
    }
}
ddrNms
#>  [1] "checkpoint"      "cycle"           "DNA"            
#>  [4] "phosphorylation" "required"        "response"       
#>  [7] "SESN1"           "DDB2"            "RAD17"          
#> [10] "checkpoint"      "cycle"           "damage"         
#> [13] "phosphorylation" "required"        "response"       
#> [16] "SESN1"           "DDB2"            "RAD17"          
#> [19] "Cancer"          "degradation"     "protein"        
#> [22] "regulates"       "Repair"          "CDT2"           
#> [25] "chromatin"       "complex"         "Damage"         
#> [28] "function"        "induced"         "DNA"            
#> [31] "DDB2"            "Excision"        "Human"          
#> [34] "mice"            "Nucleotide"      "protein"        
#> [37] "regulates"       "level"           "METTL14"        
#> [40] "checkpoint"      "DNA"             "response"       
#> [43] "RAD17"           "Cancer"          "chromatin"      
#> [46] "complex"         "function"

ddrRelated <- induced.subgraph(conet,
                 names(V(conet)) %in% unique(ddrNms))
V(ddrRelated)$SYMBOL <- names(V(ddrRelated)) %in% keys(org.Hs.eg.db,"SYMBOL")
ggraph(ddrRelated)+
    geom_edge_diagonal2(color="grey80")+
    geom_node_point(aes(color=SYMBOL))+
    geom_node_text(aes(label=name, color=SYMBOL),check_overlap=TRUE, repel=TRUE,
                   bg.color = "white", segment.color="black",
                   bg.r = .15)+
    scale_color_manual(values=c("steelblue","tomato"))+
    theme_graph()
```

<img src="06-application_example_BKPyV_files/figure-html/ddr-1.png" width="100%" style="display: block; margin: auto;" />


The network can be obtained by `returnNet=TRUE`, or is stored in returned `biotext` class object in the slot `igraphRaw`, which can be used for downstream analysis like assessment of degrees and graph-based clustering.


```r

conetDeg <- igraph::degree(conet)
conetDeg[order(conetDeg, decreasing=TRUE)] |> head(15)
#>      SESN1       DDB2        DNA      RAD17      Human 
#>         26         21         20         18         15 
#>     Cancer checkpoint   response       TP53  regulates 
#>         14         12         12         12         11 
#>    induced activation    protein     Repair     damage 
#>         11         10         10         10          9

## Graph-based clustering by igraph::cluster_fast_greedy()
## Filling the slot `communities` and assign `igraphRaw` slot `community` attribute.
btg <- btg |> graph_cluster(func=igraph::cluster_fast_greedy)

## Check the membership for merged network
data.frame(getSlot(btg, "communities")$membership,
           getSlot(btg, "communities")$names) |>
  `colnames<-`(c("membership","word")) |>
  group_by(membership) |> summarise(words=paste(word,collapse=",")) |>
  kableExtra::kable()
```

<table>
 <thead>
  <tr>
   <th style="text-align:right;"> membership </th>
   <th style="text-align:left;"> words </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> apoptosis,caspase,death,factor,TNFRSF10B,Apoptosis,carcinogenesis,Lung,promotes,Regulation,Roles </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> cellular,checkpoint,cycle,damage,DNA,phosphorylation,required,response,DDB2,RAD17,Cancer,degradation,Excision,Muscle,Nucleotide,protein,regulates,Repair,CDT2,chromatin,complex,Damage,function,Knockdown,level,METTL14,role </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> functions,growth,interacts,kinase,p53,suppressor,tumor,SESN1,TSC1,Human,Cisplatin,effects,HNSCC,increased,induced,levels,Macrophages,regulation,revealed,Study,Transcription,treatment </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> TP53,neuroblastoma,paediatric,Caspase,Death,ligand,Ligands,necroptotic,presence,Regulates,Transcriptional </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:left;"> Carcinoma,epithelialmesenchymal,head,miR3773p,Neck,Squamous,Transition </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:left;"> DNAbinding,irradiation,mice,protects </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> activation,density,inflammation,lipoproteininduced,macrophages,oxidized,Role,Sestrin1 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:left;"> Necrosis,Regulated,RIPK1mediated </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:left;"> activity,CASP8,inhibited </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:left;"> Intraflagellar,transport </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:left;"> Dimerization,procaspase8 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:left;"> acid,cyclin </td>
  </tr>
</tbody>
</table>



```r

ggraph(btg |> getSlot("igraphRaw"))+
    geom_edge_link(color="grey80", width=0.5)+
    geom_node_point(aes(color=community), size=2)+
    graphhighlight::highlight_node(glow=TRUE, filter=community==2,glow_base_size = TRUE, glow_size=0.5)+
    geom_node_text(aes(label=name, color=community),
                   check_overlap=TRUE, repel=TRUE,
                   bg.color = "white", segment.color="black",
                   bg.r = .15)+
    theme_graph()
```

<img src="06-application_example_BKPyV_files/figure-html/combineNet-1.png" width="100%" style="display: block; margin: auto;" />

Dynamic layout can be also used to compare the networks, by `graphlayouts`, for comparing the multiple graphs, especially useful for time-series analysis. See the documentation of [`layout_as_dynamic`](http://graphlayouts.schochastics.net/reference/layout_dynamic.html) for specifying the alpha, which is default to 0.5. 


```r
library(igraph)
dyn <- plotDynamic(list(abstnet, titlenet), concat="union",
                   titles=c("Abstract","Title"))
#> This graph was created by an old(er) igraph version.
#>   Call upgrade_graph() on it to use with the current igraph version
#>   For now we convert it on the fly...
#> This graph was created by an old(er) igraph version.
#>   Call upgrade_graph() on it to use with the current igraph version
#>   For now we convert it on the fly...
#> This graph was created by an old(er) igraph version.
#>   Call upgrade_graph() on it to use with the current igraph version
#>   For now we convert it on the fly...
#> This graph was created by an old(er) igraph version.
#>   Call upgrade_graph() on it to use with the current igraph version
#>   For now we convert it on the fly...
#> This graph was created by an old(er) igraph version.
#>   Call upgrade_graph() on it to use with the current igraph version
#>   For now we convert it on the fly...
#> This graph was created by an old(er) igraph version.
#>   Call upgrade_graph() on it to use with the current igraph version
#>   For now we convert it on the fly...
#> This graph was created by an old(er) igraph version.
#>   Call upgrade_graph() on it to use with the current igraph version
#>   For now we convert it on the fly...
#> This graph was created by an old(er) igraph version.
#>   Call upgrade_graph() on it to use with the current igraph version
#>   For now we convert it on the fly...
dyn
```

<img src="06-application_example_BKPyV_files/figure-html/dyn-1.png" width="100%" style="display: block; margin: auto;" />



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
#> [13] ggplot2_3.4.3        
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3           
#>   [2] ggdendro_0.1.23              
#>   [3] jsonlite_1.8.7               
#>   [4] pvclust_2.2-0                
#>   [5] magrittr_2.0.3               
#>   [6] farver_2.1.1                 
#>   [7] rmarkdown_2.25               
#>   [8] GlobalOptions_0.1.2          
#>   [9] fs_1.6.3                     
#>  [10] zlibbioc_1.47.0              
#>  [11] vctrs_0.6.3                  
#>  [12] memoise_2.0.1                
#>  [13] cyjShiny_1.0.42              
#>  [14] RCurl_1.98-1.12              
#>  [15] ggtree_3.9.1                 
#>  [16] base64enc_0.1-3              
#>  [17] htmltools_0.5.6              
#>  [18] AnnotationHub_3.9.2          
#>  [19] curl_5.0.2                   
#>  [20] gridGraphics_0.5-1           
#>  [21] sass_0.4.7                   
#>  [22] GeneSummary_0.99.6           
#>  [23] bslib_0.5.1                  
#>  [24] htmlwidgets_1.6.2            
#>  [25] plyr_1.8.8                   
#>  [26] cachem_1.0.8                 
#>  [27] mime_0.12                    
#>  [28] lifecycle_1.0.3              
#>  [29] pkgconfig_2.0.3              
#>  [30] gson_0.1.0                   
#>  [31] Matrix_1.6-3                 
#>  [32] R6_2.5.1                     
#>  [33] fastmap_1.1.1                
#>  [34] GenomeInfoDbData_1.2.10      
#>  [35] shiny_1.7.5                  
#>  [36] aplot_0.2.1                  
#>  [37] enrichplot_1.21.3            
#>  [38] digest_0.6.33                
#>  [39] colorspace_2.1-0             
#>  [40] patchwork_1.1.3              
#>  [41] RSQLite_2.3.1                
#>  [42] MPO.db_0.99.7                
#>  [43] filelock_1.0.2               
#>  [44] fansi_1.0.4                  
#>  [45] httr_1.4.7                   
#>  [46] polyclip_1.10-4              
#>  [47] HPO.db_0.99.2                
#>  [48] compiler_4.3.1               
#>  [49] bit64_4.0.5                  
#>  [50] withr_2.5.0                  
#>  [51] graphite_1.48.0              
#>  [52] BiocParallel_1.35.4          
#>  [53] viridis_0.6.4                
#>  [54] DBI_1.1.3                    
#>  [55] dendextend_1.17.1            
#>  [56] MASS_7.3-60                  
#>  [57] rappdirs_0.3.3               
#>  [58] rjson_0.2.21                 
#>  [59] HDO.db_0.99.1                
#>  [60] tools_4.3.1                  
#>  [61] scatterpie_0.2.1             
#>  [62] ape_5.7-1                    
#>  [63] interactiveDisplayBase_1.39.0
#>  [64] rentrez_1.2.3                
#>  [65] httpuv_1.6.11                
#>  [66] glue_1.6.2                   
#>  [67] nlme_3.1-163                 
#>  [68] GOSemSim_2.27.3              
#>  [69] promises_1.2.1               
#>  [70] grid_4.3.1                   
#>  [71] shadowtext_0.1.2             
#>  [72] reshape2_1.4.4               
#>  [73] fgsea_1.27.1                 
#>  [74] generics_0.1.3               
#>  [75] gtable_0.3.4                 
#>  [76] tidyr_1.3.0                  
#>  [77] bugsigdbr_1.8.1              
#>  [78] data.table_1.14.8            
#>  [79] tidygraph_1.2.3              
#>  [80] xml2_1.3.5                   
#>  [81] utf8_1.2.3                   
#>  [82] XVector_0.41.1               
#>  [83] stringr_1.5.0                
#>  [84] ggrepel_0.9.3                
#>  [85] BiocVersion_3.18.0           
#>  [86] pillar_1.9.0                 
#>  [87] yulab.utils_0.1.0            
#>  [88] later_1.3.1                  
#>  [89] splines_4.3.1                
#>  [90] dplyr_1.1.2                  
#>  [91] tweenr_2.0.2                 
#>  [92] treeio_1.25.4                
#>  [93] lattice_0.21-8               
#>  [94] BiocFileCache_2.9.1          
#>  [95] bit_4.0.5                    
#>  [96] tidyselect_1.2.0             
#>  [97] GO.db_3.17.0                 
#>  [98] tm_0.7-11                    
#>  [99] Biostrings_2.69.2            
#> [100] reactome.db_1.86.0           
#> [101] downlit_0.4.3                
#> [102] knitr_1.44                   
#> [103] gridExtra_2.3                
#> [104] NLP_0.2-1                    
#> [105] bookdown_0.35                
#> [106] xfun_0.40                    
#> [107] graphlayouts_1.0.0           
#> [108] stringi_1.7.12               
#> [109] lazyeval_0.2.2               
#> [110] ggfun_0.1.3                  
#> [111] yaml_2.3.7                   
#> [112] evaluate_0.21                
#> [113] codetools_0.2-19             
#> [114] wordcloud_2.6                
#> [115] qvalue_2.33.0                
#> [116] tibble_3.2.1                 
#> [117] BiocManager_1.30.22          
#> [118] graph_1.79.1                 
#> [119] ggplotify_0.1.2              
#> [120] cli_3.6.1                    
#> [121] xtable_1.8-4                 
#> [122] munsell_0.5.0                
#> [123] jquerylib_0.1.4              
#> [124] Rcpp_1.0.11                  
#> [125] GenomeInfoDb_1.37.4          
#> [126] dbplyr_2.3.3                 
#> [127] png_0.1-8                    
#> [128] XML_3.99-0.14                
#> [129] parallel_4.3.1               
#> [130] ellipsis_0.3.2               
#> [131] blob_1.2.4                   
#> [132] DOSE_3.27.2                  
#> [133] bitops_1.0-7                 
#> [134] tidytree_0.4.5               
#> [135] viridisLite_0.4.2            
#> [136] slam_0.1-50                  
#> [137] scales_1.2.1                 
#> [138] purrr_1.0.2                  
#> [139] crayon_1.5.2                 
#> [140] GetoptLong_1.0.5             
#> [141] rlang_1.1.1                  
#> [142] fastmatch_1.1-4              
#> [143] cowplot_1.1.1                
#> [144] KEGGREST_1.41.0
```
