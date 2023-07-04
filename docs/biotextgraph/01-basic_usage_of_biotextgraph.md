# Basic usage

Load the package and the database for converting identifiers.
In this example, we use mostly human-derived data, and use `org.Hs.eg.db`.






```r
library(biotextgraph)
library(org.Hs.eg.db)
library(ggplot2)
library(ggraph)
library(RColorBrewer)
library(ReactomePA)
library(clusterProfiler)
library(dendextend)
library(dplyr)
load(system.file("extdata", "sysdata.rda", package = "biotextgraph"))
```

## Producing networks

The main function is producing networks between words based on the co-occurrence or correlation of the words in the text, given the list of gene identifiers. `refseq` is a function that imports gene text from RefSeq and summarizes the text information. Here, we use seven ERCC genes as input. This function returns a `biotext` class object, which contains various types of information. The net slot stores the `ggraph`, which represents the visualization result of the network.


```r
## Configure input genes
inpSymbol <- c("ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","ERCC6","ERCC8")
net <- refseq(inpSymbol)
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
net@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/basic3-1.png" width="100%" style="display: block; margin: auto;" />

By default, the gene ID type is set to `SYMBOL`. The other type can be set by `keyType`. As many of the words are commonly observed, filtering based on pre-computed word frequency on whole RefSeq data is provided. You should limit word frequency by `excludeFreq`, which is default to 2000. TF-IDF on the all the summary is also precomputed, and `exclude="tfidf"` can be specified too.


```r
net <- refseq(inpSymbol, excludeFreq=1000)
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 141 words (frequency and/or tfidf)
net
#> Type: refseq
#> Number of words: 15
#> Query: ERCC1/ERCC2/ERCC3/ERCC4/ERCC5/ERCC6/ERCC8
#> Graph: V(13), E(34)
#> Degree: excision(8)/complementation(7)/Defects(7)/nucleotide(7)/incision(6)
#> 211.5 Kb
```

For visualization, The edge label corresponding to correlation or cooccurrence values can be shown by `edgeLabel=TRUE`. The number of words to be shown on plot can be specified by `numWords`. The threshold of correlation can be specified by `corThresh`. The visualized network layout can be specified by passing `layout` argument. The text color can be changed by `colorText=TRUE`. The type of edge can be specified by `edgeLink`, which is by default `TRUE` (link will be used).


```r
net <- refseq(inpSymbol, plotType="network",
                     edgeLabel=TRUE, corThresh=0.4,
                     numWords=20, colorText=TRUE, layout="kk")
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
net@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/ngramnet-1.png" width="100%" style="display: block; margin: auto;" />
One of the main questions is which words can be clustered together among the words contained in the queried gene cluster. Word clustering (`pvclust`) and identified significant clusters based on the occurrence in the text can be visualized by specifying `tag=TRUE`. For the significance threshold, `pvclAlpha` can be specified. The default parameters perform `pvclust` on the subset of dataset for words with high frequency specified by `numWords`. If one want to perform on whole matrix (which is a natural way), `tagWhole=TRUE` can be specified, although it is computationally intensive. One can pass clusters to perform parallel computing owning to `pvclust` function, by specifying `cl` as below. For tag coloring, `tagPalette` can be used. The `igraph` contained in the object can also be plotted by passing to `plot` function.


```r
net <- refseq(inpSymbol, plotType="network", corThresh=0.2,
                     numWords=20, tag=TRUE)
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> Bootstrap (r = 0.5)... Done.
#> Bootstrap (r = 0.6)... Done.
#> Bootstrap (r = 0.7)... Done.
#> Bootstrap (r = 0.8)... Done.
#> Bootstrap (r = 0.9)... Done.
#> Bootstrap (r = 1.0)... Done.
#> Bootstrap (r = 1.1)... Done.
#> Bootstrap (r = 1.2)... Done.
#> Bootstrap (r = 1.3)... Done.
#> Bootstrap (r = 1.4)... Done.
net@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/tagnet-1.png" width="100%" style="display: block; margin: auto;" />

```r
net@pvpick
#> $clusters
#> $clusters[[1]]
#> [1] "pigmentosum" "xeroderma"  
#> 
#> $clusters[[2]]
#>  [1] "activity"             "atpdependent"        
#>  [3] "complementation"      "defects"             
#>  [5] "factor"               "functions"           
#>  [7] "incision"             "interacts"           
#>  [9] "nucleotide"           "product"             
#> [11] "transcriptioncoupled"
#> 
#> $clusters[[3]]
#> [1] "dna"    "repair"
#> 
#> 
#> $edges
#> [1]  2 12 13
plot(net)
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/tagnet-2.png" width="100%" style="display: block; margin: auto;" />

The other important aim is finding important genes (queries) in the network. The genes associated with the words can be shown by specifying `genePlot=TRUE`, useful for assessing which words is associated with interesting genes. The edges connecting words to corresponding genes are shown.
One can specify `genePlotNum` for limiting the genes shown by ranking of how often a gene is associated with the high-frequency words. This can be useful for identifying important genes among the network.


```r
net <- refseq(inpSymbol, plotType="network",
                     genePlot=TRUE, corThresh=0.5,
                     tag=TRUE, edgeLink=FALSE,
                     numWords=20)
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> Bootstrap (r = 0.5)... Done.
#> Bootstrap (r = 0.6)... Done.
#> Bootstrap (r = 0.7)... Done.
#> Bootstrap (r = 0.8)... Done.
#> Bootstrap (r = 0.9)... Done.
#> Bootstrap (r = 1.0)... Done.
#> Bootstrap (r = 1.1)... Done.
#> Bootstrap (r = 1.2)... Done.
#> Bootstrap (r = 1.3)... Done.
#> Bootstrap (r = 1.4)... Done.
net@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/gplot-1.png" width="100%" style="display: block; margin: auto;" />

The associated enriched pathways (if present) can be shown by specifying `genePathPlot`, using `ggforce`. In this option, the function first performs over-representation analysis on the whole gene set, and plot enriched terms for included genes in the plot. Enrichment analysis here is performed by the library `clusterProfiler` or `ReactomePA`, and one can control which pathways to plot by `genePathPlotSig` value.


```r
library(concaveman)
library(ggforce)
net <- refseq(inpSymbol, plotType="network",
                     genePathPlot="reactome", corThresh=0.5,
                     tag=TRUE, edgeLink=FALSE,
                     genePathPlotSig=0.05, numWords=20)
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> Bootstrap (r = 0.5)... Done.
#> Bootstrap (r = 0.6)... Done.
#> Bootstrap (r = 0.7)... Done.
#> Bootstrap (r = 0.8)... Done.
#> Bootstrap (r = 0.9)... Done.
#> Bootstrap (r = 1.0)... Done.
#> Bootstrap (r = 1.1)... Done.
#> Bootstrap (r = 1.2)... Done.
#> Bootstrap (r = 1.3)... Done.
#> Bootstrap (r = 1.4)... Done.
#> Found 47 enriched term

net@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/gpplot-1.png" width="100%" style="display: block; margin: auto;" />

By default, the associated genes are plotted without colorization (grey) and without nodes. If preferred, set `colorize=TRUE` to colorize the associated genes by `geneColor` and showing the nodes by adding pseudo-frequency corresponding to the minimum frequency of words in the network. In this way, color of nodes corresponding to words are shown by the gradient of frequency, and the queried genes are shown by `geneColor`. `addFreqToGene` also adds pseudo-frequency to gene nodes and show the node, but are colorized according to the specified minimum frequency.


```r
net <- refseq(inpSymbol, plotType="network",
                     genePlot=TRUE, corThresh=0.5,
                     colorize=TRUE, geneColor="pink",
                     colorText=TRUE,
                     tag=TRUE, edgeLink=FALSE,
                     numWords=20)
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> Bootstrap (r = 0.5)... Done.
#> Bootstrap (r = 0.6)... Done.
#> Bootstrap (r = 0.7)... Done.
#> Bootstrap (r = 0.8)... Done.
#> Bootstrap (r = 0.9)... Done.
#> Bootstrap (r = 1.0)... Done.
#> Bootstrap (r = 1.1)... Done.
#> Bootstrap (r = 1.2)... Done.
#> Bootstrap (r = 1.3)... Done.
#> Bootstrap (r = 1.4)... Done.
net@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/gplot2-1.png" width="100%" style="display: block; margin: auto;" />

## Producing word clouds

The other basic usage of the package is producing a word cloud of summaries of identifiers. Specify `plotType="wc"` or use `plot_wordcloud()` function for this purpose.


```r
wc <- obtain_refseq(c("DDX41","PNKP","IRF3")) |> make_corpus() |> make_TDM() |> plot_wordcloud()
#> Input genes: 3
#>   Converted input genes: 3
gwc <- refseq(inpSymbol, plotType="wc")
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
gwc@wc
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/basic-1.png" width="100%" style="display: block; margin: auto;" />

It accepts values of the `wordcloud()` function. `numWords` specifies how many words are to be shown on word cloud. The words are ordered by their frequency, and the subset of 1:`numWords` is used to downstream visualization. The arguments for wordcloud visualization must be passed to `argList` arguments as a list. `scaleFreq` can be specified to scale the frequency when the observation count is low.


```r
gwc <- refseq(inpSymbol,
                     plotType="wc",
                     numWords=100,
                     scaleFreq=2,
                     excludeFreq=5000,
                     argList=list(
                      "random.order"=FALSE,
                      colors=RColorBrewer::brewer.pal(8, "Dark2"),
                     "rot.per"=0.4))
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 30 words (frequency and/or tfidf)
gwc@wc
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/basic2-1.png" width="100%" style="display: block; margin: auto;" />

By default, `preserve=TRUE`, which indicates the funciton tries to preserve the original cases of characters. Note that if both the lower case words and capitalized words are present, all words are converted to capitalized words, like `damage` and `Damage` would be shown as `Damage` if `preserve=TRUE`.


```r
gwc_p <- refseq(inpSymbol,
                     plotType="wc",
                     numWords=100,
                     excludeFreq=5000,
                     preserve=FALSE,
                     argList=list(
                       rot.per=0.4,
                       colors=RColorBrewer::brewer.pal(8, "Set2"),
                       random.order=FALSE
                     ))
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 30 words (frequency and/or tfidf)
gwc_p@wc
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/basic_pres_false-1.png" width="100%" style="display: block; margin: auto;" />

It also returns a data frame consisting of frequency of each term in the slot name `freqDf`.


```r
gwc
#> Type: refseq
#> Number of words: 100
#> ERCC1/ERCC2/ERCC3/ERCC4/ERCC5/ERCC6/ERCC8
#> 483.1 Kb
knitr::kable(
  head(gwc@freqDf), caption = 'Term frequencies.',
  row.names = FALSE
)
```



Table: (\#tab:table)Term frequencies.

|word          | freq|
|:-------------|----:|
|repair        |   16|
|DNA           |   11|
|syndrome      |    9|
|excision      |    8|
|transcription |    7|
|cockayne      |    6|

N-gram is supported by library `tm`, specified by `ngram`.
Default is `1`, and the example specifying `2` is shown below.


```r
gwc2 <- refseq(inpSymbol,
                      ngram=2,
                      numWords=50,
                      plotType="wc",
                      argList=list(
                       rot.per=0.4,
                       colors=RColorBrewer::brewer.pal(8, "Set2"),
                       random.order=FALSE
                      ))
#> Input genes: 7
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> Scale for size is already present.
#> Adding another scale for size, which will replace the
#> existing scale.
gwc2@wc
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/ngramwc-1.png" width="100%" style="display: block; margin: auto;" />

Using `clusterProfiler` functions, one can use enriched pathway names for visualization.
The `enrich` option can be specified for `'kegg'` or `'reactome'`, this time we specify `'reactome'`.


```r
gwc3 <- refseq(inpSymbol, plotType="wc",
                      enrich="reactome",
                      tfidf=TRUE, numWords=50)
#> Input genes: 7
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 7
#> Performing enrichment analysis
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> Scale for size is already present.
#> Adding another scale for size, which will replace the
#> existing scale.
gwc3@wc
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/ngrampath-1.png" width="100%" style="display: block; margin: auto;" />

In the word cloud, it is also possible to visualize tag information with colors. In the example below, clustering was performed for all matrices, and the results were visualized based on the colors of `tagPalette`.


```r
## Prepare the palette for tag coloring
pal <- RColorBrewer::brewer.pal(8, "Dark2") 
pal <- colorRampPalette(pal)(20)
## Cluster on whole matrix
gwclWhole <- refseq(inpSymbol,
                     numWords=50,
                     plotType="wc",
                     tag=TRUE, tagWhole=TRUE,
                     tagPalette = pal,
                     scaleFreq=5,
                     cl=snow::makeCluster(8),
                     argList=list(rot.per=0.4))
#> Input genes: 7
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> Multiscale bootstrap... Done.
#> Scale for size is already present.
#> Adding another scale for size, which will replace the
#> existing scale.
gwclWhole@pvpick
#> $clusters
#> $clusters[[1]]
#>  [1] "abnormally"        "active"           
#>  [3] "csb"               "defective"        
#>  [5] "disease"           "hereditary"       
#>  [7] "identified"        "iih"              
#>  [9] "p44"               "radiation"        
#> [11] "repeat"            "sensitive"        
#> [13] "transcriptionally" "ultraviolet"      
#> 
#> $clusters[[2]]
#>  [1] "activates"       "adjacent"        "atpase"         
#>  [4] "atpstimulated"   "dnabinding"      "downstream"     
#>  [7] "encode"          "formation"       "frame"          
#> [10] "fusion"          "individual"      "occurs"         
#> [13] "polyadenylation" "promote"         "reading"        
#> [16] "sequence"        "shares"          "sites"          
#> 
#> $clusters[[3]]
#> [1] "class" "tfiih"
#> 
#> $clusters[[4]]
#> [1] "complementation" "defects"        
#> 
#> $clusters[[5]]
#> [1] "eme1"      "ercc1"     "structure" "xp6"      
#> 
#> $clusters[[6]]
#> [1] "orf"               "piggybackderived3"
#> [3] "splice"           
#> 
#> $clusters[[7]]
#> [1] "pigmentosum" "xeroderma"  
#> 
#> $clusters[[8]]
#> [1] "interacts" "type"     
#> 
#> $clusters[[9]]
#> [1] "incision"  "including" "light"    
#> 
#> $clusters[[10]]
#>  [1] "basic"            "bivm"            
#>  [3] "cachexia"         "cancer"          
#>  [5] "cellular"         "characterized"   
#>  [7] "cognitive"        "develop"         
#>  [9] "development"      "disability"      
#> [11] "disorder"         "exists"          
#> [13] "exposure"         "function"        
#> [15] "growth"           "hypersensitivity"
#> [17] "immunoglobulin"   "increased"       
#> [19] "motif"            "neighboring"     
#> [21] "processes"        "read"            
#> [23] "referred"         "severe"          
#> [25] "singlestrand"     "skin"            
#> [27] "specific"         "susceptibility"  
#> [29] "uvinduced"        "variable"        
#> [31] "vii"              "xp7"             
#> 
#> $clusters[[11]]
#> [1] "transcriptioncoupled" "upstream"            
#> 
#> $clusters[[12]]
#>  [1] "atpdependent"        "basal"              
#>  [3] "belongs"             "btf2tfiih"          
#>  [5] "cancerprone"         "disorders"          
#>  [7] "factor"              "helicase"           
#>  [9] "helicases"           "integral"           
#> [11] "mechanism"           "rad3xpd"            
#> [13] "subfamily"           "trichothiodystrophy"
#> 
#> $clusters[[13]]
#>  [1] "alter"                     "carcinogenesis"           
#>  [3] "catalyzes"                 "cd3e"                     
#>  [5] "cerebrooculofacioskeletal" "cisplatin"                
#>  [7] "compounds"                 "crosslinks"               
#>  [9] "electrophilic"             "epsilon"                  
#> [11] "ercc4"                     "excising"                 
#> [13] "exon"                      "formed"                   
#> [15] "forms"                     "heterodimer"              
#> [17] "heterodimeric"             "induced"                  
#> [19] "interstrand"               "lesion"                   
#> [21] "lesions"                   "molecule"                 
#> [23] "overlaps"                  "pathway"                  
#> [25] "play"                      "polymorphisms"            
#> [27] "process"                   "product"                  
#> [29] "recombinational"           "required"                 
#> [31] "result"                    "strand"                   
#> [33] "xpf"                      
#> 
#> $clusters[[14]]
#> [1] "dna"    "repair"
#> 
#> 
#> $edges
#>  [1]  13  30  97  98 101 104 106 107 110 113 117 118 123 129
gwclWhole@wc
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/tagging-1.png" width="100%" style="display: block; margin: auto;" />
In this example querying ERCC genes, the term `DNA repair` is clustered as expected.

## Visualization of PubMed information.

Using `rentrez`, one can perform the same analysis on PubMed text like the article title and abstract. The function queries for the input gene symbols (or the other queries) and visualize. For typical use cases, the genes identified by showing `genePlot`, or hub genes identified in gene network analysis can be queried. The basic parameters for searching PubMed, like max number of articles retrieved and how to sort the articles can be specified by `retMax` and `sortOrder`. Be sure to obtain [an api key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) when querying heavily, and specify in `apiKey` argument. These functions including `pubmed` and `refseq` serve as wrappers for several other functions, but a detailed explanation can be found in \@ref(tidy).
 

```r
ab <- pubmed(inpSymbol[1:3], retMax=20, apiKey=apiKey, plotType="wc")
#> Scale for size is already present.
#> Adding another scale for size, which will replace the
#> existing scale.
ab@wc
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/pubmed_wc-1.png" width="100%" style="display: block; margin: auto;" />

The returned PubMed IDs are stored in `pmids` slot.


```r
ab@pmids
#>  [1] "20301571" "27051024" "27838878" "21278243" "31077069"
#>  [6] "25867436" "16835333" "32749109" "23593158" "8053936" 
#> [11] "28803404" "28088319" "28474168" "24582975" "27051038"
#> [16] "24023723" "34284736" "33125943" "25674148" "20687508"
ab@pmids |> length()
#> [1] 20
```

As fetching the same information is not desirable and time consuming, the same object can be passed to `redo` option and re-perform the analysis like tagging, or changing the visualization options. Also, using `obtain_pubmed` function, only the text can be obtained and be processed by the downstream functions. 


```r
abtag <- pubmed(redo=ab, tag=TRUE, cl=snow::makeCluster(10), apiKey=apiKey)
#> Resuming from the previous results
#> Multiscale bootstrap... Done.
abtag2 <- pubmed(redo=abtag, tag=TRUE, genePlot=TRUE,
                 plotType="network", corThresh=0.2, pre=TRUE, apiKey=apiKey)
#> Resuming from the previous results
#> Using previous pvclust results
abtag2@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/pubmed_wc_2-1.png" width="100%" style="display: block; margin: auto;" />

If only the gene symbols are to be plotted, specify `onlyGene=TRUE`.


```r
net <- pubmed(inpSymbol[1:5], plotType="network", onlyGene=TRUE, apiKey=apiKey)
#> Subsetting to the gene symbol in orgDb
net@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/onlygene-1.png" width="100%" style="display: block; margin: auto;" />

## Comparing two or more networks

One can compare two or more networks by providing list of `biotext` objects produced by text mining various databases, like `refseq`, `pubmed`, `obtain_refseq`, etc. This can be useful for assessing the similarity and dissimilarity of the various text sources, like PubMed, RefSeq, and Reactome pathway names. Additionally, performing graph-based clustering on merged networks can potentially identify groups of related terms within the overall network.


```r
cxcls <- c()
for (i in c(1,2,3,5,6,8,9,10,11,12,13,14,16)){
    cxcls <- c(cxcls, paste0("CXCL",i))
}

net1 <- refseq(inpSymbol, plotType="network",
                      corThresh=0.5, numWords=20)
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
net2 <- refseq(cxcls, plotType="network",
                      corThresh=0.5, numWords=20)
#> Input genes: 13
#>   Converted input genes: 13
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
net3 <- pubmed(redo=ab, plotType="network",
               corThresh=0.2, numWords=20)
#> Resuming from the previous results

## Not having meaningful overlaps
compareWordNet(list(net1, net2),
               titles=c("ercc","cxcl"), colNum=3)@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/compare-1.png" width="100%" style="display: block; margin: auto;" />

```r
compareWordNet(list(net1, net3),
               titles=c("ercc","ercc-PubMed"), colNum=3)@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/compare-2.png" width="100%" style="display: block; margin: auto;" />

```r
compareWordNet(list(net1, net2, net3),
               titles=c("ercc","cxcl", "ercc-PubMed"), colNum=5)@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/compare-3.png" width="100%" style="display: block; margin: auto;" />

If `tag` information is available in both gene clusters, combined tags can be visualized by `tag=TRUE`. If using `geom_mark_hull` to show the tagging information, specify `hull=TRUE` in `compareWordNet`.


```r
keggPathways <- org.Hs.egPATH2EG
mappedKeys <- mappedkeys(keggPathways)
keggList <- as.list(keggPathways[mappedKeys])

net1 <- refseq(keggList$`04110`,
                      keyType="ENTREZID",
                      plotType="network",
                      corThresh=0.3,
                      numWords=30,
                      tag=TRUE,
                      tfidf=TRUE)
#> Input genes: 124
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> Bootstrap (r = 0.5)... Done.
#> Bootstrap (r = 0.6)... Done.
#> Bootstrap (r = 0.7)... Done.
#> Bootstrap (r = 0.8)... Done.
#> Bootstrap (r = 0.9)... Done.
#> Bootstrap (r = 1.0)... Done.
#> Bootstrap (r = 1.1)... Done.
#> Bootstrap (r = 1.2)... Done.
#> Bootstrap (r = 1.3)... Done.
#> Bootstrap (r = 1.4)... Done.

net2 <- refseq(keggList$`04210`,
                      keyType="ENTREZID",
                      plotType="network",
                      corThresh=0.3,
                      numWords=30,
                      tfidf=TRUE,
                      tag=TRUE)
#> Input genes: 87
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> Bootstrap (r = 0.5)... Done.
#> Bootstrap (r = 0.6)... Done.
#> Bootstrap (r = 0.7)... Done.
#> Bootstrap (r = 0.8)... Done.
#> Bootstrap (r = 0.9)... Done.
#> Bootstrap (r = 1.0)... Done.
#> Bootstrap (r = 1.1)... Done.
#> Bootstrap (r = 1.2)... Done.
#> Bootstrap (r = 1.3)... Done.
#> Bootstrap (r = 1.4)... Done.

compareWordNet(list(net1, net2), tag=TRUE, hull=TRUE)@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/tagcomp-1.png" width="100%" style="display: block; margin: auto;" />

An application example of using text mining for transcriptome analysis of BK polyomavirus infection by combining the functions described here \@ref(app).

## Text over represenatation analysis (experimental)

For RefSeq, overrepresentation-based filtering and prioritization of words can be performed using pre-computed background.


```r
geneList <- keggList$`00785` # Lipoic acid metabolism
pvs <- textORA(geneList)
hist(pvs)
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/ora-1.png" width="100%" style="display: block; margin: auto;" />

```r
pvs[order(pvs)] |> head()
#>            lipoic              step              acid 
#>      1.600963e-14      8.016934e-08      5.300425e-07 
#> lipoateactivating  lipoatedependent            lipoyl 
#>      1.484931e-04      1.484931e-04      1.484931e-04

geneList <- keggList$`05150` # Staphylococcus aureus infection
pvs <- textORA(geneList)
hist(pvs)
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/ora-2.png" width="100%" style="display: block; margin: auto;" />

```r
pvs[order(pvs)] |> head()
#>          drb        chain   complement        class 
#> 0.000000e+00 3.669590e-82 3.266180e-74 5.496543e-46 
#>         exon         beta 
#> 6.851722e-39 1.318324e-34
```

This thresholding can be used in RefSeq visualization in the above functions.
Filter words using ORA threshold and frequency threshold by setting `ora=TRUE`.


```r
net <- refseq(inpSymbol, plotType="network",
                     ora=TRUE, edgeLink=FALSE)
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
#> Performing ORA
#> Filtered 0 words (ORA)
net@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/oranet-1.png" width="100%" style="display: block; margin: auto;" />

After obtaining the ORA results, one can plot volcano-plot like plot for the results using `plotORA` function.


```r
library(ggrepel)
plotORA(net)
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/oraplot-1.png" width="100%" style="display: block; margin: auto;" />

## Dependency analyis using udpipe

Using [`udpipe`](https://github.com/bnosac/udpipe) package ([Straka and Straková. 2017](https://aclanthology.org/K17-3009/)), one can performe dependency analysis of texts in various databases. Set `useUdpipe` to `TRUE`, and specify downloaded model to be used in `udpipeModel`.


```r
p <- biotextgraph::refseq(c("DDX41","PNKP","ERCC2"),
                                  plotType="network",
                                  useUdpipe=TRUE,
                                  udpipeModel="~/english-ewt-ud-2.5-191206.udpipe")
#> Using udpipe mode
#> Input genes: 3
#>   Converted input genes: 3
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
p@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/udpipe-1.png" width="100%" style="display: block; margin: auto;" />

## Changing the font

Font can be specified by `fontFamily`, which affects all the labels including edges and nodes.


```r
load(system.file("extdata", "sysdata.rda", package = "biotextgraph"))
degs <- d3degUpAssetta2016
## Use alien encounter fonts (http://www.hipsthetic.com/alien-encounters-free-80s-font-family/)
sysfonts::font_add(family="alien",regular="SFAlienEncounters.ttf")
showtext::showtext_auto()
p <- biotextgraph::refseq(degs,
                          plotType="network",
                          numWords=50, genePlot=TRUE,
                          fontFamily="alien",
                          colorText=TRUE)
#> Input genes: 636
#>   Converted input genes: 555
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
p@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/font-1.png" width="100%" style="display: block; margin: auto;" />

For complex networks, changing the layout is possible by the small function `changeLayout`.
Specify which layout algorithm to choose in `igrpah`.


```r
changeLayout(p, igraph::layout.graphopt)@net
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/clay-1.png" width="100%" style="display: block; margin: auto;" />

## Biofabric layout

It is possible to use the `biofabric` layout from the extracted `igraph`. In this case, one want to add text position in x-axis using `obtainTextPosition`.


```r
ex <- refseq(c("DDX41","PNKP"),
                   genePlot=TRUE)
#> Input genes: 2
#>   Converted input genes: 2
#> Filter based on GeneSummary
#> Filtered 77 words (frequency and/or tfidf)
ex <- ex@igraph |> obtainTextPosition() ## Obtain text position
ex |> 
ggraph("fabric",
  sort.by = node_rank_fabric())+
    geom_node_range() +
    geom_edge_span(end_shape = 'circle') +
    geom_node_text(aes(x=xmin-4, label=name), size=4.5)+
    shadowtext::geom_shadowtext(aes(x=center,
                       color=nodeCat,
                       y=y+1,
                       label=name),
                   size=4.5, bg.colour="white")+
    scale_color_manual(values=c("grey20", "grey80"))+
    theme_graph()
```

<img src="01-basic_usage_of_biotextgraph_files/figure-html/biof-1.png" width="100%" style="display: block; margin: auto;" />



```r
sessionInfo()
#> R version 4.2.2 (2022-10-31 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 19044)
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
#> [1] stats4    stats     graphics  grDevices utils    
#> [6] datasets  methods   base     
#> 
#> other attached packages:
#>  [1] ggforce_0.4.1         concaveman_1.1.0     
#>  [3] dplyr_1.1.0           dendextend_1.16.0    
#>  [5] clusterProfiler_4.6.2 ReactomePA_1.42.0    
#>  [7] RColorBrewer_1.1-3    ggraph_2.1.0         
#>  [9] org.Hs.eg.db_3.16.0   AnnotationDbi_1.60.0 
#> [11] IRanges_2.32.0        S4Vectors_0.36.2     
#> [13] Biobase_2.58.0        BiocGenerics_0.44.0  
#> [15] biotextgraph_0.99.0   ggplot2_3.4.1        
#> 
#> loaded via a namespace (and not attached):
#>   [1] shadowtext_0.1.2       fastmatch_1.1-3       
#>   [3] systemfonts_1.0.4      plyr_1.8.8            
#>   [5] igraph_1.4.1           lazyeval_0.2.2        
#>   [7] splines_4.2.2          BiocParallel_1.32.5   
#>   [9] GenomeInfoDb_1.34.9    digest_0.6.31         
#>  [11] yulab.utils_0.0.6      htmltools_0.5.4       
#>  [13] bugsigdbr_1.4.3        GOSemSim_2.24.0       
#>  [15] viridis_0.6.2          GO.db_3.16.0          
#>  [17] fansi_1.0.4            GeneSummary_0.99.4    
#>  [19] magrittr_2.0.3         memoise_2.0.1         
#>  [21] tm_0.7-11              Biostrings_2.66.0     
#>  [23] graphlayouts_0.8.4     pvclust_2.2-0         
#>  [25] wordcloud_2.6          sysfonts_0.8.8        
#>  [27] enrichplot_1.18.3      colorspace_2.1-0      
#>  [29] rappdirs_0.3.3         blob_1.2.3            
#>  [31] ggrepel_0.9.3          xfun_0.37             
#>  [33] crayon_1.5.2           RCurl_1.98-1.10       
#>  [35] jsonlite_1.8.4         scatterpie_0.1.8.001  
#>  [37] graph_1.76.0           ape_5.7               
#>  [39] glue_1.6.2             polyclip_1.10-4       
#>  [41] stopwords_2.3          gtable_0.3.1          
#>  [43] zlibbioc_1.44.0        XVector_0.38.0        
#>  [45] GetoptLong_1.0.5       graphite_1.44.0       
#>  [47] rentrez_1.2.3          scales_1.2.1          
#>  [49] DOSE_3.24.2            DBI_1.1.3             
#>  [51] showtextdb_3.0         Rcpp_1.0.10           
#>  [53] viridisLite_0.4.1      xtable_1.8-4          
#>  [55] tidytree_0.4.2         gridGraphics_0.5-1    
#>  [57] reactome.db_1.82.0     bit_4.0.5             
#>  [59] htmlwidgets_1.6.1      httr_1.4.5            
#>  [61] fgsea_1.24.0           ellipsis_0.3.2        
#>  [63] pkgconfig_2.0.3        XML_3.99-0.13         
#>  [65] farver_2.1.1           sass_0.4.5            
#>  [67] utf8_1.2.3             labeling_0.4.2        
#>  [69] ggplotify_0.1.0        tidyselect_1.2.0      
#>  [71] rlang_1.0.6            reshape2_1.4.4        
#>  [73] later_1.3.0            munsell_0.5.0         
#>  [75] tools_4.2.2            cachem_1.0.7          
#>  [77] downloader_0.4         cli_3.6.0             
#>  [79] generics_0.1.3         RSQLite_2.3.0         
#>  [81] gson_0.0.9             evaluate_0.20         
#>  [83] stringr_1.5.0          fastmap_1.1.1         
#>  [85] ggdendro_0.1.23        yaml_2.3.7            
#>  [87] ggtree_3.6.2           knitr_1.42            
#>  [89] bit64_4.0.5            fs_1.6.1              
#>  [91] tidygraph_1.2.3        purrr_1.0.1           
#>  [93] KEGGREST_1.38.0        showtext_0.9-5        
#>  [95] nlme_3.1-161           mime_0.12             
#>  [97] slam_0.1-50            aplot_0.1.9           
#>  [99] xml2_1.3.3             compiler_4.2.2        
#> [101] rstudioapi_0.14        png_0.1-8             
#> [103] treeio_1.22.0          tibble_3.1.8          
#> [105] tweenr_2.0.2           bslib_0.4.2           
#> [107] stringi_1.7.12         cyjShiny_1.0.34       
#> [109] highr_0.10             lattice_0.20-45       
#> [111] Matrix_1.5-3           vctrs_0.5.2           
#> [113] pillar_1.8.1           lifecycle_1.0.3       
#> [115] jquerylib_0.1.4        GlobalOptions_0.1.2   
#> [117] data.table_1.14.8      cowplot_1.1.1         
#> [119] bitops_1.0-7           httpuv_1.6.9          
#> [121] patchwork_1.1.2        qvalue_2.30.0         
#> [123] R6_2.5.1               bookdown_0.33         
#> [125] promises_1.2.0.1       gridExtra_2.3         
#> [127] codetools_0.2-18       MASS_7.3-58.1         
#> [129] rjson_0.2.21           withr_2.5.0           
#> [131] GenomeInfoDbData_1.2.9 parallel_4.2.2        
#> [133] ISOcodes_2022.09.29    ggfun_0.0.9           
#> [135] grid_4.2.2             tidyr_1.3.0           
#> [137] HDO.db_0.99.1          rmarkdown_2.20        
#> [139] downlit_0.4.2          NLP_0.2-1             
#> [141] shiny_1.7.4            base64enc_0.1-3
```
