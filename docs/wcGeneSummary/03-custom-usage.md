
# Interactive inspection of gene cluster network annotated by words




```r
library(wcGeneSummary)
library(ggplot2)
library(ggraph)
library(org.Hs.eg.db)
library(clusterProfiler)
library(RColorBrewer)
```

In this example, a Bayesian network showing the module eigengenes relationship are inferred using `boot.strength` function in `bnlearn` from the weighted gene correlation network analysis (WGCNA) results. The modules are annotated by word clouds produced by `wcGeneSummary()`, and can be exported to the format of `Cytoscape.js` or `vis.js`. In this way, module relationship can be interactively inspected with the functional implications. The other functions like `wcAbst()` can be used, however, you shold specify API keys for the function makes multiple queries.


```r
## In this example, we simulate WGCNA results.
## you can just use results from WGCNA.
## Assuming WGCNA results are stored in `mod`
mod <- wcGeneSummary::returnExample()
MEs <- mod$MEs
modColors <- mod$colors
ensg <- names(modColors)

# library(bnlearn)
library(igraph)

## Replace like boot.strength(mod$MEs, R=500, algorithm = "hc")
# dag <- model2network("[ME1][ME2|ME1]") # If using bnlearn
g <- graph_from_literal( ME1-+ME2, ME1-+ME3 )

## Convert to igraph
# g <- as.igraph(dag)

## Assign edge attributes
## Skip, if you perform boot.strength, the edge attributes can be added from the result
# el <- data.frame(as_edgelist(g))
# colnames(el) <- c("from","to")
# el <- left_join(el, bs)
# E(g)$strength <- el$strength
# E(g)$direction <- el$direction

## Node attributes
V(g)$stripName <- gsub("ME","",V(g)$name)
sizes <- table(modColors)
V(g)$size <- as.numeric(sizes[V(g)$stripName])

## Directory to save images and a script
rootDir <- "./"
netDir <- "visCyjs"
imageDir <- "images"

dir.create(paste0(rootDir, "/", netDir))
dir.create(paste0(rootDir, netDir, "/", imageDir))

images <- c()
plotType <- "bar"
numLim <- 200 # limit for gene number
for (i in V(g)$name){
    print(i)
    i <- as.numeric(gsub("ME","",i)) # strip ME

    queries <- ensg[modColors==i]
    if (length(queries)>numLim) {
        warning("Sampling random genes")
        queries <- queries[sample(1:length(queries), numLim)] ## Temporary restrict to randomly chosen genes, should be replaced to like kME values
    }
    
    ## Convert to ENTREZ
    entre <- AnnotationDbi::select(org.Hs.eg.db, keytype="ENSEMBL",
        keys = queries, columns = "ENTREZID")$ENTREZID
    
    if (plotType=="bar"){
        plt <- makeBar(entre, keyType="ENTREZID") # get barplot
    } else { ## If wordcloud
        # A <- wcGeneSummary(entre, keyType="ENTREZID",
        #                    argList=list(rot.per=0.4,
        #                                 colors=brewer.pal(10,
        #                                                   sample(row.names(RColorBrewer::brewer.pal.info), 1)),
        #                                 random.order=FALSE),
        #                    numWords=80)
        # # plt <- A@wc
        # # 
        # # ## This time use ggwordcloud()
        # plt <- ggwordcloud::ggwordcloud(A@freqDf$word, A@freqDf$freq,
        #                      shape="circle", min.freq = 1,max.words = Inf,
        #                      rot.per = 0.5, random.order = FALSE,
        #                      colors = brewer.pal(10,
        #                                          sample(row.names(RColorBrewer::brewer.pal.info), 1)))+
        #          scale_size_area(max_size = 40)
    }
    ## Save images
    ggsave(paste0(rootDir, netDir, "/", imageDir, "/", i ,".png"),
           plt, dpi=300, width=10, height=10)
    ## Store image dir
    images <- c(images, paste0(imageDir, "/", i ,".png"))
}
#> [1] "ME1"
#> Input genes: 12
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> [1] "ME2"
#> Input genes: 13
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> [1] "ME3"
#> Input genes: 7
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
V(g)$image <- images

## Node shape
if (plotType=="bar"){
    V(g)$shape <- rep("rectangle", length(V(g))) 
} else {
    V(g)$shape <- rep("circle", length(V(g)))
}

## Scale the node size
sizeMin <- 50
sizeMax <- 200
rawMin <- min(V(g)$size)
rawMax <- max(V(g)$size)
scf <- (sizeMax-sizeMin)/(rawMax-rawMin)
V(g)$size <- scf * V(g)$size + sizeMin - scf * rawMin

## Export
exportCyjs(g, rootDir, netDir)
# or, exportVisjs(g, rootDir, netDir)
```

Use like `http-server` in the directory containing a exported JavaScript, and interactively inspect the module relationship with word information. The example visualization is shown below (not by the code above, but Bayesian network of module eigenes inferred from RNA-seq dataset ofbladder cancer).

![Example visualization of a Bayesian network](https://github.com/noriakis/software/blob/main/images/wcbn.png?raw=true){width=50%}


Interactive inspection is possible using GitHub pages or the other hosting services like below.


```r
knitr::include_url("https://noriakis.github.io/cyjs_test/wordcloud")
```

<iframe src="https://noriakis.github.io/cyjs_test/wordcloud" width="100%" height="400px" data-external="1" style="border: none;"></iframe>

If you specify node attribute named `group`, and set `bubble=TRUE` in `exportCyjs` function, bubble sets will be plotted using [`cytoscape.js-bubblesets`](https://github.com/upsetjs/cytoscape.js-bubblesets), useful for inspecting the similarity betwteen the gene cluster, like the output of `pvclust` and `pvpick` on module eigengenes.


```r
V(g)$group <- c(1, 1, NA)
# exportCyjs(g, rootDir, netDir, bubble=TRUE)

## Example, not by the code above.
knitr::include_url("https://noriakis.github.io/cyjs_test/wordcloud_bubble")
```

<iframe src="https://noriakis.github.io/cyjs_test/wordcloud_bubble" width="100%" height="400px" data-external="1" style="border: none;"></iframe>


`vis.js` can be used, by exporting function `exportVisjs`.
In this example, the barplot of words are shown in the nodes.


```r
knitr::include_url("https://noriakis.github.io/cyjs_test/visjs")
```

<iframe src="https://noriakis.github.io/cyjs_test/visjs" width="100%" height="400px" data-external="1" style="border: none;"></iframe>


# Annotating gene cluster dendrogram

The relationship between gene clusters are often investigated in clustering analysis like WGCNA. As workflows involving gene clustering analysis typically plot dendrogram and heatmap of module eigengenes using `plotEigengeneNetworks`, it is useful to combine with wcGeneSummary, which plot additional word information on a dendrogram with one line.


```r
# WGCNA::plotEigengeneNetworks(mod$MEs, mod$colors, plotHeatmaps = FALSE)
plotEigengeneNetworksWithWords(MEs, modColors)
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
#> Input genes: 12
#>   Converted input genes: 12
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 13
#>   Converted input genes: 13
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
```

<img src="03-custom-usage_files/figure-html/wgcna-1.png" width="100%" style="display: block; margin: auto;" />

This calculates a dendrogram using `pvclust` internally in default. If you would like to plot segments involving only the specified gene cluster, use `candidateNodes` to specify the nodes.


```r
plotEigengeneNetworksWithWords(MEs, modColors, candidateNodes=c("ME2"))
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
#> Input genes: 12
#>   Converted input genes: 12
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 13
#>   Converted input genes: 13
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
```

<img src="03-custom-usage_files/figure-html/wgcna2-1.png" width="100%" style="display: block; margin: auto;" />


By default, the function calculates the frequency of common words across branches and plot the words that the differences between branches are large. By disabling `takeIntersect`, the function plots the frequent words for each branch.


```r
plotEigengeneNetworksWithWords(MEs,
    modColors, takeIntersect=FALSE, candidateNodes=c("ME2"))
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
#> Input genes: 12
#>   Converted input genes: 12
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 13
#>   Converted input genes: 13
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
```

<img src="03-custom-usage_files/figure-html/wgcnaInt-1.png" width="100%" style="display: block; margin: auto;" />


To plot the word cloud instead of pyramid plots, use `useWC` option. For scaling the word size, use `wcScale` option.


```r
scale4 <- plotEigengeneNetworksWithWords(MEs, modColors, useWC=TRUE, candidateNodes=c("ME2"), wcScale=4)
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
#> Input genes: 25
#>   Converted input genes: 25
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
scale15 <- plotEigengeneNetworksWithWords(MEs, modColors, useWC=TRUE, candidateNodes=c("ME2"), wcScale=15)
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
#> Input genes: 25
#>   Converted input genes: 25
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
scale4 + scale15
```

<img src="03-custom-usage_files/figure-html/wgcnawc-1.png" width="100%" style="display: block; margin: auto;" />

This uses `ggwordcloud` and a list specified by `wcArg` is passed to the function.


```r
plotEigengeneNetworksWithWords(MEs, modColors, useWC=TRUE, candidateNodes=c("ME2"), wcScale=15, wcArg=list(rot.per=0))
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
#> Input genes: 25
#>   Converted input genes: 25
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
```

<img src="03-custom-usage_files/figure-html/wgcnawc3-1.png" width="100%" style="display: block; margin: auto;" />


The horizontal plot can be specified by `horiz=TRUE`.


```r
plotEigengeneNetworksWithWords(MEs,
                               modColors,
                               useWC=TRUE,
                               candidateNodes=c("ME2"),
                               wcScale=15,
                               wcArg=list(rot.per=0),
                               horiz=TRUE)
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
#> Input genes: 25
#>   Converted input genes: 25
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
```

<img src="03-custom-usage_files/figure-html/wgcnawchor-1.png" width="100%" style="display: block; margin: auto;" />

For examining enriched pathway names in the dendrograms, specify `argList` to `wcGeneSummary`, like `list(enrich="kegg")`.


```r
plotEigengeneNetworksWithWords(MEs, modColors, type="words", argList=list(enrich="kegg"))
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
#> Input genes: 12
#>   Converted input genes: 12
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Performing enrichment analysis
#> Input genes: 13
#>   Converted input genes: 13
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Performing enrichment analysis
```

<img src="03-custom-usage_files/figure-html/wgcnapath-1.png" width="100%" style="display: block; margin: auto;" />

Other than textual information, we can simply annotate the dendrogram using enrichment analysis.
Useful for inspecting how the branches of dendrogram contains pathway information. 


```r
plotEigengeneNetworksWithWords(MEs, modColors, type="enrich")
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
```

<img src="03-custom-usage_files/figure-html/wgcna4-1.png" width="100%" style="display: block; margin: auto;" />

The column names for clusterProfiler results can be specified to `showType`.


```r
plotEigengeneNetworksWithWords(MEs, modColors, type="enrich", showType="Description")
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
```

<img src="03-custom-usage_files/figure-html/wgcnaDesc-1.png" width="100%" style="display: block; margin: auto;" />

Text sizing and wrapping can be controlled by `textSize` and `wrap`.


```r
plotEigengeneNetworksWithWords(MEs, modColors, type="enrich", showType="Description",
    textSize=1.5, wrap=30)
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
```

<img src="03-custom-usage_files/figure-html/wgcnaDesc2-1.png" width="100%" style="display: block; margin: auto;" />

If you have a specifically interested pathway, use `highlight` to highlight the names in the dendrogram.


```r
plotEigengeneNetworksWithWords(mod$MEs, mod$colors,
                               type="enrich", highlight=c("hsa04060"))
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
#> 'select()' returned 1:1 mapping between keys and
#> columns
#> 'select()' returned 1:1 mapping between keys and
#> columns
```

<img src="03-custom-usage_files/figure-html/highlight-1.png" width="100%" style="display: block; margin: auto;" />


`spacer` can control the gaps above and below the grob on the dendrogram (y-axis).
`horizontalSpacer` can be used too for x-axis.


```r
plotEigengeneNetworksWithWords(mod$MEs, mod$colors,
                               type="enrich", highlight=c("hsa04060"), spacer=0.2)
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
#> 'select()' returned 1:1 mapping between keys and
#> columns
#> 'select()' returned 1:1 mapping between keys and
#> columns
```

<img src="03-custom-usage_files/figure-html/highlight_2-1.png" width="100%" style="display: block; margin: auto;" />

Using a popular tree visualization library `ggtree` and `ggimage`, plotting of word cloud on the tip nodes is possible by specifying `tipWC=TRUE` . This saves the word cloud plots in the `imageDir`, and show these images on the tip nodes. Various positional parameters as well as aspect ratio can be passed to the function to correctly plot the images. Note that it cannot show all the tip nodes, so interesting nodes should be specified in `tipWCNodes`. In this mode, one must pass object that can be passed to `ggtree()`.


```r
hcl <- hclust(dist(t(mod$MEs)))
gr <- plotEigengeneNetworksWithWords(mod$MEs,
                                     mod$colors,
                                     dendPlot="ggtree",
                                     dhc=hcl,
                                     useWC=TRUE,
                                     tipWC=TRUE,
                                     tipWCNodes=c("ME3"),
                                     imageDir="testimage",
                                     candidateNodes=c("ME2"),
                                     wcScale=10,
                                     tipSize = .3,
                                     offset = .1,
                                     asp=1.4)
#> Annotating tip by word cloud
#>   Producing ME3
#> Input genes: 7
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Scale for size is already present.
#> Adding another scale for size, which will replace the
#> existing scale.
#> Warning in brewer.pal(10, sample(row.names(RColorBrewer::brewer.pal.info), : n too large, allowed maximum for palette Pastel1 is 9
#> Returning the palette you asked for with that many colors
#> Scale for size is already present.
#> Adding another scale for size, which will replace the
#> existing scale.
#> Input genes: 25
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 25
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Scale for size is already present.
#> Adding another scale for size, which will replace the
#> existing scale.
#> Scale for size is already present.
#> Adding another scale for size, which will replace the
#> existing scale.
gr
#> Some words could not fit on page. They have been removed.
```

<img src="03-custom-usage_files/figure-html/ggtree-1.png" width="100%" style="display: block; margin: auto;" />

By specifying `returnGlobOnly`, the grobs with the position in the dendrogram can be returned.



```r
gro <- plotEigengeneNetworksWithWords(mod$MEs, mod$colors, candidateNodes=c("ME2"),
                               returnGlobOnly=TRUE)
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
#> Input genes: 12
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 12
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 13
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 13
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
gro[[1]]$plot
#> TableGrob (20 x 47) "layout": 52 grobs
#>     z         cells       name
#> 1   5 ( 9- 9, 7- 7)   spacer-1
#> 2   7 (11-11, 7- 7)   axis-l-1
#> 3   3 (13-13, 7- 7)   spacer-1
#> 4   6 ( 9- 9, 9- 9)   axis-t-1
#> 5   1 (11-11, 9- 9)    panel-1
#> 6   9 (13-13, 9- 9)   axis-b-1
#> 7   4 ( 9- 9,11-11)   spacer-1
#> 8   8 (11-11,11-11)   axis-r-1
#> 9   2 (13-13,11-11)   spacer-1
#> 10 10 ( 8- 8, 9- 9)   xlab-t-1
#> 11 11 (14-14, 9- 9)   xlab-b-1
#> 12 12 (11-11, 6- 6)   ylab-l-1
#> 13 13 (11-11,12-12)   ylab-r-1
#> 14 14 ( 5- 5, 9- 9) subtitle-1
#> 15 15 ( 4- 4, 9- 9)    title-1
#> 16 16 (15-15, 9- 9)  caption-1
#> 17 17 ( 3- 3, 3- 3)      tag-1
#> 18 22 ( 9- 9,22-22)   spacer-2
#> 19 24 (11-11,22-22)   axis-l-2
#> 20 20 (13-13,22-22)   spacer-2
#> 21 23 ( 9- 9,24-24)   axis-t-2
#> 22 18 (11-11,24-24)    panel-2
#> 23 26 (13-13,24-24)   axis-b-2
#> 24 21 ( 9- 9,26-26)   spacer-2
#> 25 25 (11-11,26-26)   axis-r-2
#> 26 19 (13-13,26-26)   spacer-2
#> 27 27 ( 8- 8,24-24)   xlab-t-2
#> 28 28 (14-14,24-24)   xlab-b-2
#> 29 29 (11-11,21-21)   ylab-l-2
#> 30 30 (11-11,27-27)   ylab-r-2
#> 31 31 ( 5- 5,24-24) subtitle-2
#> 32 32 ( 4- 4,24-24)    title-2
#> 33 33 (15-15,24-24)  caption-2
#> 34 34 ( 3- 3,18-18)      tag-2
#> 35 39 ( 9- 9,37-37)   spacer-3
#> 36 41 (11-11,37-37)   axis-l-3
#> 37 37 (13-13,37-37)   spacer-3
#> 38 40 ( 9- 9,39-39)   axis-t-3
#> 39 35 (11-11,39-39)    panel-3
#> 40 43 (13-13,39-39)   axis-b-3
#> 41 38 ( 9- 9,41-41)   spacer-3
#> 42 42 (11-11,41-41)   axis-r-3
#> 43 36 (13-13,41-41)   spacer-3
#> 44 44 ( 8- 8,39-39)   xlab-t-3
#> 45 45 (14-14,39-39)   xlab-b-3
#> 46 46 (11-11,36-36)   ylab-l-3
#> 47 47 (11-11,42-42)   ylab-r-3
#> 48 48 ( 5- 5,39-39) subtitle-3
#> 49 49 ( 4- 4,39-39)    title-3
#> 50 50 (15-15,39-39)  caption-3
#> 51 51 ( 3- 3,33-33)      tag-3
#> 52  0 ( 1-20, 1-47) background
#>                                           grob
#> 1                               zeroGrob[NULL]
#> 2          absoluteGrob[GRID.absoluteGrob.758]
#> 3                               zeroGrob[NULL]
#> 4                               zeroGrob[NULL]
#> 5                     gTree[panel-1.gTree.756]
#> 6          absoluteGrob[GRID.absoluteGrob.757]
#> 7                               zeroGrob[NULL]
#> 8                               zeroGrob[NULL]
#> 9                               zeroGrob[NULL]
#> 10                              zeroGrob[NULL]
#> 11 zeroGrob[axis.title.x.bottom..zeroGrob.759]
#> 12   zeroGrob[axis.title.y.left..zeroGrob.760]
#> 13                              zeroGrob[NULL]
#> 14       zeroGrob[plot.subtitle..zeroGrob.762]
#> 15          zeroGrob[plot.title..zeroGrob.761]
#> 16        zeroGrob[plot.caption..zeroGrob.764]
#> 17            zeroGrob[plot.tag..zeroGrob.763]
#> 18                              zeroGrob[NULL]
#> 19         absoluteGrob[GRID.absoluteGrob.776]
#> 20                              zeroGrob[NULL]
#> 21                              zeroGrob[NULL]
#> 22                    gTree[panel-1.gTree.774]
#> 23         absoluteGrob[GRID.absoluteGrob.775]
#> 24                              zeroGrob[NULL]
#> 25                              zeroGrob[NULL]
#> 26                              zeroGrob[NULL]
#> 27                              zeroGrob[NULL]
#> 28 zeroGrob[axis.title.x.bottom..zeroGrob.777]
#> 29   zeroGrob[axis.title.y.left..zeroGrob.778]
#> 30                              zeroGrob[NULL]
#> 31       zeroGrob[plot.subtitle..zeroGrob.780]
#> 32          zeroGrob[plot.title..zeroGrob.779]
#> 33        zeroGrob[plot.caption..zeroGrob.782]
#> 34            zeroGrob[plot.tag..zeroGrob.781]
#> 35                              zeroGrob[NULL]
#> 36         absoluteGrob[GRID.absoluteGrob.795]
#> 37                              zeroGrob[NULL]
#> 38                              zeroGrob[NULL]
#> 39                    gTree[panel-1.gTree.793]
#> 40         absoluteGrob[GRID.absoluteGrob.794]
#> 41                              zeroGrob[NULL]
#> 42                              zeroGrob[NULL]
#> 43                              zeroGrob[NULL]
#> 44                              zeroGrob[NULL]
#> 45 zeroGrob[axis.title.x.bottom..zeroGrob.796]
#> 46   zeroGrob[axis.title.y.left..zeroGrob.797]
#> 47                              zeroGrob[NULL]
#> 48       zeroGrob[plot.subtitle..zeroGrob.799]
#> 49          zeroGrob[plot.title..zeroGrob.798]
#> 50        zeroGrob[plot.caption..zeroGrob.801]
#> 51            zeroGrob[plot.tag..zeroGrob.800]
#> 52             rect[plot.background..rect.816]
```

## Decorating wordclouds

If needed, wordclouds can be filtered by `ggfx` or using `shadowtext`. In this case, border is set to `FALSE` and the background will be transparent for resulting grobs. If `shadowtext` is needed, specify `bg.colour` argument. If `ggfx` is needed, specify the filter function in `useggfx` and parameters in `ggfxParams`.


```r
library(ggfx)
#> Warning: package 'ggfx' was built under R version 4.2.3
plotEigengeneNetworksWithWords(MEs, modColors, useWC=TRUE, candidateNodes=c("ME2"), wcScale=4,
    bg.colour="grey80")
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
#> Input genes: 25
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 25
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Scale for size is already present.
#> Adding another scale for size, which will replace the
#> existing scale.
#> Scale for size is already present.
#> Adding another scale for size, which will replace the
#> existing scale.
```

<img src="03-custom-usage_files/figure-html/decoword-1.png" width="100%" style="display: block; margin: auto;" />

```r
plotEigengeneNetworksWithWords(MEs, modColors, useWC=TRUE, candidateNodes=c("ME2"), wcScale=4,
    useggfx="with_outer_glow", ggfxParams=list(colour="white",expand=5))
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
#> Input genes: 25
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 25
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Scale for size is already present.
#> Adding another scale for size, which will replace the
#> existing scale.
#> Scale for size is already present.
#> Adding another scale for size, which will replace the
#> existing scale.
```

<img src="03-custom-usage_files/figure-html/decoword-2.png" width="100%" style="display: block; margin: auto;" />


# Assess the occurrence of the speicific words across gene clusters


```r
library(limma)
library(ggrepel)
query <- "DNA repair"
tab <- getGeneKEGGLinks(species="hsa")
listOfGenes <- list()
for (path in unique(tab$PathwayID)){
    listOfGenes[[path]] <- subset(tab, PathwayID==path)$GeneID
}
## Random subset! The results would be different.
frq <- findTerm(query, listOfGenes[sample(length(listOfGenes), 20)],
                split=TRUE, calc="mean",
                keyType="ENTREZID")
#> Finding query in 20 clusters ...
#> Input genes: 30
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 36
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 4
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 103
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 36
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 75
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 86
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 113
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 49
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 364
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 106
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 155
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 50
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 166
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 137
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 202
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 121
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 4
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 43
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 75
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
plt <- data.frame(t(data.frame(frq, check.names=FALSE)),
                  check.names=FALSE)

plt$name <- gsub("path:", "", rownames(plt))
p <- ggplot(plt, aes(dna, repair, label = plt[,3])) +
    geom_point(color = "red")+ 
    geom_text_repel(bg.color="white")+theme_minimal()+
    xlab("dna")+ylab("repair")
p
#> Warning: ggrepel: 6 unlabeled data points (too many
#> overlaps). Consider increasing max.overlaps
```

<img src="03-custom-usage_files/figure-html/findterm-1.png" width="100%" style="display: block; margin: auto;" />

For clustering analysis like `WGCNA`, making the list and query.


```r
query <- "antiviral response"
load("./blockwiseModule.rda")
mecolors <- bwmod$color
inputList <- names(mecolors)
names(inputList) <- paste0("ME",bwmod$color)

listOfGenes <- split(inputList, names(inputList))

frq <- findTerm(query, listOfGenes,
                split=TRUE,calc="highest",
                keyType="ENSEMBL")
#> Finding query in 17 clusters ...
#> Input genes: 12943
#> 'select()' returned 1:many mapping between keys and
#> columns
#>   Converted input genes: 9943
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 5534
#> 'select()' returned 1:many mapping between keys and
#> columns
#>   Converted input genes: 4649
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 43
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 42
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 43
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 40
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 42
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 40
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 41
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 40
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 40
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 34
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 39
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 37
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 38
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 38
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 2192
#> 'select()' returned 1:many mapping between keys and
#> columns
#>   Converted input genes: 2019
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 740
#> 'select()' returned 1:many mapping between keys and
#> columns
#>   Converted input genes: 578
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 653
#> 'select()' returned 1:many mapping between keys and
#> columns
#>   Converted input genes: 536
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 621
#> 'select()' returned 1:many mapping between keys and
#> columns
#>   Converted input genes: 534
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 591
#> 'select()' returned 1:many mapping between keys and
#> columns
#>   Converted input genes: 475
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 139
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 55
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 112
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 65
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Input genes: 89
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 87
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
plt <- data.frame(t(data.frame(frq, check.names=FALSE)),
                  check.names=FALSE)
plt$name <- row.names(plt)

p <- ggplot(plt, aes(antiviral, response, label = plt[,3])) +
  geom_point(color = "blue")+ 
  geom_text_repel(bg.color="white")+theme_minimal()+
  xlab("antiviral")+ylab("response")
p
```

<img src="03-custom-usage_files/figure-html/findtermWGCNA-1.png" width="100%" style="display: block; margin: auto;" />

# Recluster the cluster using word information


```r
simExample <- returnSim(returnExample()$color,
                        keyType="ENSEMBL", argList=list(ora=TRUE))
#> Number of clusters: 3
#> 1
#> Input genes: 12
#>   Converted input genes: 12
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Performing ORA
#> Filtered 109 words (ORA)
#> 2
#> Input genes: 13
#>   Converted input genes: 13
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Performing ORA
#> Filtered 238 words (ORA)
#> 3
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> Performing ORA
#> Filtered 148 words (ORA)
heatmap(simExample)
```

<img src="03-custom-usage_files/figure-html/textclus-1.png" width="100%" style="display: block; margin: auto;" />


```r
simExample <- returnSim(returnExample()$color,
                        keyType="ENSEMBL",
                        argList=list(tfidf=TRUE, takeMax=TRUE))
#> Number of clusters: 3
#> 1
#> Input genes: 12
#>   Converted input genes: 12
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> 2
#> Input genes: 13
#>   Converted input genes: 13
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> 3
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
heatmap(simExample)
```

<img src="03-custom-usage_files/figure-html/textclus2-1.png" width="100%" style="display: block; margin: auto;" />


```r
simExample <- returnSim(returnExample()$color,
                        keyType="ENSEMBL",
                        argList=list(tfidf=FALSE,
                            normalize=TRUE,
                            takeMean=TRUE))
#> Number of clusters: 3
#> 1
#> Input genes: 12
#>   Converted input genes: 12
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> 2
#> Input genes: 13
#>   Converted input genes: 13
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
#> 3
#> Input genes: 7
#>   Converted input genes: 7
#> Filter based on GeneSummary
#> Filtered 65 words (frequency and/or tfidf)
heatmap(simExample)
```

<img src="03-custom-usage_files/figure-html/textclus3-1.png" width="100%" style="display: block; margin: auto;" />
