

# Custom usage - pathview text

## KO queries

KEGG mapping plots produced by pathview can be populated with text information ([Luo and Brouwer. 2013](https://academic.oup.com/bioinformatics/article/29/14/1830/232698])).


```r
library(wcGeneSummary)
library(dplyr)
library(pathview);library(png)
```

Suppose an example for KEGG Orthology obtained from some microbiome-related experiments.
By using [MicrobiomeProfiler](https://github.com/YuLab-SMU/MicrobiomeProfiler), enriched pathway using over-representation analysis can be obtained as follows.


```r
kos <- c("K00640","K23304","K01760","K01697","K00549","K00899")
enriched <- MicrobiomeProfiler::enrichKO(kos)
clusterProfiler::cnetplot(enriched)
```

<img src="05-pathview_text_files/figure-html/micro-1.png" width="672" />
Now plot these KOs in the pathview.


```r
## As the output is PNG, this script is omitted.
pvInput <- rep(1,length(kos))
names(pvInput) <- kos
pathview(kos, pathway.id = "00270", species = "ko")
grid::rasterGrob(readPNG("ko00270.pathview.png"), interpolate=FALSE)
#> rastergrob[GRID.rastergrob.93]
```
![raw pathview plot](https://github.com/noriakis/software/blob/main/images/ko00270.pathview.png?raw=true){width=500px}

To use the `pathviewText()` one should specify which database to retrieve the text, `"refseq"` or `"abstract"`.
The default is "abstract", which fetches information from RefSeq database. In this case, KO number is not included in RefSeq, thus the abstract is chosen. However, querying in KO number itself is not useful. We conver the KO number to EC number using KEGG REST API, and fetches the description using `wcEC()` function, and query the enzyme name in PubMed. In this way, one must provide `searchTerm`, and `termMap`, which map the KO and EC. Also, `node.types="ortholog"` must be passed to `pathview`.


```r

ecko <- data.table::fread("https://rest.kegg.jp/link/ko/ec",
                            header=FALSE)
ecko$EC <- sapply(strsplit(ecko$V1, ":"),"[",2)
ecko$KO <- sapply(strsplit(ecko$V2, ":"), "[", 2)
filt <- ecko[ ecko$KO %in% kos, ]
ecnum <- filt$EC
termMap <- filt[,c("KO","EC")] |> `colnames<-`(c("query","number"))
ecQuery <- wcEC("../enzyme.dat", ecnum = ecnum, onlyTerm =TRUE)
#> Processing EC file
ecMap <- wcEC("../enzyme.dat", ecnum = ecnum, onlyDf =TRUE)
#> Processing EC file
termMap <- merge(termMap, ecMap[,c("number","desc")], by="number") |>
  `colnames<-`(c("number","query","description"))

## Main function

hor <- pathviewText(kos, keyType = "KO",
     target="abstract",
     pid = "00270",
     org = "ko", 
     searchTerms = ecQuery,
     termMap = termMap,
     node.types="ortholog")
#> Proceeding without API key

hor$concat
```

<img src="05-pathview_text_files/figure-html/pathview2-1.png" width="672" />

```r

## We can transpose the barplot by trans=TRUE

ver <- pathviewText(kos, keyType = "KO",
     target="abstract",
     pid = "00270",
     org = "ko", 
	   trans = TRUE,
     searchTerms = ecQuery,
     termMap = termMap,
     node.types="ortholog")
#> Proceeding without API key

ver$concat
```

<img src="05-pathview_text_files/figure-html/pathview2-2.png" width="672" />

## Gene queries

For gene queries, we can obtain the same plot using relatively simple query.


```r
query <- c("TP53","CDC45","CDC6")

ver <- pathviewText(query,
     keyType = "SYMBOL",
     target="abstract",
     pid = "04110",
     org = "hsa")
#> Converting to ENTREZID
#>   Converted input genes: 3
#> Proceeding without API key

ver$concat
```

<img src="05-pathview_text_files/figure-html/gene-1.png" width="672" />

The returned `osplot` object can be accessed as the name `text` in the list. Note that only the nodes (words) within the final network are to be plotted. If all the words are to be included, path the arguments to `argList` to include all.


```r
ver$text
#> Type: pubmed_abstract
#> Number of words: 20
#> TP53 OR CDC45 OR CDC6
#> 383 Kb

ver <- pathviewText(query,
     keyType = "SYMBOL",
     target="abstract",
     pid = "04110",
     org = "hsa", argList=list(numWords=Inf, corThresh=0))
#> Converting to ENTREZID
#> 'select()' returned 1:1 mapping between keys and
#> columns
#>   Converted input genes: 3
#> Info: Downloading xml files for hsa04110, 1/1 pathways..
#> Info: Downloading png files for hsa04110, 1/1 pathways..
#> Info: Working in directory C:/Users/nsato/Dropbox/build_wgcs/book
#> Info: Writing image file hsa04110.custom.cols.png
#> Proceeding without API key

ver$concat
```

<img src="05-pathview_text_files/figure-html/assess-1.png" width="672" />
