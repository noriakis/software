


# Usecases 


## Projecting the gene regulatory networks on KEGG map

With this package, it is possible to project inferred networks such as gene regulatory networks or KO networks inferred by other software onto KEGG maps. The following is an example of projecting a subset of KO networks within a pathway inferred by CBNplot onto the reference map of the corresponding pathway using `MicrobiomeProfiler`. Of course, it is also possible to project networks created using other methods.


```r
library(dplyr)
library(igraph)
library(tidygraph)
library(CBNplot)
library(ggkegg)
library(MicrobiomeProfiler)
data(Rat_data)
ko.res <- enrichKO(Rat_data)
exp.dat <- matrix(abs(rnorm(910)), 91, 10) %>% magrittr::set_rownames(value=Rat_data) %>% magrittr::set_colnames(value=paste0('S', seq_len(ncol(.))))
returnnet <- bngeneplot(ko.res, exp=exp.dat, pathNum=1, orgDb=NULL,returnNet = TRUE)
pg <- pathway("ko00650")
joined <- combine_with_bnlearn(pg, returnnet$str, returnnet$av)
```

Plot the resulting map. In this example, the strength estimated by CBNplot is first displayed with colored edges, and then the edges of the reference graph are drawn in black on top of it.


```r
joined |> 
  activate(nodes) |>
  mutate(convertKO=convert_id("ko")) |>
  activate(edges) |>
  filter() |>
  ggraph(x=x, y=y) +
  geom_edge_link(width=2,aes(filter=!is.na(strength), color=strength),
                 start_cap=circle(4,"mm"),
                 end_cap=circle(4,"mm"))+
  
  geom_edge_link(
    width=0.1, aes(filter=is.na(strength)),
    start_cap=circle(4,"mm"),
    end_cap=circle(4,"mm"))+
  scale_edge_color_gradient(low="blue",high="red")+
  geom_node_rect(fill="white", color="black")+
  geom_node_text(aes(label=convertKO), size=2)
```

<img src="04-usecases_files/figure-html/konetplot-1.png" width="100%" style="display: block; margin: auto;" />

## Visualizing numerical attributes from DESeq2

By providing the results of the DESeq2 package, which is often used for transcriptome analysis, it is possible to reflect numerical results in the nodes of a graph. The `assign_deseq2` function can be used for this purpose. By specifying the numerical value (e.g., `log2FoldChange`) that you want to reflect in the graph as the `column` argument, you can assign the value to the nodes. If multiple genes are hit, the `numeric_combine` argument specifies how to combine multiple values (the default is `mean`).


```r
library(DESeq2)
library(org.Hs.eg.db)
load("uro.deseq.res.rda") ## Storing DESeq() result
res
#> class: DESeqDataSet 
#> dim: 29744 26 
#> metadata(1): version
#> assays(8): counts avgTxLength ... replaceCounts
#>   replaceCooks
#> rownames(29744): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
#> rowData names(27): baseMean baseVar ... maxCooks
#>   replace
#> colnames(26): SRR14509882 SRR14509883 ... SRR14509906
#>   SRR14509907
#> colData names(27): Assay.Type AvgSpotLen ...
#>   viral_infection replaceable
vinf <- results(res, contrast=c("viral_infection","BKPyV (Dunlop) MOI=1","No infection"))

g <- pathway("hsa04110") |> mutate(deseq2=assign_deseq2(vinf),
                                   padj=assign_deseq2(vinf, column="padj"),
                                   converted_name=convert_id("hsa"))

ggraph(g, layout="manual", x=x, y=y) + 
  geom_edge_link(width=0.5, arrow = arrow(length = unit(1, 'mm')), 
                 start_cap = square(1, 'cm'),
                 end_cap = square(1.5, 'cm'), aes(color=subtype))+
  geom_node_rect(aes(fill=deseq2, filter=type=="gene"), color="black")+
  ggfx::with_outer_glow(geom_node_text(aes(label=converted_name, filter=type!="group"), size=2.5), colour="white", expand=1)+
  scale_fill_gradient(low="blue",high="red", name="LFC")+
  theme_void()
```

<img src="04-usecases_files/figure-html/deseq2-1.png" width="100%" style="display: block; margin: auto;" />

```r

## Adjusted p-values
ggraph(g, layout="manual", x=x, y=y) + 
  geom_edge_link(width=0.5, arrow = arrow(length = unit(1, 'mm')), 
                 start_cap = square(1, 'cm'),
                 end_cap = square(1.5, 'cm'), aes(color=subtype))+
  geom_node_rect(aes(fill=padj, filter=type=="gene"), color="black")+
  ggfx::with_outer_glow(geom_node_text(aes(label=converted_name, filter=type!="group"), size=2.5), colour="white", expand=1)+
  scale_fill_gradient(low="blue",high="red", name="padj")+
  theme_void()
```

<img src="04-usecases_files/figure-html/deseq2-2.png" width="100%" style="display: block; margin: auto;" />


## Visualizing multiple enrichment results

You can visualize the results of multiple enrichment analyses. Similar to using the `ggkegg` function with the `enrichResult` class, there is an `append_cp` function that can be used within the `mutate` function. By providing an `enrichResult` object to this function, if the pathway being visualized is present in the result, the gene information within the pathway can be reflected in the graph.


```r

## These are RDAs storing DEGs
load("degListRPTEC.rda")
load("degURO.rda")

library(org.Hs.eg.db);
library(clusterProfiler);
input_uro <- bitr(uroUp,
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)$ENTREZID
input_rptec <- bitr(gls$day3_up_rptec,
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)$ENTREZID

ekuro <- enrichKEGG(gene = input_uro)
ekrptec <- enrichKEGG(gene = input_rptec)

g1 <- pathway("hsa04110") |> mutate(uro=append_cp(ekuro, how="all"),
                                    rptec=append_cp(ekrptec, how="all"),
                                    converted_name=convert_id("hsa"))
ggraph(g1, layout="manual", x=x, y=y) + 
  geom_edge_link(width=0.5, arrow = arrow(length = unit(1, 'mm')), 
                 start_cap = square(1, 'cm'),
                 end_cap = square(1.5, 'cm'), aes(color=subtype))+
  geom_node_rect(aes(fill=uro, xmin=xmin, xmax=x,  filter=type=="gene"))+
  geom_node_rect(aes(fill=rptec, xmin=x, xmax=xmax, filter=type=="gene"))+
  scale_fill_manual(values=c("steelblue","tomato"), name="urothelial|rptec")+
  ggfx::with_outer_glow(geom_node_text(aes(label=converted_name, filter=type!="group"), size=2), colour="white", expand=1)+
  theme_void()
```

<img src="04-usecases_files/figure-html/multcp-1.png" width="100%" style="display: block; margin: auto;" />
