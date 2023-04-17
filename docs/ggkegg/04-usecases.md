


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
