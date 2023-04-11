# ggkegg

`ggkegg` fetches information from KEGG and parse, analyze and visualize them using `ggplot2` and `ggraph`. Combined with the other packages investigating biological functions using KEGG.


```r
# devtools::install_github("noriakis/ggkegg")
library(ggkegg)
```

One of the main aims of `ggkegg` is manupilating KEGG information in tidy ways using `tidygraph`.


```r
library(dplyr)
library(tidygraph)
parse_kgml("hsa04110") |> ## Obtain and parse the pathway
  activate(nodes) |> ## node manipulation
  mutate(convert_hsa=convert_id("hsa"),
         convert_map=convert_id("pathway")) |> ## convert IDs for organism hsa and pathway
  ggraph(x=x, y=y)+ ## ggraph plot
  geom_node_rect(aes(filter=type=="gene",
                     fill=I(bgcolor)),
                 color="black")+
  geom_node_text(aes(label=convert_hsa),
                 size=2, family="serif")
```

<img src="01-intro_files/figure-html/ggkegg-1.png" width="1344" />
