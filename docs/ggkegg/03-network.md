

# Network


```r
library(ggkegg)
library(tidygraph)
library(dplyr)
kne <- network("N00002")
kne
#> N00002
#> BCR-ABL fusion kinase to RAS-ERK signaling pathway
```

## Combining multiple networks



```r
kne <- network("N00385")  ## HCMV
kne2 <- network("N00366") ## HPV
one <- kne |> network_graph()
two <- kne2 |> network_graph()
graph_join(one, two) |> plot_kegg_network()
```

<img src="03-network_files/figure-html/network_combine-1.png" width="100%" style="display: block; margin: auto;" />
