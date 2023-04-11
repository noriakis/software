library(graphite)
keggps <- graphite::pathways("hsapiens","kegg")
gg <- igraph::graph_from_graphnel(pathwayGraph(keggps$`Cell cycle`))
E(gg)$edgeType
graphite::keggps$`Citrate cycle (TCA cycle)`
g
