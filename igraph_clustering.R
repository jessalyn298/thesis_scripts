library(igraph)
library(rgexf)

setwd("path/to/working/directory")

gexf<-read.gexf("cc_1.gexf") #gephi format file
graph<-gexf.to.igraph(gexf)
#plot.igraph(graph, vertex.size=10.0, vertex.label=NA)
cluster<-cluster_louvain(graph, weights = NULL)

#V(graph)$color <- cluster$membership+1
#graph <- set_graph_attr(graph, "layout", layout.auto(graph))
#plot(graph, vertex.label.dist=1.5, vertex.label=NA, layout=layout.auto(graph))

for(i in seq_along(cluster)) {
  Community = induced_subgraph(graph, cluster[[i]])
  V(Community)$name <- cluster[[i]] ## To preserve original node numbers
  Edgelist = as_edgelist(Community)
  FileName = paste0("Cluster", i, ".txt")
  write.table(Edgelist, FileName, row.names=FALSE, col.names=FALSE, sep=",")
}