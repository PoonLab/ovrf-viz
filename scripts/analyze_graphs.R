require(igraph)

setwd('~/git/ovrf-review/dataset_2020/plot_data/')

# read edge list
el <- read.csv("adenoviridae/graph_data.csv")
el <- el[el$edge.count >= 5, ]
