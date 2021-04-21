require(igraph)

setwd('~/git/ovrf-review/dataset_2020/plot_data/')

files <- Sys.glob('./*/nodelist.csv')

for (f in files) {
  nl <- read.csv(f)
  # read edge list
  el <- read.csv(gsub("nodelist", "edgelist", f))
  # filter edges by minimum count
  el <- el[el$edge.count >= 5, ]
  g <- graph_from_data_frame(el, vertices=nl)
  
}

