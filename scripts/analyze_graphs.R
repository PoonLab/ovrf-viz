require(igraph)

setwd('~/git/ovrf-review/dataset_2020/plot_data/')

files <- Sys.glob('./*/nodelist.csv')
require(MASS)

#f <- files[1]
for (f in files) {
  nl <- read.csv(f)
  # read edge list
  el <- read.csv(gsub("nodelist", "edgelist", f))
  ael <- el[el$edge.type=='adjacent', ]
  oel <- el[el$edge.type=='overlap', ]
  
  # add overlap counts to adjacent counts
  ael$overlap <- apply(ael, 1, function(r) {
    ocount <- oel$edge.count[oel$parent == r[1] & oel$child == r[2]]
    ifelse (length(ocount) > 0, ocount, 0)
  })
  ael$edge.count <- ael$edge.count + ael$overlap
  
  # filter adjacent edges by minimum count
  ael <- ael[ael$edge.count >= 5, ]
  
  ael$total.size <- nl$size[ael$parent] + nl$size[ael$child]
  
  ag <- graph_from_data_frame(ael, vertices=nl)
  n.nodes <- nrow(nl)
  n.adj.edges <- nrow(ael)
  
  deg <- degree(ag)
  ael$edge.degree <- deg[ael$parent] + deg[ael$child]
  
  tri <- count_triangles(g)
  ael$triangles <- tri[ael$parent] + tri[ael$child]
  
  trans <- transitivity(g, 'local')
  ael$transitivity <- trans[ael$parent] + trans[ael$child]
  
  eig <- eigen_centrality(g)$vector
  ael$centrality <- eig[ael$parent] + eig[ael$child]
  
  # stepwise AIC model selection
  fit <- glm(cbind(overlap, edge.count) ~ 
               edge.degree + total.size + triangles + transitivity + centrality, 
             data=ael, family='binomial')
  fit0 <- glm(cbind(overlap, edge.count) ~ 1, data=ael, family='binomial')
  
  step.fit <- stepAIC(fit, scope=list(lower=fit0, upper=fit), trace=F)
  #density <- n.edges / choose(n.nodes, 2)
  #triad <- transitivity(g, type='global')
  #cycles <- lapply(1:n.nodes, function(i) all_simple_paths(g, i, i))
  print(summary(step.fit))
}

