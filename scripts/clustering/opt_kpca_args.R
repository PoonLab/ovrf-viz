#!/usr/bin/env Rscript
sayHello <- function(){
   print('hello')
}

args<-commandArgs(TRUE)

sayHello()


set.seed(5)
print(args)
name <- strsplit(args[1], split="/")[[1]][5]
out<-paste0(args[3],'clust_',name)

########################################
# Import k-mer counts and headers
########################################
km <- read.csv(args[1], header=F, row.names=1)
km <- 1-as.matrix(km)
headers <- read.csv(args[2], header=T)

########################################
# t-SNE and Hierarchical Clustering
########################################
# t-stochastic neighbour embedding results in more consistent cluster
# sizes
require(Rtsne)
res <- Rtsne(km, is_distance=T, verbose=T, dims=2)
hc2 <- hclust(dist(res$Y), method='ward.D2')

#diagnostics.hclust(hc2, x=seq(30, 100, length.out=20))

########################################
# Optimization with cutree
########################################

# Choose cutoff that minimizes number of multiple cluster assignments per genome
# At one extreme (all ORFs in one clusters), this number is maximized
# At the other extreme (every ORF is a cluster), this number is minimized
# so we want to simultaneously maximize cluster assignment across genomes

acc <- headers$accession  # Genomes labeled by accession number
acc <- as.character(acc)
n.acc <- length(unique(acc))

obj.func <- function(h) {
  clusters <- cutree(hc2, h=h)
  # proportion of ORFs in genome with unique cluster IDs
  x <- sapply(split(clusters, acc), function(x) {
    tab <- table(x)
    sum(tab==1) / length(tab)
  })
  # mean proportion of genomes carrying a given cluster ID
  y <- sapply(split(acc, clusters), length) / n.acc
  (mean(x)-mean(y))^2
}

opt <- optimize(obj.func, c(0, 100))

clusters <- cutree(hc2, h=opt$minimum)

########################################
# Save clustering results as .csv file
########################################
info <- cbind(headers, clusters)
write.csv(info, out)
