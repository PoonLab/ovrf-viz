set.seed(5)

########################################
# Import k-mer counts and headers
########################################
setwd('/home/laura/Projects/ProtClust/data')
km <- read.csv('kmer-distance.csv', header=F, row.names=1)
km <- 1-as.matrix(km)
headers <- read.csv('proteins-header.txt', header=T)

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
# Function to generate color palettes
########################################
gg2.cols <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

########################################
# Plot clustering results
########################################
pal <- gg2.cols(n=max(clusters))
par(mfrow=c(1,1))
plot(res$Y, type='n')
text(res$Y, label=clusters, col=pal[clusters], cex=0.8)

########################################
# Save clustering results as .csv file
########################################
info <- cbind(headers, clusters)
write.csv(info, 'protein-clusters-info.csv')


