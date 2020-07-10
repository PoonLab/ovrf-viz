setwd('~/Projects/ovrf-review/data/adenoviridae/')
set.seed(5)
# to generate color palettes
gg2.cols <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# analyze k-mer intersection distance matrix
km <- read.csv('distance_matrix_kmer_adeno.csv', header=F, row.names=1)
km <- 1-as.matrix(km)
headers <- read.csv('header_adeno_kmer.csv', header=T)

# annotate plot with some gene names
proteins <- c('E1A', 'E1B 19K', 'E1B 55K', 'IX', 'IVa2', 'polymerase', 
              'pTP', '52K', 'pIIIa', 'III', 'pVII', '^V$')
pal <- gg2.cols(length(proteins))

# t-stochastic neighbour embedding results in more consistent cluster
# sizes
require(Rtsne)
res <- Rtsne(km, is_distance=T, verbose=T, dims=2)
hc2 <- hclust(dist(res$Y), method='ward.D2')

diagnostics.hclust(hc2, x=seq(30, 100, length.out=20))

## Trying a different, optimization-based approach:

# Choose cutoff that minimizes number of multiple cluster assignments
# per genome (accession)
# At one extreme (all ORFs in one clusters), this number is maximized
# At the other extreme (every ORF is a cluster), this number is minimized
#  so we want to simultaneously maximize cluster assignment across genomes
acc <- headers$accession
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
#clusters <- cutree(hc, h=12)
pal <- gg2.cols(n=max(clusters))
par(mfrow=c(1,1))
plot(res$Y, type='n')
text(res$Y, label=clusters, col=pal[clusters], cex=0.8)
info <- cbind(headers, clusters)
write.csv(info, 'opt_kpca_adeno_cluster.csv')

foo <- lapply(split(headers$gene.name, clusters), function(x) {
  tab <- table(as.character(x))
  sort(tab, decreasing=TRUE)
})

# Drawing all genomes
plot(rect(1, 1, 10, 2))

