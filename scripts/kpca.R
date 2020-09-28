setwd('~/git/ovrf-review/data/')

# to generate color palettes
gg2.cols <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# analyze k-mer intersection distance matrix
km <- read.csv('adeno.inter-matrix.csv', header=F, row.names=1)
km <- 1-as.matrix(km)
headers <- read.csv('adeno.headers.csv', header=T)


# look at PCA of distance matrix
pca <- prcomp(as.dist(km))
par(mar=c(5,5,1,1))
plot(pca$x, pch=16, cex=0.5, col='grey')

# annotate plot with some gene names
proteins <- c('E1A', 'E1B 19K', 'E1B 55K', 'IX', 'IVa2', 'polymerase', 
              'pTP', '52K', 'pIIIa', 'III', 'pVII', '^V$')
pal <- gg2.cols(length(proteins))
for (i in 1:12) {
  points(pca$x[grepl(proteins[i], headers$gene.name), 1:2], 
         pch=21, bg=pal[i])  
}

# try visualizing in 3D
require(rgl)
plot3d(pca$x)
for (i in 1:12) {
  spheres3d(pca$x[grepl(proteins[i], headers$gene.name), 1:3], 
            col=pal[i], radius=0.1)
}


# applying hierarchical clustering directly on the distance
# matrix results in two clusters that are much larger than others
hc <- hclust(as.dist(km), method='ward.D2')


diagnostics.hclust <- function(hc, x=c(NA)) {
  par(mfrow=c(1,2))
  
  if (any(is.na(x))) {
    x <- seq(0, max(hc$height), length.out=20)  
  }
  y <- lapply(x, function(xx) as.integer(table(cutree(hc, h=xx))))
  boxplot(y, xlab='Cutoff', ylab='Cluster sizes', log='y', xaxt='n')
  axis(side=1, at=1:length(y), labels = round(x, 1))
  # actual number of genomes
  abline(h=length(unique(headers$accession)), col='red')
  
  y <- sapply(x, function(xx) length(table(cutree(hc, h=xx))))
  plot(x, y, log='y', type='b',
       xlab='Cutoff (height)', ylab='Number of clusters')
  # average number of ORF annotations per genome
  y <- table(headers$accession)
  abline(h=mean(y), col='red')
  abline(h=quantile(y, 0.25), col='red', lty=2)
  abline(h=quantile(y, 0.75), col='red', lty=2)  
  
  par(mfrow=c(1,1))  # reset parameters
}

diagnostics.hclust(hc)

# t-stochastic neighbour embedding results in more consistent cluster
# sizes
require(Rtsne)
res <- Rtsne(km, is_distance=T, verbose=T, dims=2)
hc2 <- hclust(dist(res$Y), method='ward.D2')

diagnostics.hclust(hc2, x=seq(30, 100, length.out=20))

# look at the raw hclust results
par(mfrow=c(1,2))
plot(hc, label=F, main='Direct')
plot(hc2, label=F, main='t-SNE')

# based on diagonstic plot
clusters <- cutree(hc2, h=65)
pal <- gg2.cols(n=max(clusters))
par(mfrow=c(1,1))
plot(res$Y, type='n')
text(res$Y, label=clusters, col=pal[clusters], cex=0.8)


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


