# analyze k-mer distance matrix
km <- read.csv('~/git/ovrf-review/data/adeno.i2-matrix.csv', 
               header=F, row.names=1)
#rnames <- km[,1]
rnames <- row.names(km)
#km <- km[,-1]  # drop row names column - repeated entries


require(kernlab)
kp <- kpca(as.kernelMatrix(as.matrix(km)))

barplot(eig(kp)[1:20] / sum(eig(kp)))

gg2.cols <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

pal <- gg2.cols(12)[c(1,4,7,10,2,5,12,8,11,3,6,9)]

par(mar=c(5,5,1,1))
plot(rotated(kp), pch=16, cex=0.5, col='grey')

proteins <- c('E1A', 'E1B 19K', 'E1B 55K', 'IX', 'IVa2', 'polymerase', 
              'pTP', '52K', 'pIIIa', 'III', 'pVII', '^V$')
for (i in 1:12) {
  points(rotated(kp)[grepl(proteins[i], rnames), 1:2], 
         pch=21, bg=pal[i])  
}

# try visualizing in 3D
require(rgl)
plot3d(rotated(kp)[,10:12])
for (i in 1:12) {
  spheres3d(rotated(kp)[grepl(proteins[i], rnames), 1:3], col=pal[i], radius=0.5)
}



require(Rtsne)
res <- Rtsne(1-km, is_distance=T, verbose=T, dims=2)

plot(res$Y, pch=16, col='grey')
for (i in 1:12) {
  points(res$Y[grepl(proteins[i], rnames), ], pch=21, bg=pal[i])
}

res <- Rtsne(1-km, is_distance=T, verbose=T, dims=3)
plot3d(res$Y, pch=16, col='grey')
for (i in 1:12) {
  #points(res$Y[grepl(proteins[i], rnames), ], pch=21, bg=pal[i])
  spheres3d(res$Y[grepl(proteins[i], rnames), 1:3], col=pal[i], radius=0.7)
}
spheres3d(res$Y[grepl('hypothetical', rnames), ], col='black', radius=0.3)



# try agglomerative clustering
hc <- hclust(as.dist(1-as.matrix(km)), method='ward.D2')

pdf(file='temp.pdf', width=100, height=16)
plot(hc, cex=0.5)
dev.off()

# 72 accessions
acc <- sapply(rnames, function(x) strsplit(x, "\\.")[[1]][1])
hist(table(acc), col='grey', border='white', main=NA, breaks=20)  
summary(as.integer(table(acc)))  # median 32 proteins per genome


# we want to get the peak as close to 72 as possible
#clusters <- cutree(hc, h=2.5)
clusters <- cutree(hc, k=32)
hist(table(clusters), col='grey', border='white', breaks=20, main=NA)
abline(v=72, lty=2)


temp <- split(as.character(
  sapply(rnames, function(x) strsplit(x, "\\.")[[1]][2])
  ), clusters)



res <- Rtsne(1-km, is_distance=T, verbose=T, dims=2)
pal <- sample(gg2.cols(32), 32)
plot(res$Y, type='n')
for (i in 1:32) {
  text(res$Y[clusters==i, ], label=i, col=pal[i], cex=0.75)
}
