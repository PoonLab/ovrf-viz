# analyze k-mer distance matrix
km <- read.csv('~/git/ovrf-review/data/adeno.inter-matrix.csv', header=F, row.names=1)
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
plot3d(rotated(kp))
for (i in 1:12) {
  spheres3d(rotated(kp)[grepl(proteins[i], rnames), 1:3], col=pal[i], radius=0.5)
}



require(Rtsne)
res <- Rtsne(1-km, is_distance=T, verbose=T, dims=3)

plot(res$Y, pch=16, col='grey')
for (i in 1:12) {
  points(res$Y[grepl(proteins[i], rnames), ], pch=21, bg=pal[i])
}

plot3d(res$Y, pch=16, col='grey')
for (i in 1:12) {
  #points(res$Y[grepl(proteins[i], rnames), ], pch=21, bg=pal[i])
  spheres3d(res$Y[grepl(proteins[i], rnames), 1:3], col=pal[i], radius=0.7)
}
spheres3d(res$Y[grepl('hypothetical', rnames), ], col='black', radius=0.3)
