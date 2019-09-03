setwd('~/git/ovrf-review/data')
virus <- read.csv('viruses.csv')

virus$Genome.length <- as.integer(gsub(' nt$', '', virus$Genome.length))
summary(virus)

virus$Number.of.proteins[virus$Number.of.proteins=='-'] <- NA
virus$Number.of.proteins <- as.integer(as.character(virus$Number.of.proteins))
summary(virus)

table(virus$Family)
boxplot(split(log10(virus$Genome.length), virus$Family))
boxplot(split(log10(virus$Number.of.proteins), virus$Family))


# listing of reading frames in all genomes
orfs <- read.csv('orfs-fixed.csv')


overlaps <- read.csv('ovrfs-reduced.csv')


# 1. how many overlaps per genome?
levels(overlaps$accn) <- unique(orfs$accno)
temp <- sapply(split(overlaps$overlap, overlaps$accn), length)
noverlaps <- data.frame(accn=names(temp), count=temp)

# sum total number of overlaps per genome
index <- sapply(noverlaps$accn, function(accn) {
  which(grepl(as.character(accn), virus$Accession))
})

noverlaps$genome <- virus$Genome[index]

temp <- sapply(split(noverlaps$count, noverlaps$genome), sum)
index <- match(virus$Genome, names(temp))
virus$n.overlaps <- temp[index]


set.seed(5)
pal <- rainbow(n=20, v=0.8, alpha=0.5)
pal <- sample(pal, size=20, replace=F)

x <- virus$Genome.length
y <- virus$n.overlaps
y[y==0] <- 0.5
y <- jitter(y)


par(mar=c(5,5,5,1), xpd=F)

plot(x, y, log='xy', type='n', 
     xlab='Genome length (nt)', ylab='Number of overlaps',
     xaxt='n', yaxt='n')
axis(side=1, at=10^(1:7), label=10^(1:7))
axis(side=2, at=c(0.5, 5, 50, 500), label=c(0, 5, 50, 500), las=2)

taxa <- c(
  'Siphoviridae', 'Myoviridae', 'Geminiviridae', 'Podoviridae', 'Circoviridae',
  'Potyviridae', 'Papillomaviridae', 'Rhabdoviridae', 'Parvoviridae', 'Picornaviridae',
  'Flaviviridae', 'Tolecusatellitidae', 'Polyomaviridae', 'Betaflexiviridae', 'Genomoviridae',
  'Anelloiviridae', 'Herpesviridae', 'Alphasatellitidae', 'Retroviridae', 'Herellelviridae'
)
points(x[!is.element(virus$Family, taxa)], y[!is.element(virus$Family, taxa)], 
       col=rgb(0.5,0.5,0.5,0.3), cex=0.8, pch=16)
for (i in 1:15) {
  points(x[virus$Family==taxa[i]], y[virus$Family==taxa[i]], 
         bg=pal[i], col='black', pch=(21:25)[i%%5])
}

par(xpd=NA)
legend(x=3e4, y=5000, xjust=0.5, legend=taxa, pch=rep(21:25, 4), pt.bg=pal, 
       bty='n', cex=0.7, ncol=3, x.intersp=0.5)
par(xpd=FALSE)

