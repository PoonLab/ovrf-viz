setwd('~/git/ovrf-review/data')
virus <- read.csv('viruses.csv')

# TODO: exclude the following accession numbers:
# NC_042059


virus$Genome.length <- as.integer(gsub(' nt$', '', virus$Genome.length))
#summary(virus)

virus$Number.of.proteins[virus$Number.of.proteins=='-'] <- NA
virus$Number.of.proteins <- as.integer(as.character(
  virus$Number.of.proteins))
summary(virus)


# the number of ORFs increases linearly with genome size
par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(virus$Genome.length, virus$Number.of.proteins, 
     xlim=c(0, 5e5), ylim=c(0, 1000))
plot(virus$Genome.length, virus$Number.of.proteins, log='x')



# listing of reading frames in all genomes
orfs <- read.csv('orfs-fixed.csv')


# multiple overlaps in genes with "introns" are combined into one entry
overlaps <- read.csv('ovrfs-reduced.csv')



# 1. how many overlaps per genome?
# levels(overlaps$accn) <- unique(orfs$accno)  # FIXME: this doesn't work

# add levels without overlaps
overlaps$accn <- factor(
  overlaps$accn, levels=c(
    levels(overlaps$accn), 
    setdiff(levels(orfs$accno), levels(overlaps$accn))
    ))

temp <- sapply(split(overlaps$overlap, overlaps$accn), length)
noverlaps <- data.frame(accn=names(temp), count=temp)

noverlaps$mean.olen <- sapply(split(overlaps$overlap, overlaps$accn), mean)


# sum total number of overlaps per genome
virus$first.acc <- sapply(virus$Accession, function(x) {
  xx <- as.character(x)
  if (grepl("^\\[", xx)) {
    gsub("^\\['([A-Z0-9_]+)'.+", "\\1", xx)
  } else {
    xx
  }
})

index <- match(virus$first.acc, noverlaps$accn)
virus$n.overlaps <- noverlaps$count[index]

# carry over Genome to <noverlaps>
all.acc <- sapply(1:nrow(virus), function(i) {
  xx <- as.character(virus$Accession[i])
  if (grepl("^\\[", xx)) {
    x2 <- gsub("[\\]'\\[]", "", xx, perl=T)
    strsplit(x2, ", ")[[1]]
  } else {
    xx
  }
})
all.acc <- data.frame(acc=unlist(all.acc), genome=rep(virus$Genome, times=sapply(all.acc, length)))

index2 <- match(noverlaps$accn, all.acc$acc)
noverlaps$genome <- all.acc$genome[index2]

temp <- sapply(split(noverlaps$mean.olen/noverlaps$count, noverlaps$genome), mean)
index3 <- match(virus$Genome, names(temp))
virus$len.overlaps <- ifelse(virus$n.overlaps==0, NA, temp[index3] * virus$n.overlaps)


set.seed(6)
pal <- rainbow(n=20, v=0.7)
pal <- sample(pal, size=20, replace=F)

x <- virus$Number.of.proteins
y <- virus$n.overlaps
plot(x,y)
y[y==0] <- 0.5
y <- jitter(y)


pdf(file='viruses.pdf', width=6, height=6)
par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(x, y, log='xy', 
     cex=ifelse(is.na(virus$len.overlaps), 0.1, sqrt(virus$len.overlaps)/10),
     col=ifelse(is.element(virus$Family, 
                           c('Siphoviridae', 'Myoviridae', 'Podaviridae', 'Herelleviridae')), 
                'salmon', 'black'),
     xlab='Number of ORFs', ylab='Number of overlaps')
z <- which(virus$Genome=='Hepatitis B virus')
points(x[z], y[z], pch=16, col='red', cex=2)
abline(a=0, b=1, col=rgb(0,0,0,0.5), lty=2)
dev.off()

# =========

taxa <- c(
  'Siphoviridae', 'Myoviridae', 'Geminiviridae', 'Podoviridae', 'Circoviridae',
  'Potyviridae', 'Papillomaviridae', 'Rhabdoviridae', 'Parvoviridae', 'Picornaviridae',
  'Flaviviridae', 'Tolecusatellitidae', 'Polyomaviridae', 'Betaflexiviridae', 'Genomoviridae',
  'Anelloviridae', 'Herpesviridae', 'Alphasatellitidae', 'Retroviridae', 'Herelleviridae'
)


png(file='nOverlapsPerORF.png', width=10*600, height=8*600, res=600)

par(mar=c(0,0,0,0), xpd=F, mfrow=c(4,5))

for (i in 1:length(taxa)) {
  tx <- taxa[i]
  plot(x[virus$Family != tx], y[virus$Family != tx], log='xy',
       cex=0.5, pch=16, col=rgb(.5,.5,.5, .2), yaxt='n')
  #axis(side=1, at=10^(1:7), label=10^(1:7))
  axis(side=2, at=c(0.5, 5, 50, 500), label=c(0, 5, 50, 500), las=2)
  
  points(x[virus$Family == tx], y[virus$Family == tx], 
         pch=(21:25)[(i-1)%%5+1], bg=pal[i], col='white', cex=1.5)
  text(x=1, y=300, label=tx, cex=1.5, adj=0)
}

dev.off()

