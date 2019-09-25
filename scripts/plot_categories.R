library('dplyr')

setwd('~/Projects/ovrf-review/data/')
# All virus information 
virus <- read.csv('virus_with_info.csv')
# all overlapping genes in all genomes (multiple overlaps in genes with "introns" are combined into one entry)
overlaps <- read.csv('ovrfs-reduced.csv')
# listing of reading frames in all genomes
orfs <- read.csv('orfs-fixed.csv')


# gsub() replaces all matches of a string, if the parameter is a string vector, returns a string vector f the same length and with the same attributes
virus$Genome.length <- as.integer(gsub(' nt$', '', virus$Genome.length))

virus$Number.of.proteins[virus$Number.of.proteins=='-'] <- NA
virus$Number.of.proteins <- as.integer(as.character(
  virus$Number.of.proteins))


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

#START PLOT
pdf(file='viruses_colored2.pdf', width=8, height=11)
par(mfrow=c(3,2), mar=c(5,5,1,1))
############################
x <- virus$Genome.length
y <- virus$n.overlaps
y[y==0] <- 0.5
molecules= c("DNA", "RNA")
plot(x,y, log = 'xy', type='n', xlab="Genome length", ylab = "Number of overlaps")
'%ni%' <- Negate('%in%')
points( x[virus$Molecule.type%ni%molecules], y[virus$Molecule.type%ni%molecules], col='yellow2', pch=15)
points( x[virus$Molecule.type%in%'DNA'], y[virus$Molecule.type%in%'DNA'], col='yellowgreen', pch=16)
points( x[virus$Molecule.type%in%'RNA'], y[virus$Molecule.type%in%'RNA'], col='plum4', pch=17)
legend("topleft", legend=c('Other',"DNA", "RNA"), fill =c('yellow2', 'yellowgreen','plum4'))


###############################

a <- virus$Number.of.proteins
b <- virus$n.overlaps
b[b==0] <- 0.5
b <- jitter(b)
a <- jitter(a)

###########################
# Plot according to Molecule
plot(a,b, log='xy', type='n', col=as.numeric(virus$Molecule.type), cex=1, pch=as.numeric(virus$Topology), 
     xlab = 'Number of ORFs', ylab = 'Number of overlaps')
molecules= c("DNA", "RNA")
'%ni%' <- Negate('%in%')
points( a[virus$Molecule.type%ni%molecules], b[virus$Molecule.type%ni%molecules], col='yellow2', pch=15)
points( a[virus$Molecule.type%in%'DNA'], b[virus$Molecule.type%in%'DNA'], col='yellowgreen', pch=16)
points( a[virus$Molecule.type%in%'RNA'], b[virus$Molecule.type%in%'RNA'], col='plum4', pch=17)
legend("topleft", legend=c('Other',"DNA", "RNA"), fill =c('yellow2', 'yellowgreen','plum4'))

###########################
#Plot according to Host
plot(a,b, log='xy', type='n', col=as.numeric(virus$Host), cex=1, pch=as.numeric(virus$Topology), 
     xlab = 'Number of ORFs', ylab = 'Number of overlaps')

points( a[virus$Host%in%'bacteria'], b[virus$Host%in%'bacteria'], col='yellowgreen', pch=4)
points( a[virus$Host%in%'vertebrates'], b[virus$Host%in%'vertebrates'], col='plum4', pch=5)
points( a[virus$Host%in%'plants'], b[virus$Host%in%'plants'], col='yellow2', pch=8)
points( a[virus$Host%in%'invertebrates'], b[virus$Host%in%'invertebrates'], col='darkorange', pch=16)
points( a[virus$Host%in%'fungi'], b[virus$Host%in%'fungi'], col='firebrick', pch=18)
points( a[virus$Host%in%'archaea'], b[virus$Host%in%'archaea'], col='palevioletred', pch=19)
points( a[virus$Host%in%'protozoa'], b[virus$Host%in%'protozoa'], col='dodgerblue3', pch=17)

legend('topleft', legend= c('Bacteria', 'Vertebrates', 'Plants', 'Invertebrates', 'Fungi', 'Archea', 'Protozoa'), 
       fill = c('yellowgreen', 'plum4', 'yellow2', 'darkorange', 'firebrick', 'palevioletred', 'dodgerblue3'))
#############################
# By topology
plot(a,b, log='xy', type='n', col=as.numeric(virus$Molecule.type), cex=1, pch=as.numeric(virus$Topology), 
     xlab = 'Number of ORFs', ylab = 'Number of overlaps')
points( a[virus$Topology%in%'linear'], b[virus$Topology%in%'linear'], col='yellowgreen', pch=16)
points( a[virus$Topology%in%'circular'], b[virus$Topology%in%'circular'], col='plum4', pch=17)
legend("topleft", legend=c('Linear',"Circular"), fill=c('yellowgreen','plum4'))


############################
# BY MOLECULE
x<- virus$Genome.length
y <-virus$len.overlaps
plot(x,y, log = 'xy', type='n', xlab = 'Genome length', ylab='Overalps lenght')

points( x[virus$Molecule.type%in%'DNA'], y[virus$Molecule.type%in%'DNA'], col='yellowgreen', pch=16)
points( x[virus$Molecule.type%in%'RNA'], y[virus$Molecule.type%in%'RNA'], col='plum4', pch=17)
points( x[virus$Molecule.type%ni%molecules], y[virus$Molecule.type%ni%molecules], col='yellow2', pch=5)
legend("topleft", legend=c('Other',"DNA", "RNA"), fill =c('yellow2', 'yellowgreen','plum4'))
abline(a=0, b=1, col="red", lty=2, lwd=3)

#BY HOST
plot(x,y, log = 'xy', type='n', xlab = 'Genome length', ylab='Overalps lenght')
points( x[virus$Host%in%'bacteria'], y[virus$Host%in%'bacteria'], col='yellowgreen', pch=4)
points( x[virus$Host%in%'vertebrates'], y[virus$Host%in%'vertebrates'], col='plum4', pch=5)
points( x[virus$Host%in%'plants'], y[virus$Host%in%'plants'], col='yellow2', pch=8)
points( x[virus$Host%in%'invertebrates'], y[virus$Host%in%'invertebrates'], col='darkorange', pch=16)
points( x[virus$Host%in%'fungi'], y[virus$Host%in%'fungi'], col='firebrick', pch=18)
points( x[virus$Host%in%'archaea'], y[virus$Host%in%'archaea'], col='palevioletred', pch=19)
points( x[virus$Host%in%'protozoa'], y[virus$Host%in%'protozoa'], col='dodgerblue3', pch=17)
abline(a=0, b=1, col="red", lty=2, lwd=3)
legend('topleft', legend= c('Bacteria', 'Vertebrates', 'Plants', 'Invertebrates', 'Fungi', 'Archea', 'Protozoa'), 
       fill = c('yellowgreen', 'plum4', 'yellow2', 'darkorange', 'firebrick', 'palevioletred', 'dodgerblue3'))

dev.off()
############################


# Checking extreme entries
too_long <- virus %>% filter(virus$Number.of.proteins > 500)
too_short<- virus %>% filter(virus$Number.of.proteins <= 1)

long_overlap <- virus%>% filter(virus$len.overlaps >1000)

# Checking outliers
one_protein <- virus[which(virus$Number.of.proteins == 1),]
one_p_taxonomy <- one_protein$Taxonomy

no_special_ch <- function(x) {
  x<-as.character(x)
  xx<-gsub("\'|\\[|\\]", "", x)
  y <- strsplit(xx, ", ")
  return(y)
}
out <-sapply(one_p_taxonomy, no_special_ch)


x<-as.character(one_p_taxonomy[5])
x<-gsub("\'|\\[|\\]", "", xx)

sapply(1:nrow(virus), function(i) {
  row <- virus[i,]
  (row$n.overlaps>row$Genome.length)
})

subset(virus, n.overlaps>Genome.length)

sub <- subset(virus, n.overlaps < 1)
sub <- subset(virus, n.overlaps > Number.of.proteins)


