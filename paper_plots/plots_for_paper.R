library(RColorBrewer)

setwd('/home/lmunoz/Projects/ovrf-review/data/Baltimore')
virus <- read.csv('species_file_baltimore2.csv')

# 1. Baltimore classification of the data base 
sub <- subset(virus, family = 'Unknown')
my_table <- table(sub$baltimore.class)
lbls <- paste(names(my_table), "\n", my_table, sep="")
pie(my_table, border = "black", labels = lbls, main = "Baltimore classification of data base")


# 2. Summary of overlapping regions in the viral genomes
# a. Number of ORFS vs Number of overlaps
a <- virus$Number.of.proteins
b <- virus$n.overlaps
b[b==0] <- 0.5
b <- jitter(b)
a <- jitter(a)
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
par(lwd = 2)
plot(a,b, log='xy', type='n', col=as.numeric(virus$Molecule.type), cex=1, pch=as.numeric(virus$Topology), 
     xlab = '[log10] Number of ORFs', ylab = '[log10] Number of overlaps', main = "Virus data Base")

points( a[virus$baltimore.class%in%'ds_DNA'], b[virus$baltimore.class%in%'ds_DNA'], col='yellowgreen', pch=5)
points( a[virus$baltimore.class%in%'Unknown'], b[virus$baltimore.class%in%'Unknown'], col='plum4', pch=17)
points( a[virus$baltimore.class%in%'ds_RNA'], b[virus$baltimore.class%in%'ds_RNA'], col='yellow2', pch=8)
points( a[virus$baltimore.class%in%'ss_DNA'], b[virus$baltimore.class%in%'ss_DNA'], col='darkorange', pch=18)
points( a[virus$baltimore.class%in%'ss_RNA_-'], b[virus$baltimore.class%in%'ss_RNA_-'], col='firebrick', pch=19)
points( a[virus$baltimore.class%in%'ss_RNA_+'], b[virus$baltimore.class%in%'ss_RNA_+'], col='palevioletred', pch=17)
points( a[virus$baltimore.class%in%'RT_viruses'], b[virus$baltimore.class%in%'RT_viruses'], col='dodgerblue3', pch=16)
points( a[virus$baltimore.class%in%'circular_ss_RNA'], b[virus$baltimore.class%in%'circular_ss_RNA'], col="#49DA9A", pch=4)

legend('topleft', legend= c('dsDNA', 'Unknown', 'dsRNA', 'ssDNA', '(-) ssRNA', '(+) ssRNA', 'RTviruses', '(circular) ssRNA'), 
       fill = c( 'yellowgreen', 'plum4', 'yellow2', 'darkorange', 'firebrick', 'palevioletred', 'dodgerblue3', "#49DA9A"))

# b. Genome length vs mean nucleotides involved in an overlap
more_than_one <- virus
x <- more_than_one$Genome.length
y <- more_than_one$len.overlaps
y <- jitter(y)
x <- jitter(x)

plot(x,y, type='n', log = 'xy', col=as.numeric(more_than_one$Molecule.type), cex=2, pch=as.numeric(more_than_one$Topology), 
     xlab = '[log10] Genome Length', ylab = '[log10] Overlap Lenght', main = "Overlap lenght")
points( x[more_than_one$baltimore.class%in%'ss_RNA_+'], y[more_than_one$baltimore.class%in%'ss_RNA_+'], col='palevioletred', pch=15, cex = 0.8)
points( x[more_than_one$baltimore.class%in%'ds_DNA'], y[more_than_one$baltimore.class%in%'ds_DNA'], col='yellowgreen', pch=5, cex = 0.8)
points( x[more_than_one$baltimore.class%in%'Unknown'], y[more_than_one$baltimore.class%in%'Unknown'], col='plum4', pch=1, cex = 0.8)
points( x[more_than_one$baltimore.class%in%'ds_RNA'], y[more_than_one$baltimore.class%in%'ds_RNA'], col='yellow2', pch=8, cex = 0.8)
points( x[more_than_one$baltimore.class%in%'ss_DNA'], y[more_than_one$baltimore.class%in%'ss_DNA'], col='darkorange', pch=18, cex = 0.8)
points( x[more_than_one$baltimore.class%in%'ss_RNA_-'], y[more_than_one$baltimore.class%in%'ss_RNA_-'], col='firebrick', pch=20, cex = 0.8)
points( x[more_than_one$baltimore.class%in%'RT_viruses'], y[more_than_one$baltimore.class%in%'RT_viruses'], col='dodgerblue3', pch=2, cex = 0.8)
points( x[more_than_one$baltimore.class%in%'circular_ss_RNA'], y[more_than_one$baltimore.class%in%'circular_ss_RNA'], col="#49DA9A", pch=4, cex = 0.8)

legend('topright', legend= c('(+) ssRNA', 'dsDNA', 'Unknown', 'dsRNA', 'ssDNA', '(-) ssRNA', 'RTviruses', '(circular) ssRNA'), 
       fill = c('palevioletred', 'yellowgreen', 'plum4', 'yellow2', 'darkorange', 'firebrick', 'dodgerblue3', "#49DA9A"))

# 3. 


