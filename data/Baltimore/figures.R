# File with Baltimore info
setwd('~/Projects/ovrf-review/data/Baltimore')
virus <- read.csv('species_file_baltimore2.csv')
plot(table(virus_edited$baltimore.class))


a <- virus$Number.of.proteins
b <- virus$n.overlaps
b[b==0] <- 0.5
b <- jitter(b)
a <- jitter(a)


# FantasticFox1 = c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20")
#GrandBudapest2 = c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
#fill = c('yellowgreen', 'plum4', 'yellow2', 'darkorange', 'firebrick', 'palevioletred', 'dodgerblue3'))
###################################################
# Plot according to Molecule PAGE 1
###################################################

pdf(file='full_virus.pdf', width=12, height=9, title = "Analysis of overlappin genes in viral genomes")
par(mar = c(4,4,3,4))

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
par(lwd = 2)

plot(a,b, log='xy', type='n', col=as.numeric(virus$Molecule.type), cex=1, pch=as.numeric(virus$Topology), 
     xlab = 'Number of ORFs', ylab = 'Number of overlaps', main = "Virus data Base")

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

par(lwd = 2)
my_table <- table(virus$baltimore.class)
lbls <- paste(names(my_table), "\n", my_table, sep="")
pie(my_table, border = "black", labels = lbls, col = c('plum4', 'yellowgreen', 'yellow2', 'darkorange', 'firebrick', 'palevioletred', 'dodgerblue3', "#49DA9A")
    , main = "Baltimore classification of data base")
#pie(my_table, border = NA, density = 80, labels = lbls, col = c('plum4', 'yellowgreen', 'yellow2', 'darkorange', 'firebrick', 'palevioletred', 'dodgerblue3', "#49DA9A"))

#dotchart(my_table, cex = 0.8, col = c('plum4', 'yellowgreen', 'yellow2', 'darkorange', 'firebrick', 'palevioletred', 'dodgerblue3', "#49DA9A"),
         #main = "Baltimore classification of data base")

unknown <- virus[which(virus$baltimore.class=='Unknown'),]
unk_table <- table(unknown$Molecule.type)
lbls <- paste(names(my_table), "\n", my_table, sep="")
dotchart(unk_table, col=c('yellow2', 'plum4', 'yellowgreen', 'darkorange', 'dodgerblue3', "#49DA9A", 'palevioletred'),
         main = "Molecules of unknown Baltimore viruses", labels = lbls)


pal <- c("#f6f078", "#01d28e", "#434982", "#730068")
#pal <- c("#f1d4d4", "#ddb6c6", "#ac8daf", "#484c7f")
more_than_one <- subset(virus, n.overlaps > 0)

subset1 <- subset(virus, n.overlaps == 0)
subset2 <- subset(virus, n.overlaps == 1)
subset3 <- subset(virus, n.overlaps == 2)
subset4 <- subset(virus, n.overlaps == 3)
subset5 <- subset(virus, n.overlaps == 4)
subset6 <- subset(virus, n.overlaps == 5)
subset7 <- subset(virus, n.overlaps == 6)
subset8 <- subset(virus, n.overlaps == 7)
subset9 <- subset(virus, n.overlaps == 8)
subset10 <- subset(virus, n.overlaps == 9)
subset11 <- subset(virus, 10 < n.overlaps & n.overlaps <= 50)
subset12 <- subset(virus, 50<n.overlaps & n.overlaps<=100)
subset13 <- subset(virus, 100<n.overlaps & n.overlaps<=500)
subset14 <- subset(virus, 500<n.overlaps)

matrix_overlap <- matrix(c(dim(subset1)[1], dim(subset2)[1], dim(subset3)[1], dim(subset4)[1], dim(subset5)[1], dim(subset6)[1], 
                           dim(subset7)[1], dim(subset8)[1], dim(subset9)[1], dim(subset10)[1], dim(subset11)[1], dim(subset12)[1], dim(subset13)[1], dim(subset14)[1]),
                         ncol = 14)
colnames(matrix_overlap) <- c("0","1", "2", "3", "4", "5", "6", "7", "8", "9", "10 to 50", "50 to 100", "100 to 500", "> 500")
dotchart(matrix_overlap[1,], main = "Number of overlaps", cex = 0.8, col=pal)


####
# Outliers

# Entries with no overlap 
no_ov <- subset(virus, n.overlaps==0)
# No overlap with only one protein
one_prot <- subset(virus, Number.of.proteins==1)
# No overlaps, more than one protein
more_than_12 <- subset(no_ov, Number.of.proteins > 12)

###################################################
# Viruses with no overlaps PAGE 2
###################################################
#pal <- c("#f1d4d4", "#ddb6c6", "#ac8daf", "#484c7f")
pal<- c("#293462", "#00818a", "#ec9b3b", "#f7be16")

# Displaying information for viruses with no overlaps
par(mfrow=c(2,2))
pie(table(no_ov$family), main = "Families of viruses with no overlaps", cex = 0.8, col=pal)
dotchart(table(no_ov$Number.of.proteins), main = "Number of proteins", cex = 0.8, col=pal)
pie(table(one_prot$family), main = "Families with one protein", cex = 0.8, col=pal)
pie(table(more_than_12$family), main =  "Families with more than 12 proteins", cex = 0.8, col=pal)
#mtext("Viruses with no overlaps", side=3, outer=TRUE, line=-2, cex = 1.2)


###################################################
# Viruses with at least one overlap PAGE 3
###################################################
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))

# Displaying information for viruses with more then zero overlap 
more_than_one <- subset(virus, n.overlaps > 0)

pie(table(more_than_one$family), main="Families with at least one overlap", cex = 0.8, col=pal)
#pie(table(more_than_one$n.overlaps), main = "Number of overlaps", cex = 0.8, col=pal)

#pie(table(more_than_one$Number.of.proteins), main = "Number of proteins", cex = 0.8, col=pal)


subset2 <- subset(more_than_one, Number.of.proteins == 1)
subset3 <- subset(more_than_one, Number.of.proteins == 2)
subset4 <- subset(more_than_one, Number.of.proteins == 3)
subset5 <- subset(more_than_one, Number.of.proteins == 4)
subset6 <- subset(more_than_one, Number.of.proteins == 5)
subset7 <- subset(more_than_one, Number.of.proteins == 6)
subset8 <- subset(more_than_one, Number.of.proteins == 7)
subset9 <- subset(more_than_one, Number.of.proteins == 8)
subset10 <- subset(more_than_one, Number.of.proteins == 9)
subset11 <- subset(more_than_one, 10 < Number.of.proteins & Number.of.proteins <= 50)
subset12 <- subset(more_than_one, 50<Number.of.proteins & Number.of.proteins<=100)
subset13 <- subset(more_than_one, 100<Number.of.proteins & Number.of.proteins<=500)
subset14 <- subset(more_than_one, 500<Number.of.proteins)

matrix_proteins <- matrix(c(dim(subset2)[1], dim(subset3)[1], dim(subset4)[1], dim(subset5)[1], dim(subset6)[1], 
                           dim(subset7)[1], dim(subset8)[1], dim(subset9)[1], dim(subset10)[1], dim(subset11)[1], dim(subset12)[1], dim(subset13)[1], dim(subset14)[1]),
                         ncol = 13)
colnames(matrix_proteins) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10 to 50", "50 to 100", "100 to 500", "> 500")
dotchart(matrix_proteins[1,], main = "Number of Proteins", cex = 0.8, col=pal)



# Medium len of the overlap
#hist(more_than_one$len.overlaps)

subset2 <- subset(more_than_one, 1 <= len.overlaps & len.overlaps <= 2)
subset3 <- subset(more_than_one, 2 < len.overlaps & len.overlaps <= 3)
subset4 <- subset(more_than_one, 3 < len.overlaps & len.overlaps <= 4)
subset5 <- subset(more_than_one, 4 < len.overlaps & len.overlaps <= 5)
subset6 <- subset(more_than_one, 5 < len.overlaps & len.overlaps <= 6)
subset7 <- subset(more_than_one, 6 < len.overlaps & len.overlaps <= 7)
subset8 <- subset(more_than_one, 7 < len.overlaps & len.overlaps <= 8)
subset9 <- subset(more_than_one, 8 < len.overlaps & len.overlaps <= 9)
subset10 <- subset(more_than_one, 9 < len.overlaps & len.overlaps <= 10)
subset11 <- subset(more_than_one, 10 < len.overlaps & len.overlaps <= 11)
subset12 <- subset(more_than_one, 11 < len.overlaps & len.overlaps <= 50)
subset13 <- subset(more_than_one, 50 < len.overlaps & len.overlaps <= 100)
subset14 <- subset(more_than_one, 100 <len.overlaps & len.overlaps<=500)
subset15 <- subset(more_than_one, 500 <len.overlaps & len.overlaps<=1000)
subset16 <- subset(more_than_one, 1000 <len.overlaps)

matrix_proteins <- matrix(c(dim(subset2)[1], dim(subset3)[1], dim(subset4)[1], dim(subset5)[1], dim(subset6)[1], 
                            dim(subset7)[1], dim(subset8)[1], dim(subset9)[1], dim(subset10)[1], dim(subset11)[1], dim(subset12)[1],
                            dim(subset13)[1], dim(subset14)[1], dim(subset15)[1], dim(subset16)[1]), ncol = 15)
colnames(matrix_proteins) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11 to 50", "50 to 100", "100 to 500", "500to 1k", " > 1k ")
dotchart(matrix_proteins[1,], main = "Overlap lenght", cex = 0.8, col=pal)

#mtext("Viruses with at least one overlap", side=3, outer=TRUE, line=-2, cex = 1.2)

###################################################
# Overlaps longer than 3000 nucleotides PAGE 4
###################################################
#par(mfrow=c(1,2), mar=c(5,5,1,3))
#pdf(file='long_overlaps.pdf', width=8, height=11)

pal<- c("#293462", "#00818a", "#ec9b3b", "#f7be16")

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE),
       widths=c(5,5), heights=c(5,5))
long_ov<-subset(virus, len.overlaps > 1000)
pie(table(long_ov$family), main = "Families with overlaps larger than 1000 nt", col = pal, cex = 0.8)
dotchart(table(long_ov$n.overlaps), main = "Number of overlaps", col = pal)


x <- long_ov$Genome.length
y <- long_ov$len.overlaps
y <- jitter(y)
x <- jitter(x)

# Plot according to Molecule
par(lwd = 2)
plot(x,y, type='n', col=as.numeric(long_ov$Molecule.type), cex=1.5, pch=as.numeric(long_ov$Topology), 
     xlab = 'Genome Length', ylab = 'Medium overlap Lenght', main = "Long overlaps")
points( x[long_ov$baltimore.class%in%'ss_RNA_+'], y[long_ov$baltimore.class%in%'ss_RNA_+'], col='palevioletred', pch=17, cex = 1.5)
points( x[long_ov$baltimore.class%in%'ds_DNA'], y[long_ov$baltimore.class%in%'ds_DNA'], col='yellowgreen', pch=5, cex = 1.5)
points( x[long_ov$baltimore.class%in%'Unknown'], y[long_ov$baltimore.class%in%'Unknown'], col='plum4', pch=19, cex = 1.5)
points( x[long_ov$baltimore.class%in%'ds_RNA'], y[long_ov$baltimore.class%in%'ds_RNA'], col='yellow2', pch=8, cex = 1.5)
points( x[long_ov$baltimore.class%in%'ss_DNA'], y[long_ov$baltimore.class%in%'ss_DNA'], col='darkorange', pch=18, cex = 1.5)
points( x[long_ov$baltimore.class%in%'ss_RNA_-'], y[long_ov$baltimore.class%in%'ss_RNA_-'], col='firebrick', pch=19, cex = 1.5)
points( x[long_ov$baltimore.class%in%'RT_viruses'], y[long_ov$baltimore.class%in%'RT_viruses'], col='dodgerblue3', pch=16, cex = 1.5)
points( x[long_ov$baltimore.class%in%'circular_ss_RNA'], y[long_ov$baltimore.class%in%'circular_ss_RNA'], col="#49DA9A", pch=4, cex = 1.5)

legend('topleft', legend= c('(+) ssRNA', 'dsDNA', 'Unknown', 'dsRNA', 'ssDNA', '(-) ssRNA', 'RTviruses', '(circular) ssRNA'), 
       fill = c('palevioletred', 'yellowgreen', 'plum4', 'yellow2', 'darkorange', 'firebrick', 'dodgerblue3', "#49DA9A"))

dev.off()

#########
