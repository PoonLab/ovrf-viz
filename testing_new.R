# File with Baltimore info
setwd('~/Projects/ovrf-review/data/Baltimore')
virus <- read.csv('species_file_baltimore2.csv')

# 1. Piechart with virus count according to Baltimore classification
pdf(file='/home/lmunoz/Projects/ovrf-review/Figures/pie_chart.pdf', width=12, height=9)
par(lwd=5)
my_table <- sort(table(virus$baltimore.class), decreasing = TRUE)
lbls <- paste(names(my_table), "\n", my_table, sep="")
pie(my_table, border = "black", labels = lbls, col = pal3
    , main = "Baltimore classification of data base")
dev.off()

pal3<- c("#1D6996", "#38A6A5", "#0F8554", "#73AF48", "#EDAD08", "#E17C05", "#CC503E", "#94346E")
#2. Plots of families versus ovrfs

par(lwd=2)
a <- virus$Number.of.proteins
b <- virus$n.overlaps
b[b==0] <- 0.5
b <- jitter(b)
a <- jitter(a)

pdf(file='/home/lmunoz/Projects/ovrf-review/Figures/Number_of_orfs.pdf', width=9, height=9)
plot(a,b, log='xy', type='n', col=as.numeric(virus$Molecule.type), cex=1, pch=as.numeric(virus$Topology), 
     xlab = 'Number of ORFs', ylab = 'Number of overlaps', main = "Virus data Base")

points( a[virus$baltimore.class%in%'ds_DNA'], b[virus$baltimore.class%in%'ds_DNA'], col=pal3[1], pch=5, cex = 1.5)
points( a[virus$baltimore.class%in%'Unknown'], b[virus$baltimore.class%in%'Unknown'], col=pal3[2], pch=8, cex = 1.5)
points( a[virus$baltimore.class%in%'ss_DNA'], b[virus$baltimore.class%in%'ss_DNA'], col=pal3[4], pch=15, cex = 1.2)
points( a[virus$baltimore.class%in%'ss_RNA_+'], b[virus$baltimore.class%in%'ss_RNA_+'], col=pal3[3], pch=17, cex = 1.2)
points( a[virus$baltimore.class%in%'ss_RNA_-'], b[virus$baltimore.class%in%'ss_RNA_-'], col=pal3[5], pch=20, cex = 1.5)
points( a[virus$baltimore.class%in%'ds_RNA'], b[virus$baltimore.class%in%'ds_RNA'], col=pal3[6], pch=18, cex = 1.5)
points( a[virus$baltimore.class%in%'RT_viruses'], b[virus$baltimore.class%in%'RT_viruses'], col=pal3[7], pch=20, cex = 1.5)
points( a[virus$baltimore.class%in%'circular_ss_RNA'], b[virus$baltimore.class%in%'circular_ss_RNA'], col=pal3[8], pch=4, cex = 1.5)

legend('topleft', legend= c('dsDNA', 'Unknown', 'ssDNA', '(+) ssRNA', '(-) ssRNA', 'dsRNA',  'RTviruses', '(circular) ssRNA'), 
       fill = pal3)
dev.off()

#3. Plots of length 
more_than_one <- subset(virus, n.overlaps > 0)
x <- more_than_one$Genome.length
y <- more_than_one$len.overlaps
y <- jitter(y)
x <- jitter(x)

pdf(file='/home/lmunoz/Projects/ovrf-review/Figures/ovrfs_length.pdf', width=9, height=9)
plot(x,y, log='xy', type='n', col=as.numeric(more_than_one$Molecule.type), cex=1, pch=as.numeric(more_than_one$Topology), 
     xlab = 'Number of ORFs', ylab = 'Number of overlaps', main = "more_than_one data Base")

points( x[more_than_one$baltimore.class%in%'ds_DNA'], y[more_than_one$baltimore.class%in%'ds_DNA'], col=pal3[1], pch=5, cex = 1.5)
points( x[more_than_one$baltimore.class%in%'Unknown'], y[more_than_one$baltimore.class%in%'Unknown'], col=pal3[2], pch=8, cex = 1.5)
points( x[more_than_one$baltimore.class%in%'ss_DNA'], y[more_than_one$baltimore.class%in%'ss_DNA'], col=pal3[4], pch=15, cex = 1.2)
points( x[more_than_one$baltimore.class%in%'ss_RNA_+'], y[more_than_one$baltimore.class%in%'ss_RNA_+'], col=pal3[3], pch=17, cex = 1.2)
points( x[more_than_one$baltimore.class%in%'ss_RNA_-'], y[more_than_one$baltimore.class%in%'ss_RNA_-'], col=pal3[5], pch=20, cex = 1.5)
points( x[more_than_one$baltimore.class%in%'ds_RNA'], y[more_than_one$baltimore.class%in%'ds_RNA'], col=pal3[6], pch=18, cex = 1.5)
points( x[more_than_one$baltimore.class%in%'RT_viruses'], y[more_than_one$baltimore.class%in%'RT_viruses'], col=pal3[7], pch=20, cex = 1.5)
points( x[more_than_one$baltimore.class%in%'circular_ss_RNA'], y[more_than_one$baltimore.class%in%'circular_ss_RNA'], col=pal3[8], pch=4, cex = 1.5)

legend('topleft', legend= c('dsDNA', 'Unknown', 'ssDNA', '(+) ssRNA', '(-) ssRNA', 'dsRNA',  'RTviruses', '(circular) ssRNA'), 
       fill = pal3)

dev.off()

#4. 

# Correcting ringeplot

library(ggplot2)
library(ggridges)
library(gridExtra)
#by length

pdf(file='/home/lmunoz/Projects/ovrf-review/Figures/Ridgeplots.pdf', width=9, height=6)
virus3<- virus[which(virus$n.overlaps!=0),]

ov.len<- ggplot(virus3, aes(x = log10(len.overlaps), y = baltimore.class, group = baltimore.class, fill = baltimore.class)) + 
  geom_density_ridges()+
  scale_fill_manual(values=pal3) +
  ggtitle("Average overlap length")+
  theme_bw()+
  labs(y="Category", x = "[log10] Average overlap length (nts)")


n.ov <-ggplot(virus3, aes(x = log10(n.overlaps), y = baltimore.class, group = baltimore.class, fill = baltimore.class)) + 
  geom_density_ridges( 
    aes(point_color = baltimore.class, point_fill = baltimore.class, point_shape = baltimore.class),
    alpha = .2, point_alpha = 1, jittered_points = TRUE
  )+
  scale_point_color_hue(l = 40) +
  scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22, 23, 24, 25, 21, 26))+
  ggtitle("Number of overlapping regions")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(y="Category", x = "[log10] Number of overlapping regions")


grid.arrange(ov.len, n.ov, nrow = 2)

dev.off()

### 
# Analizing the overlaps

pal <- c("#588C7E",
"#F2E394",
"#F2AE72",
"#D96459",
"#8C4646")

pal<- c("#eaa9bd","#dd88ac","#ca699d","#b14d8e","#91357d")

pal <- c("#F7894F", "#D6697B", "#BC69EB", "#5C5FD6", "#75C7E6")

pal2 <- c("#009392","#39b185","#9ccb86","#e9e29c","#eeb479","#e88471","#cf597e")

pal2 <- c("#39b185","#9ccb86","#e9e29c","#eeb479","#e88471")
  
pal3<- c("#1D6996", "#38A6A5", "#0F8554", "#73AF48", "#EDAD08", "#E17C05", "#CC503E", "#94346E")


#5F4690,#1D6996,#38A6A5,#0F8554,#73AF48,#EDAD08,#E17C05,#CC503E,#94346E,#6F4070,#994E95,#666666

setwd('~/Projects/ovrf-review/data/')
# all overlapping genes in all genomes (multiple overlaps in genes with "introns" are combined into one entry)
overlaps <- read.csv('ovrfs-reduced.csv')

pdf(file='/home/lmunoz/Projects/ovrf-review/Figures/shift.pdf', width=9, height=6)

pal <- c("#F7894F", "#D6697B", "#BC69EB", "#5C5FD6", "#75C7E6")

ov<-ggplot(overlaps, aes(x = log10(overlap), y = shift , group = shift, fill=as.factor(shift))) + 
  geom_density_ridges()+
  scale_fill_manual(values=pal) +
  ggtitle("Average overlap length")+
  #theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme_bw()+
  labs(y="Frame shift", x = "[log10] Overlap length (nts)")

ov.len<- ggplot(virus3, aes(x = log10(len.overlaps), y = baltimore.class, group = baltimore.class, fill = baltimore.class)) + 
  geom_density_ridges()+
  scale_fill_manual(values=pal3) +
  ggtitle("Average overlap length")+
  theme_bw()+
  labs(y="Category", x = "[log10] Average overlap length (nts)")

grid.arrange(ov, ov.len, nrow = 2)

dev.off()


pdf(file='/home/lmunoz/Projects/ovrf-review/Figures/pie_shift.pdf', width=9, height=6)
my_table <- table(overlaps$shift)
lbls <- paste(names(my_table), "\n", my_table, sep="")
pie<-pie(my_table, border = "black", labels = lbls, col = pal  , main = "Frame shifts")

dev.off()


###
# Virus, overlaps
index <- match(overlaps$accn, virus$Accession)
overlaps$baltimore.class <- virus$baltimore.class[index]

counts <- sort(table(overlaps$shift, overlaps$baltimore.class))
matrix <- as.matrix(counts)
#TODO legend with most abundant families
barplot(matrix, col = pal, offset = 0)
legend("right",inset=c(-0.2,0), legend=rownames(matrix), fill= pal)

plus2<- subset(overlaps, shift == 1)
test<- subset(overlaps, baltimore.class == "dsRNA")

pdf(file='/home/lmunoz/Projects/ovrf-review/Figures/barplots.pdf', width=9, height=6)
par(mfrow=c(2,3))

n<-c("#18596B", "#569DB0", "#6FCAE3", "#F5D485", "#D9AA5F")


balt<- names(table(overlaps$baltimore.class))[2:7]
shift<- names(table(overlaps$shift))
blah <- lapply(balt, function(i){
  shift<- names(table(overlaps$shift))
  x <- subset(overlaps, baltimore.class == i)
  #barplot(sort(table(x$shift), decreasing = TRUE), col=pal3)
  barplot(table(x$shift), col=n, main = i)
})

dev.off()



library(plyr)

new_p<-c("#4D7A76", "#3EA89E", "#68DB8F", "#E082A8", "#A83D9F")
par(mfrow=c(3,4))
names <- names(table(overlaps$shift))
blah <- lapply(names, function(i){
  x <- subset(overlaps, shift == i)
  #barplot(sort(table(x$shift), decreasing = TRUE), col=pal3)
  bp<- barplot(sort(table(x$baltimore.class)), col=new_p, main = i, axes=FALSE)
  
  axis(1, at=bp, labels = names)
})

dev.off()

## Mosaic plot
mosaicplot(table(overlaps$baltimore.class, overlaps$shift))
