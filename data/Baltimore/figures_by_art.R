# File with Baltimore info
setwd('/home/lmunoz/Projects/ovrf-review/data/Baltimore')
virus <- read.csv('species_file_baltimore2.csv')

hist(log10(virus$len.overlaps))
sum(virus$len.overlaps==0, na.rm=TRUE)

#Density plot
x <- log10(virus$len.overlaps)
plot(density(x, na.rm=T), type='n', main='', ylim=c(0,1.5))
pal <- rainbow(7, v=0.8)
count <- 1
for (bc in unique(virus$baltimore.class)) {
  if (bc == 'circular_ss_RNA') next;
  lines(density(x[virus$baltimore.class==bc], na.rm=T, bw=0.2), 
        col=pal[count], lwd=2)
  count <- count + 1
}


## collect mean overlap length by family
tab <- table(virus$family) 
family <- names(tab)[tab > 100]
family <- family[!grepl("unclassified|Unknown", family)]

virus2 <- virus[is.element(virus$family, family), ]
virus2$baltimore.class <- factor(virus2$baltimore.class)
virus2$family <- factor(virus2$family)
# TODO: get equiprobable cut points
cut.pts <- cut2(log10(virus2$len.overlaps), g=5, onlycuts=T)

m <- sapply(split(virus2, virus2$family), function(fam) {
  h <- hist(log10(fam$len.overlaps), breaks=cut.pts, plot=F)
  y <- c(sum(is.na(fam$len.overlaps)), h$counts)
  y/sum(y)
})
stars(t(m))

##################################
# Correcting ringeplot

library(ggplot2)
library(ggridges)
library(gridExtra)
#by length
virus3<- virus[which(virus$n.overlaps!=0),]

ov.len<- ggplot(virus3, aes(x = log10(len.overlaps), y = baltimore.class, group = baltimore.class, fill = baltimore.class)) + 
   geom_density_ridges( 
     aes(point_color = baltimore.class, point_fill = baltimore.class, point_shape = baltimore.class),
                        alpha = .3, point_alpha = 1, jittered_points = TRUE
   )+
  scale_point_color_hue(l = 40) +
  scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22, 23, 24, 25, 21, 26))+
  ggtitle("Average overlap length")+
  theme(plot.title = element_text(hjust = 0.5))+
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



# Blah
sub <- subset(virus, family = 'Unknown')
sub2<-(subset(sub, Genome.length > 5000 & Genome.length < 12000))

my_table <- table(sub$Number.of.proteins)
lbls <- paste(names(my_table), "\n", my_table, sep="")
pie(my_table, border = "black", labels = lbls, main = "Baltimore classification of data base")
