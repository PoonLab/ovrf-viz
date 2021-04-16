library(RColorBrewer)
options(scipen=999)
library(extrafont)
font_import()

######################################################
# Preparing information and color palette 
######################################################
setwd('/home/lmunoz/Projects/ovrf-review/data/Baltimore')
virus <- read.csv('species_file_baltimore2.csv', stringsAsFactors = F)

#Re naming Unknown according to molecule type
virus$baltimore.class[which(virus$baltimore.class == "Unknown" & grepl("RNA", virus$Molecule.type))] <- "unknown_RNA"
virus$baltimore.class[which(virus$baltimore.class == "Unknown" & grepl("DNA", virus$Molecule.type))] <- "unknown_DNA"

# Setting color palette
bal.class <- c("ds_RNA", "ss_RNA_+", "ss_RNA_-", "unknown_RNA", "ds_DNA", "ss_DNA", "unknown_DNA", "RT_viruses", "circular_ss_RNA")
col<- c("#1d6996", "#38a6a5", "#57a1b4", "#87bfdb", "#edad08", "#e17c05", "#cc503e", "#808080", "#73af48")
pal <- as.data.frame(cbind(bal.class, col))

#######################################################
# Barplot database distribution
#######################################################
# tbl <- as.data.frame(table(virus$baltimore.class))
# tbl$bal.class <- factor(tbl$Var1, levels=c("ds_DNA", "ss_DNA", "unknown_DNA", "ds_RNA", "ss_RNA_+", "ss_RNA_-", "unknown_RNA", "RT_viruses" ))
# p <- ggplot(data=tbl, aes(x=Var1, y=Freq, fill=Var1)) +
#   geom_bar(stat="identity", width = 0.5)+
#   scale_fill_manual(values = c("circular_ss_RNA" = "#73af48", "ds_RNA" = "#1d6996", "ss_RNA_+" = "#38a6a5", "ss_RNA_-" = "#57a1b4", "unknown_RNA" = "#87bfdb",
#                                "ds_DNA" = "#edad08",  "ss_DNA" = "#e17c05", "unknown_DNA" = "#cc503e", "RT_viruses" = "#808080", "circular_ss_RNA" = "#73af48" ))+
#   
#   geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.3)+
#   labs(fill = "Baltimore class")+ 
#   theme_bw()+
#   ggtitle("Database distribution")
# p

# Try to create my own
#vector of colors based on baltimore class
bal.col <- merge(tbl[,c("bal.class", "Freq")], pal, by="bal.class") # freq of baltimore classification with assigned colors 
my_order <- rev(c("ds_DNA", "ss_DNA", "unknown_DNA", "ds_RNA", "ss_RNA_+", "ss_RNA_-", "unknown_RNA", "RT_viruses" ))
bal.col <- bal.col[match(my_order,bal.col$bal.class),]

pdf(file="/home/lmunoz/Projects/ovrf-review/paper_plots/database_distribution.jpg")
par(mar = c(6, 10, 4, 4))
barplot(height=bal.col$Freq, main = "Database distribution", xlab = "Number of entries",
        names = bal.col$bal.class, 
        col = as.character(bal.col$col),
        horiz = T,
        las=1, 
        xlim=c(0,4000), 
        border = F
        )
text(y = bp, x =bal.col$Freq, pos=4, labels = bal.col$Freq)
box(lwd=2, col="gray80")
dev.off()

#######################################################
# Plotting overlapping summary
#######################################################
loadfonts()
pdf(file="general_trends.pdf", family="Amiri", width=14, height=6.5)
# a. Number of ORFS vs Number of overlaps
x <- virus$Number.of.proteins
y <- virus$n.overlaps
y[y==0] <- 0.5
y <- jitter(y)
x <- jitter(x)

#layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
par(mfrow=c(1,2))
par(mar = c(6, 6, 4, 2))
plot(x,y, log='xy', type='n', col=as.numeric(virus$Molecule.type), cex=1, pch=as.numeric(virus$Topology), 
     xlab = '[log10] Number of ORFs', ylab = '[log10] Number of overlaps', main = "Virus data Base", cex.lab = 1.5)
lwd <- 3
cex <- 1
points( x[virus$baltimore.class%in%'ds_DNA'], y[virus$baltimore.class%in%'ds_DNA'], col=as.character(pal$col[pal$bal.class=="ds_DNA"]), pch=17, cex = cex)
points( x[virus$baltimore.class%in%'ss_RNA_+'], y[virus$baltimore.class%in%'ss_RNA_+'], col=as.character(pal$col[pal$bal.class=="ss_RNA_+"]), pch=4, cex = cex, lwd = lwd)
points( x[virus$baltimore.class%in%'ss_DNA'], y[virus$baltimore.class%in%'ss_DNA'], col=as.character(pal$col[pal$bal.class=="ss_DNA"]), pch=18, cex = cex)
points( x[virus$baltimore.class%in%'unknown_RNA'], y[virus$baltimore.class%in%'unknown_RNA'], col=as.character(pal$col[pal$bal.class=="unknown_RNA"]), pch=5, cex = cex, lwd = lwd)
points( x[virus$baltimore.class%in%'ds_RNA'], y[virus$baltimore.class%in%'ds_RNA'], col=as.character(pal$col[pal$bal.class=="ds_RNA"]), pch=8, cex = cex, lwd = 2)
points( x[virus$baltimore.class%in%'ss_RNA_-'], y[virus$baltimore.class%in%'ss_RNA_-'], col=as.character(pal$col[pal$bal.class=="ss_RNA_-"]), pch=6, cex = cex)
points( x[virus$baltimore.class%in%'unknown_DNA'], y[virus$baltimore.class%in%'unknown_DNA'], col=as.character(pal$col[pal$bal.class=="unknown_DNA"]), pch=20, cex = cex)
points( x[virus$baltimore.class%in%'RT_viruses'], y[virus$baltimore.class%in%'RT_viruses'], col=as.character(pal$col[pal$bal.class=="RT_viruses"]), pch=1, cex = cex, lwd = lwd)
points( x[virus$baltimore.class%in%'circular_ss_RNA'], y[virus$baltimore.class%in%'circular_ss_RNA'], col=as.character(pal$col[pal$bal.class=='circular_ss_RNA']), pch=4, cex = cex)


# b. Genome length vs mean nucleotides involved in an overlap
x <- virus$Genome.length
y <- virus$len.overlaps
y[y==0] <- 0.5
y <- jitter(y)
x <- jitter(x)
par(xpd = T, mar = c(6, 5, 4, 11))
plot(x,y, type='n', log = 'xy', col=as.numeric(virus$Molecule.type), cex=1, pch=as.numeric(virus$Topology), cex.lab = 1.5,
     xlab = '[log10] Genome Length (nt)', ylab = '[log10] Mean overlap lenght (nt)', main = "Overlap lenght")

cex <- 1

points( x[virus$baltimore.class%in%'ds_DNA'], y[virus$baltimore.class%in%'ds_DNA'], col=as.character(pal$col[pal$bal.class=="ds_DNA"]), pch=17, cex = cex)
points( x[virus$baltimore.class%in%'ss_RNA_+'], y[virus$baltimore.class%in%'ss_RNA_+'], col=as.character(pal$col[pal$bal.class=="ss_RNA_+"]), pch=4, cex = cex, lwd = lwd)
points( x[virus$baltimore.class%in%'ss_DNA'], y[virus$baltimore.class%in%'ss_DNA'], col=as.character(pal$col[pal$bal.class=="ss_DNA"]), pch=18, cex = cex)
points( x[virus$baltimore.class%in%'unknown_RNA'], y[virus$baltimore.class%in%'unknown_RNA'], col=as.character(pal$col[pal$bal.class=="unknown_RNA"]), pch=5, cex = cex, lwd = lwd)
points( x[virus$baltimore.class%in%'ds_RNA'], y[virus$baltimore.class%in%'ds_RNA'], col=as.character(pal$col[pal$bal.class=="ds_RNA"]), pch=8, cex = cex, lwd = 2)
points( x[virus$baltimore.class%in%'ss_RNA_-'], y[virus$baltimore.class%in%'ss_RNA_-'], col=as.character(pal$col[pal$bal.class=="ss_RNA_-"]), pch=6, cex = cex)
points( x[virus$baltimore.class%in%'unknown_DNA'], y[virus$baltimore.class%in%'unknown_DNA'], col=as.character(pal$col[pal$bal.class=="unknown_DNA"]), pch=20, cex = cex)
points( x[virus$baltimore.class%in%'RT_viruses'], y[virus$baltimore.class%in%'RT_viruses'], col=as.character(pal$col[pal$bal.class=="RT_viruses"]), pch=1, cex = cex, lwd = lwd)
points( x[virus$baltimore.class%in%'circular_ss_RNA'], y[virus$baltimore.class%in%'circular_ss_RNA'], col=as.character(pal$col[pal$bal.class=='circular_ss_RNA']), pch=4, cex = cex)

#source("http://www.math.mcmaster.ca/bolker/R/misc/legendx.R") # To adjust the size of boxes in the legend
legend(5000000, 1000, legend = pal$bal.class[1:8], fill = as.character(pal$col), box.cex=c(1,1), ncol = 1, cex = 1.2, border= F, bty="n", y.intersp = 1.5)
dev.off()

#######################################################
# Ridgeplot average overlap length
#######################################################
library(ggplot2)
library(ggridges)
library(gridExtra)
virus <- virus[which(virus$baltimore.class!= "circular_ss_RNA"),]

# Ordering baltimore groups
level_order <- factor(virus$baltimore.class, level = c("ds_RNA", "ss_RNA_+", "ss_RNA_-", "unknown_RNA", "ds_DNA", "ss_DNA", "unknown_DNA", "RT_viruses", "circular_ss_RNA"))
level_order <- factor(virus$baltimore.class, level = c("circular_ss_RNA","RT_viruses", "ds_DNA", "ss_DNA", "unknown_DNA", "ds_RNA", "ss_RNA_+", "ss_RNA_-", "unknown_RNA" ))

pdf(file='ridge_plot.pdf', width=8, height=11)
# By length
ov.len<- ggplot(virus, aes(x = log10(len.overlaps), y = level_order, group = baltimore.class, fill = baltimore.class)) + 
  geom_density_ridges(alpha = 1)+
  ggtitle("Average overlap length per Baltimore classification")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(y="Category", x = "[log10] Mean overlap length (nts)") +
  scale_fill_manual(values = c("circular_ss_RNA" = "#73af48", "ds_RNA" = "#1d6996", "ss_RNA_+" = "#38a6a5", "ss_RNA_-" = "#57a1b4", "unknown_RNA" = "#87bfdb",
                               "ds_DNA" = "#edad08",  "ss_DNA" = "#e17c05", "unknown_DNA" = "#cc503e", "RT_viruses" = "#808080", "circular_ss_RNA" = "#73af48" )) +
  #scale_color_manual(values = as.character(pal$col))+ 
  theme_bw()

ov.len
dev.off()

# organizing colors
for (i in 1:nrow(pal)){
  molecule <- as.character(pal$balt.cat[i])
  color <- as.character(pal$col[i])
  tog <- paste(molecule, color, sep=" = ")
  string <- append(string, tog)
}


#######################################################
# Frame shift ridgeplot
#######################################################
setwd('~/Projects/ovrf-review/data/')
overlaps <- read.csv('ovrfs-reduced.csv')
orfs <- read.csv('orfs-fixed.csv')

# add levels without overlaps
overlaps$accn <- factor(
  overlaps$accn, levels=c(
    levels(overlaps$accn),
    setdiff(levels(orfs$accno), levels(overlaps$accn))
  ))
temp <- sapply(split(overlaps$overlap, overlaps$accn), length)
noverlaps <- data.frame(accn=names(temp), count=temp)
noverlaps$mean.olen <- sapply(split(overlaps$overlap, overlaps$accn), mean)

# Ridgeplot
plot <- ggplot(overlaps, aes(x = log10(overlap), y = shift, group = shift, fill = as.factor(shift))) +
  geom_density_ridges()+
  ggtitle("Overlap length per frame shift")+
  scale_fill_manual(values = c("#eccb77", "#c6e377", "#729d39", "#36622b", "#1e3f2b")) +
  theme_bw()
plot
 

#######################################################
# Maptree for frameshifts
#######################################################
library(treemap)
library(plyr)
require(grid)
colnames(overlaps)[1] <- "Accession"
total <- merge(x=overlaps, y=virus[,c("Accession", "baltimore.class")], by="Accession", all.x =T)
tbl <- as.data.frame(table(total$shift, total$baltimore.class))
colnames(tbl) <- c("shift", "balt", "Freq")
new_p<-c("#9d2503", "#df8543", "#f6e0a7", "#90a95c","#5a715c")

pdf(file="treemap2.pdf")

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2)))
sub <- subset(tbl, balt == "ds_DNA")
treemap(sub, index = "shift", vSize = "Freq", type="index", 
        palette=new_p, title="ds_DNA", fontsize.title=15, fontsize.labels=14,
        vp = vplayout(1,1))

sub <- subset(tbl, balt == "ds_RNA")
treemap(sub, index = "shift", vSize = "Freq", type="index", 
        palette=new_p, title="ds_RNA", fontsize.title=15, fontsize.labels=14,
        vp = vplayout(1,2))

sub <- subset(tbl, balt == "ss_DNA")
treemap(sub, index = "shift", vSize = "Freq", type="index", 
        palette=new_p, title="ss_DNA", fontsize.title=15, fontsize.labels=14,
        vp = vplayout(2,1))

sub <- subset(tbl, balt == "ss_RNA_-")
treemap(sub, index = "shift", vSize = "Freq", type="index", 
        palette=new_p, title="ss_RNA_-", fontsize.title=15, fontsize.labels=14,
        vp = vplayout(2,2))

sub <- subset(tbl, balt == "ss_RNA_+")
treemap(sub, index = "shift", vSize = "Freq", type="index", 
        palette=new_p, title="ss_RNA_+", fontsize.title=15, fontsize.labels=14,
        vp = vplayout(3,1))

sub <- subset(tbl, balt == "RT_viruses")
treemap(sub, index = "shift", vSize = "Freq", type="index", 
        palette=new_p, title="RT_viruses", fontsize.title=15, fontsize.labels=14,
        vp = vplayout(3,2))

dev.off()

#################################################
# Further analysis for long overlaps
################################################
long<-subset(overlaps, overlap>=50)


