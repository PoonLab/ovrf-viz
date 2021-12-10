library(RColorBrewer)
library(ggplot2)
options(scipen=999)
library(plyr)
library(dplyr)

######################################################
# Preparing information and color palette 
######################################################
setwd('/home/lmunoz/Projects/ovrf-review/dataset_2020')
virus <- read.csv('no_cero_out_virus_df_inspection.csv', stringsAsFactors = F)
#virus <- read.csv('50_NOISE_no_cero_out_virus_df_inspection.csv', stringsAsFactors = F)

# DO NOT USE, PLOT IS SKEWED.
# If NA, the virus do not have overlapping proteins, therefore == 0 
virus[which(is.na(virus$n.overlaps)),]$n.overlaps <- 0
#virus[which(is.na(virus$len.overlaps)),]$len.overlaps <- 0
#virus[which(is.na(virus$total.overlaps)),]$total.overlaps <- 0

#Re naming Unknown according to molecule type
virus$baltimore.class[which(virus$baltimore.class == "Unknown" & grepl("RNA", virus$Molecule))] <- "unknown_RNA"
virus$baltimore.class[which(virus$baltimore.class == "Unknown" & grepl("DNA", virus$Molecule))] <- "unknown_DNA"

# Setting color palette
bal.class <- c("ds_RNA", "ss_RNA_+", "ss_RNA_-", "unknown_RNA", "ds_DNA", "ss_DNA", "unknown_DNA", "RT_viruses", "circular_ss_RNA")
col<- c("#1d6996", "#38a6a5", "#57a1b4", "#87bfdb", "#edad08", "#e17c05", "#cc503e", "#808080", "#73af48")
shape <- c(8, 4, 6, 5, 17, 18, 20, 1, 4)
pal <- as.data.frame(cbind(bal.class, col))
sh <- cbind.data.frame(bal.class, shape)

#######################################################
# Barplot database distribution
#######################################################
tbl <- as.data.frame(table(virus$baltimore.class))
tbl$bal.class <- factor(tbl$Var1, levels=c("ds_DNA", "ss_DNA", "unknown_DNA", "ds_RNA", "ss_RNA_+", "ss_RNA_-", "unknown_RNA", "RT_viruses" ))
# Vertical barplot
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

# Horizontal, paper barplot
#vector of colors based on baltimore class
bal.col <- merge(tbl[,c("bal.class", "Freq")], pal, by="bal.class") # freq of baltimore classification with assigned colors 
my_order <- rev(c("ds_DNA", "ss_DNA", "unknown_DNA", "ds_RNA", "ss_RNA_+", "ss_RNA_-", "unknown_RNA", "RT_viruses" ))
bal.col <- bal.col[match(my_order,bal.col$bal.class),]

##pdf(file="/home/lmunoz/Projects/ovrf-review/paper_plots/database_distribution.jpg")
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
#dev.off()

# Proteins or overlaps by baltimore class
prot<- as.data.frame(tapply(virus$Proteins, virus$baltimore.class, FUN=sum, na.rm=T))  
over<-tapply(virus$n.overlaps, virus$baltimore.class, FUN=sum, na.rm=T)
merge(prot, over, by=baltimore.class)
count<- ddply(virus, "baltimore.class", summarise,
              overlaps=sum(n.overlaps, na.rm=TRUE),
              proteins=sum(Proteins, na.rm= TRUE))
all<-melt(count)

ggplot(all, aes(x=baltimore.class, y=value, fill=factor(variable)))+
  geom_bar(stat="identity",position="dodge")

#######################################################
# Plotting overlapping summary
#######################################################

#pdf(file="/home/lmunoz/Projects/ovrf-review/dataset_2020/general_trends_no_cero.pdf", width=14, height=6.5)
# a. Number of ORFS vs Number of overlaps
set.seed(5)
x <- virus$Proteins
y <- virus$n.overlaps
y[y==0] <- 0.5
x <- jitter(x)
y<-ifelse(y==0.5, jitter(y), y)

#layout(matrix(c(1,2), nrow=1, ncol=2, byrow = TRUE))
par(mfrow=c(1,2))
par(mar = c(5, 5, 4, 1))

plot(x,y, log='xy', type='n', col=as.numeric(virus$Molecule.type), cex.main=1.9, pch=as.numeric(virus$Topology), 
     xlab = 'Number of ORFs', ylab = 'Number of overlaps', main = "Virus data Base", cex.lab=1.6, cex.axis=1.5)
lwd <- 2
cex <- 1
cex_rna<-0.8
points( x[virus$baltimore.class%in%'ds_DNA'], y[virus$baltimore.class%in%'ds_DNA'], col=as.character(pal$col[pal$bal.class=="ds_DNA"]), pch=17, cex = cex, lwd = lwd)
points( x[virus$baltimore.class%in%'ss_RNA_+'], y[virus$baltimore.class%in%'ss_RNA_+'], col=as.character(pal$col[pal$bal.class=="ss_RNA_+"]), pch=4, cex = cex_rna, lwd = lwd)
points( x[virus$baltimore.class%in%'ss_DNA'], y[virus$baltimore.class%in%'ss_DNA'], col=as.character(pal$col[pal$bal.class=="ss_DNA"]), pch=18, cex = cex, lwd = lwd)
points( x[virus$baltimore.class%in%'unknown_RNA'], y[virus$baltimore.class%in%'unknown_RNA'], col=as.character(pal$col[pal$bal.class=="unknown_RNA"]), pch=5, cex = cex_rna, lwd = lwd)
points( x[virus$baltimore.class%in%'ss_RNA_-'], y[virus$baltimore.class%in%'ss_RNA_-'], col=as.character(pal$col[pal$bal.class=="ss_RNA_-"]), pch=6, cex = cex_rna, lwd=lwd)
points( x[virus$baltimore.class%in%'ds_RNA'], y[virus$baltimore.class%in%'ds_RNA'], col=as.character(pal$col[pal$bal.class=="ds_RNA"]), pch=8, cex = cex_rna, lwd = 2)
points( x[virus$baltimore.class%in%'unknown_DNA'], y[virus$baltimore.class%in%'unknown_DNA'], col=as.character(pal$col[pal$bal.class=="unknown_DNA"]), pch=20, cex = cex, lwd = lwd)
points( x[virus$baltimore.class%in%'RT_viruses'], y[virus$baltimore.class%in%'RT_viruses'], col=as.character(pal$col[pal$bal.class=="RT_viruses"]), pch=1, cex = cex, lwd = lwd)
#points( x[virus$baltimore.class%in%'circular_ss_RNA'], y[virus$baltimore.class%in%'circular_ss_RNA'], col=as.character(pal$col[pal$bal.class=='circular_ss_RNA']), pch=4, cex = cex)
#abline(lm(y~x), col="red", cex=8) # regression line (y~x) 

# b. Genome length vs mean nucleotides involved in an overlap
y <- virus$Length
x <- virus$len.overlaps
y[y==0] <- 0.5
par(xpd = T, mar = c(5, 5, 4, 1))
plot(x,y, type='n', log = 'xy', col=as.numeric(virus$Molecule.type), cex.main=1.9, pch=as.numeric(virus$Topology),
     xlab = 'Overlap length (nt)', ylab = 'Genome Length (nt)', main = "Overlap length", cex.lab = 1.6, cex.axis=1.5)

points( x[virus$baltimore.class%in%'ds_DNA'], y[virus$baltimore.class%in%'ds_DNA'], col=as.character(pal$col[pal$bal.class=="ds_DNA"]), pch=17, cex = cex, lwd = lwd)
points( x[virus$baltimore.class%in%'ss_RNA_+'], y[virus$baltimore.class%in%'ss_RNA_+'], col=as.character(pal$col[pal$bal.class=="ss_RNA_+"]), pch=4, cex = cex_rna, lwd = lwd)
points( x[virus$baltimore.class%in%'ss_DNA'], y[virus$baltimore.class%in%'ss_DNA'], col=as.character(pal$col[pal$bal.class=="ss_DNA"]), pch=18, cex = cex, lwd = lwd)
points( x[virus$baltimore.class%in%'unknown_RNA'], y[virus$baltimore.class%in%'unknown_RNA'], col=as.character(pal$col[pal$bal.class=="unknown_RNA"]), pch=5, cex = cex_rna, lwd = lwd)
points( x[virus$baltimore.class%in%'ss_RNA_-'], y[virus$baltimore.class%in%'ss_RNA_-'], col=as.character(pal$col[pal$bal.class=="ss_RNA_-"]), pch=6, cex = cex_rna, lwd = lwd)
points( x[virus$baltimore.class%in%'ds_RNA'], y[virus$baltimore.class%in%'ds_RNA'], col=as.character(pal$col[pal$bal.class=="ds_RNA"]), pch=8, cex = cex_rna, lwd = 2)
points( x[virus$baltimore.class%in%'unknown_DNA'], y[virus$baltimore.class%in%'unknown_DNA'], col=as.character(pal$col[pal$bal.class=="unknown_DNA"]), pch=20, cex = cex)
points( x[virus$baltimore.class%in%'RT_viruses'], y[virus$baltimore.class%in%'RT_viruses'], col=as.character(pal$col[pal$bal.class=="RT_viruses"]), pch=1, cex = cex, lwd = lwd)
points( x[virus$baltimore.class%in%'circular_ss_RNA'], y[virus$baltimore.class%in%'circular_ss_RNA'], col=as.character(pal$col[pal$bal.class=='circular_ss_RNA']), pch=4, cex = cex)

source("http://www.math.mcmaster.ca/bolker/R/misc/legendx.R") # To adjust the size of boxes in the legend
legend(5000000, 1000000, legend = pal$bal.class[1:8], col = as.character(pal$col), pch = sh$shape, ncol = 1, cex = 2, border= F, bty="n", y.intersp = 1.5, lwd=3)
#dev.off()

#######################################################
# Ridgeplot average overlap length
#######################################################
i
library(ggplot2)
library(ggridges)
library(gridExtra)
virus <- virus[which(virus$baltimore.class!= "circular_ss_RNA"),]

# Ordering baltimore groups
level_order <- factor(virus$baltimore.class, level = c("ds_RNA", "ss_RNA_+", "ss_RNA_-", "unknown_RNA", "ds_DNA", "ss_DNA", "unknown_DNA", "RT_viruses", "circular_ss_RNA"))
level_order <- factor(virus$baltimore.class, level = c("circular_ss_RNA","RT_viruses", "ds_DNA", "ss_DNA", "unknown_DNA", "ds_RNA", "ss_RNA_+", "ss_RNA_-", "unknown_RNA" ))

#pdf(file='ridge_plot.#pdf', width=8, height=11)
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
#dev.off()

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
setwd('/home/lmunoz/Projects/ovrf-review/dataset_2020/')
overlaps1 <- read.csv('ovrfs_regular_ex.csv')
overlaps <- read.csv('2_overlapping_march2020.csv', colClasses='character')
overlaps$loc1 <- as.integer(overlaps$loc1)
overlaps$dir1 <- as.factor(overlaps$dir1)
overlaps$loc2 <- as.integer(overlaps$loc2)
overlaps$dir2 <- as.factor(overlaps$dir2)
overlaps$seqlen1 <- as.integer(overlaps$seqlen1)
overlaps$seqlen2 <- as.integer(overlaps$seqlen2)
overlaps$overlap <- as.integer(overlaps$overlap)
overlaps$shift <- factor(overlaps$shift, 
                         levels=c('-2', '-1', '-0', '+0', '+1', '+2'))

orfs <- read.csv('total_orfs.csv')

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
plot <- ggplot(overlaps, aes(x = log10(overlap), y = as.character(shift), group = as.character(shift), fill = as.factor(shift))) +
  geom_density_ridges()+
  ggtitle("Overlap length per frame shift")+
  #scale_fill_manual(values = c("#eccb77", "#c6e377", "#729d39", "#36622b", "#1e3f2b")) +
  theme_bw()
plot