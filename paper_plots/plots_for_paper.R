library(RColorBrewer)
options(scipen=999)

######################################################
# Preparing information and color palette 
######################################################
setwd('/home/lmunoz/Projects/ovrf-review/data/Baltimore')
virus <- read.csv('species_file_baltimore2.csv', stringsAsFactors = F)

#Re naming Unknown according to molecule type
virus$baltimore.class[which(virus$baltimore.class == "Unknown" & grepl("RNA", virus$Molecule.type))] <- "unknown_RNA"
virus$baltimore.class[which(virus$baltimore.class == "Unknown" & grepl("DNA", virus$Molecule.type))] <- "unknown_DNA"

# Setting color palette
balt.cat <- c("ds_RNA", "ss_RNA_+", "ss_RNA_-", "unknown_RNA", "ds_DNA", "ss_DNA", "unknown_DNA", "RT_viruses", "circular_ss_RNA")
col<- c("#1d6996", "#38a6a5", "#57a1b4", "#87bfdb", "#edad08", "#e17c05", "#cc503e", "#808080", "#73af48")
pal <- as.data.frame(cbind(balt.cat, col))


#######################################################
# Plotting overlapping summary
#######################################################
# a. Number of ORFS vs Number of overlaps
x <- virus$Number.of.proteins
y <- virus$n.overlaps
y[y==0] <- 0.5
y <- jitter(y)
x <- jitter(x)
#layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
par(mfrow=c(1,2))
plot(x,y, log='xy', type='n', col=as.numeric(virus$Molecule.type), cex=1, pch=as.numeric(virus$Topology), 
     xlab = '[log10] Number of ORFs', ylab = '[log10] Number of overlaps', main = "Virus data Base")

cex <- 0.8
points( x[virus$baltimore.class%in%'ds_DNA'], y[virus$baltimore.class%in%'ds_DNA'], col=as.character(pal$col[pal$balt.cat=="ds_DNA"]), pch=17, cex = cex)
points( x[virus$baltimore.class%in%'ss_RNA_+'], y[virus$baltimore.class%in%'ss_RNA_+'], col=as.character(pal$col[pal$balt.cat=="ss_RNA_+"]), pch=4, cex = cex)
points( x[virus$baltimore.class%in%'ss_DNA'], y[virus$baltimore.class%in%'ss_DNA'], col=as.character(pal$col[pal$balt.cat=="ss_DNA"]), pch=18, cex = cex)
points( x[virus$baltimore.class%in%'unknown_RNA'], y[virus$baltimore.class%in%'unknown_RNA'], col=as.character(pal$col[pal$balt.cat=="unknown_RNA"]), pch=5, cex = cex)
points( x[virus$baltimore.class%in%'ds_RNA'], y[virus$baltimore.class%in%'ds_RNA'], col=as.character(pal$col[pal$balt.cat=="ds_RNA"]), pch=8, cex = cex)
points( x[virus$baltimore.class%in%'ss_RNA_-'], y[virus$baltimore.class%in%'ss_RNA_-'], col=as.character(pal$col[pal$balt.cat=="ss_RNA_-"]), pch=6, cex = cex)
points( x[virus$baltimore.class%in%'unknown_DNA'], y[virus$baltimore.class%in%'unknown_DNA'], col=as.character(pal$col[pal$balt.cat=="unknown_DNA"]), pch=20, cex = cex)
points( x[virus$baltimore.class%in%'RT_viruses'], y[virus$baltimore.class%in%'RT_viruses'], col=as.character(pal$col[pal$balt.cat=="RT_viruses"]), pch=1, cex = cex)
points( x[virus$baltimore.class%in%'circular_ss_RNA'], y[virus$baltimore.class%in%'circular_ss_RNA'], col=as.character(pal$col[pal$balt.cat=='circular_ss_RNA']), pch=4, cex = cex)


legend("topleft", legend = pal$balt.cat, fill = as.character(pal$col), ncol = 2, cex = 0.7)

# b. Genome length vs mean nucleotides involved in an overlap
x <- virus$Genome.length
y <- virus$len.overlaps
y <- jitter(y)
x <- jitter(x)
plot(x,y, type='n', log = 'xy', col=as.numeric(more_than_one$Molecule.type), cex=2, pch=as.numeric(more_than_one$Topology), 
     xlab = '[log10] Genome Length (nt)', ylab = '[log10] Mean overlap lenght (nt)', main = "Overlap lenght")

cex <- 0.8
points( x[virus$baltimore.class%in%'ds_DNA'], y[virus$baltimore.class%in%'ds_DNA'], col=as.character(pal$col[pal$balt.cat=="ds_DNA"]), pch=17, cex = cex)
points( x[virus$baltimore.class%in%'ss_RNA_+'], y[virus$baltimore.class%in%'ss_RNA_+'], col=as.character(pal$col[pal$balt.cat=="ss_RNA_+"]), pch=4, cex = cex)
points( x[virus$baltimore.class%in%'ss_DNA'], y[virus$baltimore.class%in%'ss_DNA'], col=as.character(pal$col[pal$balt.cat=="ss_DNA"]), pch=18, cex = cex)
points( x[virus$baltimore.class%in%'unknown_RNA'], y[virus$baltimore.class%in%'unknown_RNA'], col=as.character(pal$col[pal$balt.cat=="unknown_RNA"]), pch=5, cex = cex)
points( x[virus$baltimore.class%in%'ds_RNA'], y[virus$baltimore.class%in%'ds_RNA'], col=as.character(pal$col[pal$balt.cat=="ds_RNA"]), pch=8, cex = cex)
points( x[virus$baltimore.class%in%'ss_RNA_-'], y[virus$baltimore.class%in%'ss_RNA_-'], col=as.character(pal$col[pal$balt.cat=="ss_RNA_-"]), pch=6, cex = cex)
points( x[virus$baltimore.class%in%'unknown_DNA'], y[virus$baltimore.class%in%'unknown_DNA'], col=as.character(pal$col[pal$balt.cat=="unknown_DNA"]), pch=20, cex = cex)
points( x[virus$baltimore.class%in%'RT_viruses'], y[virus$baltimore.class%in%'RT_viruses'], col=as.character(pal$col[pal$balt.cat=="RT_viruses"]), pch=1, cex = cex)
points( x[virus$baltimore.class%in%'circular_ss_RNA'], y[virus$baltimore.class%in%'circular_ss_RNA'], col=as.character(pal$col[pal$balt.cat=='circular_ss_RNA']), pch=4, cex = cex)



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


# organizing colors
for (i in 1:nrow(pal)){
  molecule <- as.character(pal$balt.cat[i])
  color <- as.character(pal$col[i])
  tog <- paste(molecule, color, sep=" = ")
  string <- append(string, tog)
}


#######################################################
# Frame shift analysis
#######################################################
setwd("/home/lmunoz/Projects/ovrf-review/scripts")
orfs <- read.csv("orfs_with_nucleotide.csv")
colnames(orfs)[1]<- "Accession"
total <- merge(x = orfs[,c("Accession", "strand", "coords")], y = virus[,c("Accession", "baltimore.class")], by="Accession")

strand <- ggplot(orfs, aes())
