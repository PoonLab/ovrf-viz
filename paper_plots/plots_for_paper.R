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
# Barplot
#######################################################
table <- as.data.frame(table(virus$baltimore.class))
table$Var1 <- factor(table$Var1, levels=c("ds_DNA", "ss_DNA", "unknown_DNA", "ds_RNA", "ss_RNA_+", "ss_RNA_-", "unknown_RNA", "RT_viruses" ))
p <- ggplot(table, aes(x=Var1, y=Freq, fill=Var1)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("circular_ss_RNA" = "#73af48", "ds_RNA" = "#1d6996", "ss_RNA_+" = "#38a6a5", "ss_RNA_-" = "#57a1b4", "unknown_RNA" = "#87bfdb",
                               "ds_DNA" = "#edad08",  "ss_DNA" = "#e17c05", "unknown_DNA" = "#cc503e", "RT_viruses" = "#808080", "circular_ss_RNA" = "#73af48" ))+
  
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25)+
  labs(fill = "Baltimore class")+ 
  theme_bw()+
  ggtitle("Database distribution")

p
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
# Frameshift barplot with Baltimore class
#######################################################
colnames(overlaps)[1] <- "Accession"

total <- merge(x=overlaps, y=virus[,c("Accession", "baltimore.class")], by="Accession", all.x =T)
table <- table(total$baltimore.class, total$shift)
plot <- ggplot(table, aes(x=, y=baltimore.class, group = baltimore.class, fill = baltimore.class))


#result <- within(total, {count<-ave(baltimore.class, shift, FUN=function(x) length(x))})
library(dplyr)

