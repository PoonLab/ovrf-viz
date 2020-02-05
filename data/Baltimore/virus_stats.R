setwd('~/Projects/ovrf-review/data/Baltimore')
virus <- read.csv('species_file_baltimore2.csv')

library(ggplot2)
library(gridExtra)
library(dplyr)

# OVERLAP LENGTH 
tapply(virus$len.overlaps, virus$baltimore.class, summary)

# Boxplot
ggplot(virus, aes(baltimore.class, len.overlaps)) + geom_boxplot()

# Histograms
p1 <- ggplot(virus, aes(len.overlaps)) + 
  geom_histogram(fill = "white", color = "grey30") +
  facet_wrap(~ baltimore.class)

p2 <- ggplot(virus, aes(len.overlaps)) + 
  geom_histogram(fill = "white", color = "grey30") +
  facet_wrap(~ baltimore.class) +
  scale_x_log10()

grid.arrange(p1, p2, nrow=2)

rna<- subset(virus, xor(baltimore.class == "ds_DNA", baltimore.class == "ss_RNA_-"))

rna<- subset(virus, xor(Molecule.type == "DNA", Molecule.type == "RNA"))

wilcox.test(len.overlaps ~ Molecule.type, data = rna)


dna<- subset(virus, Molecule.type == "DNA")
rna_hist <- hist(log(rna$len.overlaps))

# ANOVA of overlap length between baltimore classes
sum <- aggregate(virus$len.overlaps, by=virus["baltimore.class"], FUN=summary, na.rm = T)
res.aov <- aov(log(len.overlaps) ~ baltimore.class, data = virus) 
summary(res.aov) # Is significant
TukeyHSD(res.aov)
plot(res.aov, 1) # Outliers possibly affecting normality and homogenity if variance

# Grapping dsDNA info
dsdna<- subset(virus, baltimore.class == "ds_DNA")
test <-grepl("Caudovirales", dsdna$Taxonomy) #Which are Bacteriophages
table(test)

#Long overlaps
subset(virus, len.overlaps > 5000)


# NUMBER OF OVERLAPS
tapply(virus$n.overlaps, virus$baltimore.class, summary)
tapply(virus$n.overlaps, virus$baltimore.class, sum, na.rm = T)

library(ggpubr)
ggplot(virus, aes(baltimore.class, log(n.overlaps))) + geom_boxplot()

p1 <- ggplot(virus, aes(n.overlaps)) + 
  geom_histogram(fill = "white", color = "grey30") +
  facet_wrap(~ baltimore.class)

p2 <- ggplot(virus, aes(n.overlaps)) + 
  geom_histogram(fill = "white", color = "grey30") +
  facet_wrap(~ baltimore.class) +
  scale_x_log10()

grid.arrange(p1, p2, nrow=2)

sum <- aggregate(virus$n.overlaps, by=virus["baltimore.class"], FUN=summary, na.rm = T)

kruskal.test(log(n.overlaps) ~ baltimore.class, data=virus) # There is a significant difference between groups

sum2 <- aggregate(virus$n.overlaps, by=virus["baltimore.class"], FUN=summary, na.rm = T)
pairwise.wilcox.test(virus$n.overlaps, virus$baltimore.class, p.adjust.method = "BH")


# Summary by group
group_by(virus, baltimore.class) %>%
  summarise(
    count = n(), 
    mean = mean(n.overlaps, na.rm= T),
    sd = sd(n.overlaps, na.rm = T),
    median = median(n.overlaps, na.rm = T),
    IQR=IQR(n.overlaps, na.rm = T)
  )
