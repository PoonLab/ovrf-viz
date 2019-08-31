setwd('~/git/ovrf-review/data')
virus <- read.csv('viruses.csv')

virus$Genome.length <- as.integer(gsub(' nt$', '', virus$Genome.length))
summary(virus)

virus$Number.of.proteins[virus$Number.of.proteins=='-'] <- NA
virus$Number.of.proteins <- as.integer(as.character(virus$Number.of.proteins))
summary(virus)

table(virus$Family)
boxplot(split(log10(virus$Genome.length), virus$Family))
boxplot(split(log10(virus$Number.of.proteins), virus$Family))


orfs <- read.csv('orfs-fixed.csv')


overlaps <- read.csv('find_ovrfs.csv')


# 1. how many overlaps per genome?
levels(overlaps$accn) <- unique(orfs$accno)
temp <- sapply(split(overlaps$overlap, overlaps$accn), length)
noverlaps <- data.frame(accn=names(temp), count=temp)

foo <- sapply(noverlaps$accn[1:10], function(accn) {
  which(grepl(as.character(accn), virus$Accession))
})

noverlaps$genome <- 