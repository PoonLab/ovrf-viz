# File with Baltimore info
setwd('~/Projects/ovrf-review/data/Baltimore')
virus <- read.csv('species_file_baltimore2.csv')


#stacked bar plot
more_than_zero <- subset(virus, n.overlaps > 0)
pal <- c("#f1d4d4", "#ddb6c6", "#ac8daf", "#484c7f")
counts <- table(more_than_zero$family, more_than_zero$baltimore.class)
barplot(counts, col = pal, space =0.8)

#Stores the Column names in l, so we can iterate column by column
l <- dimnames(counts)[[2]]

#Applies a function to every member of some list 
#(ie. Takes inputs from l, and does the same thing, with each col heading)
blah <- lapply(l, function(i){
  #Get counts given the header
  x <- counts[,i]
  #Get the largest three counts
  pos <- tail(order(x),3)
  #Get the family names of those three counts
  family <- names(x[pos])
  return(family)
})

specific_fam <- blah[[2]][3]
counts[specific_fam,]




# Ridgeplot
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
