# File with Baltimore info
setwd('~/Projects/ovrf-review/data/Baltimore')
virus <- read.csv('species_file_baltimore2.csv')


#stacked bar plot
more_than_zero <- subset(virus, len.overlaps > 0)
pal <- c("#f1d4d4", "#ddb6c6", "#ac8daf", "#484c7f")
counts <- table(more_than_zero$family, more_than_zero$baltimore.class)
matrix <- as.matrix(counts)
barplot(matrix, col = pal, legend = rownames(matrix))

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

# TODO: How do I relate the names of the most abundant families with the actual barplot?




########################################
# Ridgeplot for number of overlaps
########################################
#require(devtools)
devtools::install_github("ArtPoon/ggfree")

library(ggplot2)
library(ggridges)

test <- virus
test$subset <- NA

test$subset[test$n.overlaps < 1] <- "Zero overlaps"
test$subset[1 <= test$n.overlaps & test$n.overlaps < 10] <- "1 to 10"
test$subset[10 <= test$n.overlaps & test$n.overlaps < 20] <- "10 to 20"
test$subset[20 <= test$n.overlaps & test$n.overlaps < 30] <- "20 to 30"
test$subset[30 <= test$n.overlaps & test$n.overlaps < 40] <- "30 to 40"
test$subset[40 <= test$n.overlaps & test$n.overlaps < 50] <- "40 to 50"
test$subset[50 <= test$n.overlaps & test$n.overlaps < 100] <- "50 to 100"
test$subset[100 <= test$n.overlaps & test$n.overlaps < 500] <- "100 to 500"
test$subset[500 < test$n.overlaps] <- "> 500"

test <- subset(test, !is.na(n.overlaps))

ggplot(test, aes(x = n.overlaps, y = subset, group = subset)) + 
  geom_density_ridges(
    jittered_points = TRUE, quantile_lines = TRUE, scale = 0.9, alpha = 0.7,
    vline_size = 1, vline_color = "red", vline_alpha = 1,
    point_size = 0.4, point_alpha = 1,
    position = position_raincloud(adjust_vlines = TRUE)
  )


ggplot(test, aes(x = n.overlaps, y = subset, group = subset)) + 
  geom_density_ridges(
    aes(point_color = subset, point_fill = subset, point_shape = subset),
    alpha = .2, point_alpha = 1, jittered_points = TRUE
  ) +
  scale_point_color_hue(l = 40) +
  scale_discrete_manual(aesthetics = "point_shape", values = c(18,17,19, 20, 21, 22, 23, 25, 25))

#Small overlaps
test1 <- subset(test, n.overlaps < 100)

splitted<-split(test$n.overlaps, test$subset)

ggplot(test1, aes(x = n.overlaps, y = subset, group = subset)) + 
  geom_density_ridges(scale = 3)


# Large overlaps
test2 <- subset(test, n.overlaps > 500)


ggplot(test2, aes(x = n.overlaps, y = subset, group = subset)) + 
  geom_density_ridges(scale = 9)


########################################
# Ridgeplot for Medium overlap length
########################################
more_than_zero$ov.len.class[more_than_zero$len.overlaps < 1] <- "Zero overlaps"
more_than_zero$ov.len.class[1 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 10] <- "1 to 10"
more_than_zero$ov.len.class[10 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 20] <- "10 to 20"
more_than_zero$ov.len.class[20 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 30] <- "20 to 30"
more_than_zero$ov.len.class[30 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 40] <- "30 to 40"
more_than_zero$ov.len.class[40 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 50] <- "40 to 50"
more_than_zero$ov.len.class[50 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 100] <- "50 to 100"
more_than_zero$ov.len.class[100 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 200] <- "100 to 200"
more_than_zero$ov.len.class[200 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 300] <- "200 to 300"
more_than_zero$ov.len.class[300 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 400] <- "300 to 400"
more_than_zero$ov.len.class[400 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 500] <- "400 to 500"
more_than_zero$ov.len.class[500 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 600] <- "500 to 600"
more_than_zero$ov.len.class[600 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 700] <- "600 to 700"
more_than_zero$ov.len.class[700 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 800] <- "700 to 800"
more_than_zero$ov.len.class[800 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 900] <- "800 to 900"
more_than_zero$ov.len.class[900 <= more_than_zero$len.overlaps & more_than_zero$len.overlaps < 1000] <- "900 to 1000"
more_than_zero$ov.len.class[1000 < more_than_zero$len.overlaps] <- "> 1000"



ggplot(more_than_zero, aes(x = len.overlaps, y = ov.len.class, group = ov.len.class)) + 
  geom_density_ridges()


ggplot(more_than_zero, aes(x = len.overlaps, y = ov.len.class, group = ov.len.class)) + 
  geom_density_ridges(
    aes(point_color = ov.len.class, point_fill = ov.len.class, point_shape = ov.len.class),
    alpha = .2, point_alpha = 1, jittered_points = TRUE
  ) +
  scale_point_color_hue(l = 40) +
  scale_discrete_manual(aesthetics = "point_shape", values = c(18,17,19, 20, 21, 22, 23, 25, 25))

small <- subset(more_than_zero, len.overlaps < 100 & len.overlaps != 0)

ggplot(small, aes(x = len.overlaps, y = ov.len.class, group = ov.len.class)) + 
  geom_density_ridges(
    aes(point_color = ov.len.class, point_fill = ov.len.class, point_shape = ov.len.class),
    alpha = .2, point_alpha = 1, jittered_points = TRUE
  ) +
  scale_point_color_hue(l = 40) +
  scale_discrete_manual(aesthetics = "point_shape", values = c(18,17,19, 20, 21, 22, 23, 25, 25))

large <- subset(more_than_zero, len.overlaps > 100 & len.overlaps <1000)
ggplot(large, aes(x = len.overlaps, y = ov.len.class, group = ov.len.class)) + 
  geom_density_ridges(
    aes(point_color = ov.len.class, point_fill = ov.len.class, point_shape = ov.len.class),
    alpha = .2, point_alpha = 1, jittered_points = TRUE
  ) +
  scale_point_color_hue(l = 40) +
  scale_discrete_manual(aesthetics = "point_shape", values = c(18,17,19, 20, 21, 22, 23, 25, 25, 16))

