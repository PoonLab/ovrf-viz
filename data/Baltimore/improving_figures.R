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

# TODO: How do I relate the names of the most abundant families with the actual barplot?


# Ridgeplot
require(devtools)
devtools::install_github("ArtPoon/ggfree")


subset1 <- subset(virus, n.overlaps < 1)
subset2 <- subset(virus, 1 <= n.overlaps & n.overlaps < 10)
subset3 <- subset(virus, 10 <= n.overlaps & n.overlaps < 20)
subset4 <- subset(virus, 20 <= n.overlaps & n.overlaps < 30)
subset5 <- subset(virus, 30 <= n.overlaps & n.overlaps < 40)
subset6 <- subset(virus, 40 <= n.overlaps & n.overlaps < 50)
subset7 <- subset(virus, 50<n.overlaps & n.overlaps<=100)
subset8 <- subset(virus, 100<n.overlaps & n.overlaps<=500)
subset9 <- subset(virus, 500<n.overlaps)

hist(subset2$n.overlaps)

ridgeplot <- function(x, xlim=NA, labels=NA, yaxt='s', xlab=NA, ylab=NA, step=0.2,
                      col=NA, fill=NA, lwd=1, density.args=list(),
                      add.grid=F, grid.args=list(), extend.lines=TRUE, add=FALSE, ...) {
  
  # check inputs
  if (is.list(x)) {
    # check that list contains numeric vectors
    if (!all(sapply(x, is.numeric))) {
      stop("List 'x' must contain only numeric vectors")
    }
  } else {
    stop("Unsupported class of argument x (must be list)")
  }
  
  n <- length(x)  # number of groups
  
  if (any(is.na(labels))) {
    # extract labels from list 'x'
    labels <- names(x)
  }
  
  # parse optional colour vector arguments
  if (all(is.na(col))) {
    pal <- rep('black', n)  # default
  } else {
    # reuse values in <col> if less than n
    pal <- rep(col, length.out=n)
  }
  if (all(is.na(fill))) {
    bg <- rep(NA, n)
  } else {
    bg <- rep(fill, length.out=n)
  }
  
  # generate kernel densities
  kdens <- list()
  for (i in 1:length(x)) {
    xi <- x[[i]]
    if (length(xi) == 0 && is.numeric(xi)) {
      # is numeric(0)
      kdens[[names(x)[i]]] <- NA
    } else {
      kdens[[names(x)[i]]] <- do.call('density', c(list(x=xi), density.args))
    }
  }
  
  if (!add) {
    # determine plot ranges from densities
    all.x <- c(sapply(kdens, function(k) k$x))
    all.y <- c(sapply(1:n, function(i) kdens[[i]]$y + i*step))
    
    # generate plot region
    if (any(is.na(xlim))) {
      xlim <- range(all.x)
    }
    plot(NA, xlim=xlim, ylim=range(all.y), xlab=xlab, ylab=ylab, yaxt='n', ...)
    
    # override y-axis labels
    if (yaxt != 'n') axis(side=2, at=seq(step, n*step, step), labels=labels, las=2)
    
    # optional grid
    if (add.grid) {
      do.call("add.grid", grid.args)
    }
  }
  
  # arrangement of densities
  ordering <- seq(n, 1, -1)
  if (step < 0) {
    ordering <- 1:n
  }
  for (i in ordering) {
    kd <- kdens[[i]]
    if (length(kd) == 1 && is.na(kd)) next  # skip missing entry
    y <- kd$y + i*step
    
    polygon(kd$x, y, col=bg[i], border=NA)
    lines(kd$x, y, col=pal[i], lwd=lwd)
    
    # extend density curve to full horizontal range
    if (extend.lines) {
      segments(x0=min(all.x), x1=min(kd$x), y0=min(y), y1=min(y), col=pal[i], lwd=lwd)
      segments(x0=max(kd$x), x1=max(all.x), y0=min(y), y1=min(y), col=pal[i], lwd=lwd)
    }
  }
}


ridgeplot(subset2$n.overlaps)
