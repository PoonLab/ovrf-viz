require(igraph)

setwd('~/git/ovrf-review/dataset_2020/plot_data/')

files <- Sys.glob('./*/nodelist.csv')
require(MASS)

parse.edgelist <- function(f) {
  #vir.name <- gsub("viridae", "", strsplit(f, "/")[[1]][2])
  nl <- read.csv(f)
  # read edge list
  el <- read.csv(gsub("nodelist", "edgelist", f))
  el$parent <- as.factor(el$parent)
  el$child <- as.factor(el$child)
  ael <- el[el$edge.type=='adjacent', ]
  oel <- el[el$edge.type=='overlap', ]
  
  # add overlap counts to adjacent counts
  ael$overlap <- apply(ael, 1, function(r) {
    ocount <- oel$edge.count[oel$parent == r[1] & oel$child == r[2]]
    ifelse (length(ocount) > 0, ocount, 0)
  })
  ael$edge.count <- ael$edge.count + ael$overlap
  
  # filter adjacent edges by minimum count
  ael <- ael[ael$edge.count >= 5, ]
  
  #max.size <- max(nl$size)
  ael$total.size <- nl$size[as.integer(ael$parent)] + nl$size[as.integer(ael$child)]
  
  ag <- graph_from_data_frame(ael)
  
  deg <- degree(ag)  # named vector
  ael$edge.degree <- deg[ael$parent] + deg[ael$child]
  # https://www.sciencedirect.com/science/article/pii/S0012365X07003251
  #ael$edge.degree <- apply(ael, 1, function(r) {
  #  length(union(
  #    neighbors(ag, as.integer(r[1])), 
  #    neighbors(ag, as.integer(r[2]))
  #    )) - 2
  #})
  
  tri <- count_triangles(ag)
  names(tri) <- names(deg)
  ael$triangles <- tri[ael$parent] + tri[ael$child]
  
  trans <- transitivity(ag, 'local')
  names(trans) <- names(deg)
  ael$transitivity <- trans[ael$parent] + trans[ael$child]
  
  eig <- eigen_centrality(ag)$vector
  ael$centrality <- eig[ael$parent] + eig[ael$child]
  
  # relabel nodes
  #ael$parent <- paste(vir.name, ael$parent, sep="")
  #ael$child <- paste(vir.name, ael$child, sep="")
  ael
}

#edgelists <- lapply(files[c(1, 3, 6, 9)], parse.edgelist)
#ael <- do.call("rbind", edgelists)

require(VGAM)

#res <- do.call("rbind", 
families <- c('Adenoviridae', 'Coronaviridae', 'Geminiviridae', 'Papillomaviridae', 'Rhabdoviridae')
res <- list()
count <- 1
for (i in c(1, 3, 4, 6, 9)) {
  ael <- parse.edgelist(files[i])
  fit <- vglm(cbind(overlap, edge.count) ~ total.size + edge.degree + 
                triangles + transitivity + centrality, zibinomialff, 
              data=ael, trace=F)
  fit0 <- vglm(cbind(overlap, edge.count) ~ 1, zibinomialff, data=ael, trace=F)
  print(i)
  fit.aic <- step4vglm(fit, scope=list(lower=fit0, upper=fit), direction='both', trace=F)
  #fit <- glm(cbind(overlap, edge.count) ~ total.size + edge.degree + triangles + 
  #             transitivity + centrality, data=ael, family='quasibinomial')
  
  ci <- as.data.frame(confint(fit.aic))
  ci$est <- fit.aic@coefficients
  ci$coef <- row.names(ci)
  ci$family <- rep(families[count], nrow(ci))
  res[[i]] <- ci
  count <- count + 1
}


#res <- do.call("rbind", res)
#row.names(res) <- NULL

require(ggfree)
pal <- hcl.colors(5, palette="Dark3")
pal2 <- hcl.colors(5, palette='Dark3', alpha=0.5)

pdf(file='~/papers/ovrf-viz/regression.pdf', width=5, height=5)
par(mfrow=c(5,1), mar=c(2,3,1,6), cex=1, xpd=F)
for (term in c('total.size', 'edge.degree', 'triangles', 'transitivity', 'centrality')) {
  temp <- t(sapply(c(1, 3, 4, 6, 9), function(i) {
    r <- res[[i]]
    if (!is.element(term, r$coef)) {
      rep(NA, 3)
    } else {
      as.numeric(r[which(r$coef==term), 1:3])
    }
  }))
  
  x <- temp[,3]
  lim.x <- max(abs(temp), na.rm=T)
  
  plot(x, 1:5, xlim=c(-lim.x, lim.x), bty='n', ylim=c(0.5, 5.5), yaxt='n', 
       ylab=NA, cex.axis=0.7, mgp=c(3,0.5,0))
  text(x=-1.2*lim.x, y=3, label=gsub('\\.', ' ', term), 
       xpd=NA, adj=0.5, cex=0.75, srt=90)
  abline(v=0, col='grey')
  z <- temp[,1]>0 | temp[,2]<0  # significant
  segments(x0=temp[,1], x1=temp[,2], y0=1:5, lwd=ifelse(z, 5, 3), 
           col=ifelse(z, pal, pal2))
  points(x, 1:5, pch=21, bg='white', col=ifelse(z, pal, pal2), lwd=2, cex=1.2)
}
legend(x=1.1*lim.x, y=35, xpd=NA, bty='n', cex=0.7, fill=rev(pal), 
       y.intersp=1.5, x.intersp=0.5,
       legend=rev(families), border=NA)
dev.off()

#plot(jitter(ael$transitivity), ael$overlap / ael$edge.count)
#summary(glm(cbind(overlap, edge.count) ~ edge.degree, data=ael, family='quasibinomial'))

# stepwise AIC model selection

#fit0 <- glm(cbind(overlap, edge.count) ~ 1, data=ael, family='binomial')
#step.fit <- stepAIC(fit, scope=list(lower=fit0, upper=fit), trace=F)



