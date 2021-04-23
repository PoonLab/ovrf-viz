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
  
  max.size <- max(nl$size)
  ael$total.size <- nl$size[ael$parent] + nl$size[ael$child]
  
  ag <- graph_from_data_frame(ael)
  
  deg <- degree(ag)
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

res <- do.call("rbind", lapply(c(1, 3, 6, 9), function(i) {
  ael <- parse.edgelist(files[i])
  fit <- glm(cbind(overlap, edge.count) ~ 
               total.size + edge.degree + triangles + transitivity + centrality, 
             data=ael, family='quasibinomial')
  ci <- as.data.frame(confint(fit))
  ci$est <- fit$coefficients
  ci
}))

res$effect <- gsub("[0-9]$", "", row.names(res))
res$family <- rep(c('Adenoviridae', 'Coronaviridae', 'Papillomaviridae', 'Rhabdoviridae'), 
                   each=6)

require(ggfree)
pal <- hcl.colors(4, palette="Dark3")
pal2 <- hcl.colors(4, palette='Dark3', alpha=0.5)

pdf(file='~/papers/ovrf-viz/regression.pdf', width=5, height=5)
par(mfrow=c(5,1), mar=c(3,3,0,6), cex=1, xpd=F)
for (term in unique(res$effect)) {
  if (term == "(Intercept)") next
  temp <- res[res$effect==term, ]
  x <- temp$est
  lim.x <- max(abs(temp[,1:3]))
  plot(x, 1:4, xlim=c(-lim.x, lim.x), bty='n', ylim=c(0,5), yaxt='n', 
       ylab=NA, cex.axis=0.7, mgp=c(3,0.5,0))
  text(x=-0.95*lim.x, y=-1, label=gsub('\\.', ' ', term), 
       xpd=NA, adj=1, cex=0.8)
  abline(v=0, col='grey')
  z <- temp[,1]>0 | temp[,2]<0  # significant
  segments(x0=temp[,1], x1=temp[,2], y0=1:4, lwd=ifelse(z, 5, 3), 
           col=ifelse(z, pal, pal2))
  points(x, 1:4, pch=21, bg='white')
}
legend(x=lim.x, y=32, xpd=NA, bty='n', cex=0.8, fill=rev(pal), 
       y.intersp=1.5, x.intersp=0.5,
       legend=c('Adenoviridae', 'Coronaviridae', 'Papillomaviridae', 
                'Rhabdoviridae'), border=NA)
dev.off()

#plot(jitter(ael$transitivity), ael$overlap / ael$edge.count)
#summary(glm(cbind(overlap, edge.count) ~ edge.degree, data=ael, family='quasibinomial'))

# stepwise AIC model selection

#fit0 <- glm(cbind(overlap, edge.count) ~ 1, data=ael, family='binomial')
#step.fit <- stepAIC(fit, scope=list(lower=fit0, upper=fit), trace=F)



