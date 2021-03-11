setwd('/home/laura/Projects/ovrf-review/data/kmer_distance_all')


filenames <- list.files(full.names=TRUE)  
All <- lapply(filenames,function(i){
  read.csv(i, header=FALSE, row.names= 1)
})

for(i in 1:length(All)) {
  d <- All[i][[1]]
  up <- upper.tri(d)
  m<- mean(d[up])
  print(filenames[i])
  print(m)
}
