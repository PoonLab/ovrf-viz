########################################################
# Plot frame shift distribution
########################################################
setwd('/home/lmunoz/Projects/ovrf-review/dataset_2020/')
overlaps <- read.csv('2_overlapping_march2020.csv', colClasses='character')

overlaps$loc1 <- as.integer(overlaps$loc1)
overlaps$dir1 <- as.factor(overlaps$dir1)
overlaps$loc2 <- as.integer(overlaps$loc2)
overlaps$dir2 <- as.factor(overlaps$dir2)
overlaps$seqlen1 <- as.integer(overlaps$seqlen1)
overlaps$seqlen2 <- as.integer(overlaps$seqlen2)
overlaps$overlap <- as.integer(overlaps$overlap)
overlaps$shift <- factor(overlaps$shift, 
                         levels=c('-2', '-1', '-0', '+0', '+1', '+2'))


# Remove entries where number of nucleotides is not consequent with frame shift
non_splicing_overlaps <- overlaps %>%
  filter(
    ((shift == "+0" | shift == "-0") & overlap%%3 == 0) |
    ((shift == "+1" | shift == "-1") & overlap%%3 == 2) |
    ((shift == "+2" | shift == "-2") & overlap%%3 == 1)
  )


# Ridgeplot
plot <- ggplot(non_splicing_overlaps, aes(x = log10(overlap), y = shift, group = shift, fill = shift)) +
  geom_density_ridges()+
  ggtitle("non_splicing_overlaps length per frame shift")+
  scale_fill_manual(values = c("#c1462b", "#de8f42" ,"#f1dba2", "#cacfa6", "#769845" ,"#32621f")) +
  theme_bw()
plot

#palete<-c("#32621f", "#769845", "#cacfa6", "#f1dba2", "#de8f42", "#c1462b")
# +2, +1, +0, +1, -1, -2

library(treemap)
library(plyr)
require(grid)
colnames(non_splicing_overlaps)[1] <- "Accession"
colnames(virus)[3] <- "Accession"
total <- merge(x=non_splicing_overlaps, y=virus[,c("Accession", "baltimore.class")], by="Accession", all.x =T)
tbl <- as.data.frame(table(total$shift, total$baltimore.class))
colnames(tbl) <- c("shift", "balt", "Freq")
#new_p<-c("#9d2503", "#df8543", "#f6e0a7", "#90a95c","#5a715c" )
new_p<-c("#c1462b", "#de8f42" ,"#f1dba2", "#cacfa6", "#769845" ,"#32621f")

#pdf(file="treemap2.#pdf")

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2)))
sub <- subset(tbl, balt == "ds_DNA")
treemap(sub, index = "shift", vSize = "Freq", type="index", 
        palette=new_p, title="ds_DNA", fontsize.title=15, fontsize.labels=14,
        vp = vplayout(1,1))

sub <- subset(tbl, balt == "ds_RNA")
treemap(sub, index = "shift", vSize = "Freq", type="index", 
        palette=new_p, title="ds_RNA", fontsize.title=15, fontsize.labels=14,
        vp = vplayout(1,2))

sub <- subset(tbl, balt == "ss_DNA")
treemap(sub, index = "shift", vSize = "Freq", type="index", 
        palette=new_p, title="ss_DNA", fontsize.title=15, fontsize.labels=14,
        vp = vplayout(2,1))

sub <- subset(tbl, balt == "ss_RNA_-")
treemap(sub, index = "shift", vSize = "Freq", type="index", 
        palette=new_p, title="ss_RNA_-", fontsize.title=15, fontsize.labels=14,
        vp = vplayout(2,2))

sub <- subset(tbl, balt == "ss_RNA_+")
treemap(sub, index = "shift", vSize = "Freq", type="index", 
        palette=new_p, title="ss_RNA_+", fontsize.title=15, fontsize.labels=14,
        vp = vplayout(3,1))

sub <- subset(tbl, balt == "RT_viruses")
treemap(sub, index = "shift", vSize = "Freq", type="index", 
        palette=new_p, title="RT_viruses", fontsize.title=15, fontsize.labels=14,
        vp = vplayout(3,2))

########################################################
# Import and merge data
########################################################
setwd('/home/lmunoz/Projects/ovrf-review/dataset_2020/')
overlaps <- non_splicing_overlaps
orfs <- read.csv('total_orfs.csv')
virus <- read.csv('try_baltimore_class.csv')

virus$Proteins[virus$Proteins=='-'] <- NA
virus$Proteins <- as.integer(as.character(
  virus$Proteins))
summary(virus)


# the number of ORFs increases linearly with genome size
plot(virus$Length, virus$Proteins, log='x')

temp <- sapply(split(overlaps$overlap, overlaps$Accession), length)
noverlaps <- data.frame(accn=names(temp), count=temp)

noverlaps$mean.olen <- sapply(split(overlaps$overlap, overlaps$Accession), mean)
noverlaps$total.olen <- sapply(split(overlaps$overlap, overlaps$Accession), sum)

# sum total number of overlaps per genome
virus$first.acc <- sapply(virus$Representative, function(x) {
  xx <- as.character(x)
  if (grepl("^\\[", xx)) {
    gsub("^\\['([A-Z0-9_]+)'.+", "\\1", xx)
  } else {
    xx
  }
})

index <- match(virus$Representative, noverlaps$accn)
virus$n.overlaps <- noverlaps$count[index] # 
virus$total.overlaps<-noverlaps$total.olen[index] #Total number of overlaps


# carry over Genome to <noverlaps>
all.acc <- sapply(1:nrow(virus), function(i) {
  xx <- as.character(virus$Representative[i])
  if (grepl("^\\[", xx)) {
    x2 <- gsub("[\\]'\\[]", "", xx, perl=T)
    strsplit(x2, ", ")[[1]]
  } else {
    xx
  }
})

all.acc <- data.frame(acc=unlist(all.acc), genome=rep(virus$Taxonomy.name, times=sapply(all.acc, length)))

index2 <- match(noverlaps$accn, all.acc$acc)
noverlaps$genome <- all.acc$genome[index2]

temp <- sapply(split(noverlaps$mean.olen/noverlaps$count, noverlaps$genome), mean)
index3 <- match(virus$Taxonomy.name, names(temp))
virus$len.overlaps <- ifelse(virus$n.overlaps==0, NA, temp[index3] * virus$n.overlaps)

write.csv(virus, "out_virus_df_inspection.csv")
write.csv(total, "overlap_info.csv")

##################################################
# Virus stats Dr Poon
##################################################

virus <- read.csv('/home/lmunoz/Projects/ovrf-review/dataset_2020/out_virus_df_inspection.csv')
virus$rel.ovrf <- virus$n.overlaps / virus$Proteins

plot(virus$len.overlaps, virus$Length, log='xy', 
     col=rainbow(8)[as.factor(virus$baltimore.class)])

cor.test(log(virus$len.overlaps), log(virus$Length))
cor.test(log(virus$n.overlaps), log(virus$Proteins))

for (sub in split(virus, virus$baltimore.class)) {
  if (all(sub$baltimore.class=='circular_ss_RNA')) {
    next
  }
  res <- cor.test(log(sub$len.overlaps), log(sub$Length))  
  print(res)
}

temp <- virus[virus$baltimore.class != 'circular_ss_RNA', ]
by(temp, temp$baltimore.class, function(x) cor.test(log(x$len.overlaps), log(x$Length)))


# significant variation among Baltimore classes
fit <- glm(rel.ovrf ~ baltimore.class, data=virus[virus$baltimore.class != "circular_ss_RNA",])
anova(fit, test='F')

summary(virus[grepl("DNA", virus$baltimore.class),])
summary(virus[grepl("RNA", virus$baltimore.class),])

##################################################
# Frame shift stats (total)
##################################################
virus <- read.csv('out_virus_df_inspection.csv', stringsAsFactors = F)

#Re naming Unknown according to molecule type
virus$baltimore.class[which(virus$baltimore.class == "Unknown" & grepl("RNA", virus$Molecule))] <- "unknown_RNA"
virus$baltimore.class[which(virus$baltimore.class == "Unknown" & grepl("DNA", virus$Molecule))] <- "unknown_DNA"

colnames(non_splicing_overlaps)[1] <- "Accession"
colnames(virus)[3] <- "Accession"
total <- merge(x=non_splicing_overlaps, y=virus[,c("Accession", "baltimore.class")], by="Accession", all.x =T)

plus_two<-(subset(total, shift =="+2"))
summary(plus_two)

minus_two<-(subset(total, shift =="-2"))
summary(minus_two)

plus_one<-subset(total, shift =="+1")
minus_one<-subset(total, shift=="-1")

plus_cero<-subset(total, shift=="+0")
minus_cero<-subset(total, shift=="-0")
