library(dplyr)
library(treemap)
library(plyr)
library(ggplot2)
library(ggridges)
require(grid)
########################################################
# Plot frame shift distribution
########################################################
setwd('/home/lmunoz/Projects/ovrf-review/dataset_2020/')
#overlaps <- read.csv('2_overlapping_march2020.csv', colClasses='character')
overlaps <- read.csv('overlaps_with_noise.csv', colClasses='character')

overlaps$extreme_left1 <- as.integer(overlaps$extreme_left1)
overlaps$extreme_left2 <- as.integer(overlaps$extreme_left2)
overlaps$dir1 <- as.factor(overlaps$dir1)
overlaps$extreme_right1 <- as.integer(overlaps$extreme_right1)
overlaps$extreme_right2 <- as.integer(overlaps$extreme_right2)
overlaps$seqlen1 <- as.integer(overlaps$seqlen1)
overlaps$seqlen2 <- as.integer(overlaps$seqlen2)
overlaps$overlap <- as.integer(overlaps$overlap)
overlaps$shift <- factor(overlaps$shift, 
                         levels=c('-2', '-1', '-0', '+0', '+1', '+2'))
non_filtered_ovps<- overlaps

### Find duplicates with different overlap lengths and merge them (remove 1426 entries)
overlaps2 <- overlaps %>% group_by(across(c(-overlap))) %>% summarise(overlap = sum(overlap))

# Remove 5241 entries where number of nucleotides is not consequent with frame shift
non_splicing_overlaps <- overlaps %>%
  filter(
    ((shift == "+0" | shift == "-0") & overlap%%3 == 0) |
    ((shift == "+1" | shift == "-1") & overlap%%3 == 2) |
    ((shift == "+2" | shift == "-2") & overlap%%3 == 1)
  )

# Remove 49 entries with circular and alternative splacing genomes
circular <- which(non_splicing_overlaps$extreme_left1==non_splicing_overlaps$extreme_left2
                  & non_splicing_overlaps$extreme_right1==non_splicing_overlaps$extreme_right2)

non_splicing_overlaps <- non_splicing_overlaps[-circular,]

overlaps <- non_splicing_overlaps

# Ridgeplot
plot <- ggplot(non_splicing_overlaps, aes(x = log10(overlap), y = shift, group = shift, fill = shift)) +
  geom_density_ridges()+
  ggtitle("non_splicing_overlaps length per frame shift")+
  scale_fill_manual(values = c("#c1462b", "#de8f42" ,"#f1dba2", "#cacfa6", "#769845" ,"#32621f")) +
  theme_bw()
plot

#palete<-c("#32621f", "#769845", "#cacfa6", "#f1dba2", "#de8f42", "#c1462b")
# +2, +1, +0, +1, -1, -2

# Introduce Noise: Randomly select rows and change coordinates 
n<-3
modified_ov <- overlaps %>% 
                        slice_sample(n=n)

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
# Summary of overlap trends
########################################################
nss_RNA <- total[which(total$baltimore.class == "ss_RNA_-"),]
pss_RNA <- total[which(total$baltimore.class == "ss_RNA_+"),]
ds_RNA <- total[which(total$baltimore.class == "ds_RNA"),]
rt <- total[which(total$baltimore.class == "RT_viruses"),]

# Cero frameshift
pcero <- total[which(total$shift == "+0"),]
replicases<-cero[which(grepl("replica", cero$prod1)),]
capsids<-cero[which(grepl("capsid", cero$prod1)),]
structural <-cero[which(grepl("struct", cero$prod1)),]
polymerases <-cero[which(grepl("polyme", cero$prod1)),]

# Antisense
anti<- total[which(grepl("-", total$shift)),]
none <- anti[which(grepl("1", anti$shift)),]
ntwo <- anti[which(grepl("2", anti$shift)),]
ncero <- anti[which(grepl("0", anti$shift)),]

res <- subset(pss_RNA, shift=="+0")
tsorted<-ntwo[order(ntwo$overlap, decreasing=TRUE),]
nRNA <- total[which(grepl("RNA", total$baltimore.class) & grepl("-", total$shift)),]
total[which(grepl("RNA", total$baltimore.class) & total$shift=="-1"),]

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

temp <- sapply(split(overlaps$overlap, overlaps$accn), length)


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

#write.csv(virus, "out_virus_df_inspection.csv")
#write.csv(total, "overlap_info.csv")

# Calculate number of proteins for each virus based on orf data frame
nprot <- data.frame(count(orfs, "accno"))
prot_index <- match(nprot$accno, virus$Representative)
virus_combined <- virus %>% left_join(nprot, by=c("Representative" = "accno"))
plot(virus$Length, virus$freq, log='x')
 
virus_combined$Proteins<- virus_combined$freq


##################################################
# Virus stats Dr Poon
##################################################
#setwd('/home/lmunoz/Projects/ovrf-review/dataset_2020')
setwd('~/git/ovrf-review/dataset_2020')
virus <- read.csv('out_virus_df_inspection.csv', stringsAsFactors = F)

virus$n.overlaps[which(is.na(virus$n.overlaps))] <- 0
#virus$len.overlaps[which(is.na(virus$len.overlaps))] <- 0.5
#virus$total.overlaps[which(is.na(virus$total.overlaps))] <- 0

virus$rel.ovrf <- virus$n.overlaps / virus$Proteins
virus$rel.ovrf <- virus$n.overlaps / virus_combined$Proteins
virus$baltimore.class[which(virus$baltimore.class == "Unknown" & grepl("RNA", virus$Molecule))] <- "unknown_RNA"
virus$baltimore.class[which(virus$baltimore.class == "Unknown" & grepl("DNA", virus$Molecule))] <- "unknown_DNA"


cor.test(virus$n.overlaps, virus$Proteins, method='spearman')

plot(virus$len.overlaps, virus$Length, log='xy')
cor.test(virus$len.overlaps, virus$Length, method='spearman')

for (sub in split(virus, virus$baltimore.class)) {
  if (all(sub$baltimore.class=='circular_ss_RNA')) {
    next
  }
  res <- cor.test(sub$len.overlaps, sub$Length, method='spearman')
  print(unique(sub$baltimore.class))
  print(res)
}


pdf(file='~/papers/ovrf-viz/sm1.pdf', width=8, height=8)
par(mfrow=c(3,3), mar=c(5,5,1,1))
for (sub in split(virus, virus$baltimore.class)) {
  if (all(sub$baltimore.class=='circular_ss_RNA')) {
    next
  }
  plot(sub$len.overlaps, sub$Length, log='xy', 
       xlab='Mean overlap length (nt)', ylab='Genome length (nt)',
       type='n')
  tab <- table(sub$family)
  pal <- hcl.colors(length(tab[tab>1]), 'Spectral')
  set.seed(1)
  pal <- sample(pal, length(pal))
  counter <- 1
  for (ssub in split(sub, sub$family)) {
    if (nrow(ssub) > 1) {
      points(ssub$len.overlaps, ssub$Length, bg=pal[counter], pch=21, cex=0.8)
      counter <- counter + 1
    } else {
      points(ssub$len.overlaps, ssub$Length, col='grey', cex=0.5)
    }
  }
  res <- cor.test(sub$len.overlaps, sub$Length, method='spearman')
  u <- par("usr")
  text(x=10^(u[1] + 0.05*(u[2]-u[1])), 
       y=10^(u[3] + 0.1*(u[4]-u[3])), 
       label=paste("rho", round(res$estimate, 2), "\nP=", round(res$p.value, 2)), 
       adj=0)
  #legend(x=min(sub$len.overlaps)+100, y=min(sub$Length)+1, pch=21, 
  #       legend=names(sort(tab, decreasing = T))[1:3], yjust=0)
  
  title(main=unique(sub$baltimore.class), adj=0)
}
dev.off()


dsDNA <- virus[virus$baltimore.class=='ds_DNA',]
for (sub in split(dsDNA, dsDNA$family)) {
  if (nrow(sub) > 100) {
    print(unique(sub$family))
    print(cor.test(sub$Length, sub$len.overlaps, method='spearman'))
  }
}


require(lme4)
fit <- lmer(log(len.overlaps) ~ log(Length) + (1 | family), 
            data=virus[virus$baltimore.class != 'circular_ss_RNA', ])


temp <- virus[virus$baltimore.class != 'circular_ss_RNA', ]
by(temp, temp$baltimore.class, function(x) cor.test(log(x$len.overlaps), log(x$Length)))


# significant variation among Baltimore classes
fit <- glm(rel.ovrf ~ baltimore.class, data=virus[virus$baltimore.class != "circular_ss_RNA",])
anova(fit, test='F')

summary(virus[grepl("ds_DNA", virus$baltimore.class),])
summary(virus[grepl("DNA", virus$baltimore.class),])
summary(virus[grepl("RNA", virus$baltimore.class),])
summary(virus[grepl("ss_RNA_+", virus$baltimore.class),])
summary(virus[grepl("ss_RNA_-", virus$baltimore.class),])

wilcox.test(virus$Length[virus$n.overlaps==0], virus$Length[virus$n.overlaps > 0])
t.test(log(virus$Length[virus$n.overlaps==0]), log(virus$Length[virus$n.overlaps > 0]))
boxplot(log(virus$Length[virus$n.overlaps==0]), log(virus$Length[virus$n.overlaps > 0]))

##################################################
# Frame shift stats (total)
##################################################
setwd('/home/lmunoz/Projects/ovrf-review/dataset_2020/')
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

####################################################
# Checking for Zero entries
####################################################

zero <- total[which(total$shift == "+0"),]
table(zero$extreme_left1==zero$extreme_left2, zero$extreme_right1==zero$extreme_right2)
head(zero[(zero$extreme_left1==zero$extreme_left2 & zero$extreme_right1==zero$extreme_right2),])
cases<-zero[which((zero$extreme_left1==zero$extreme_left2 & zero$extreme_right1==zero$extreme_right2)),]


extreme_left1 == extreme_left2 == 0


