# Inspecting data frames
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

# Removing small overlaps 
# non_splicing_overlaps <- overlaps %>%
#   filter(
#     ((shift == "+0" | shift == "-0") & overlap%%3 == 0) |
#       ((shift == "+1" | shift == "-1") & overlap%%3 == 2) |
#       ((shift == "+2" | shift == "-2") & overlap%%3 == 1),
#   )


orfs <- read.csv('total_orfs.csv')

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

