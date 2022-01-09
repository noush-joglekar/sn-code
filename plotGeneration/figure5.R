# Intro ----------------------
# /bin/R
# By Anoushka Joglekar 10.2021
# Figures for snisor-seq paper
# Figure 5

# Setup ----------------------
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(ggsignif)

# Read input -----------------
oneCT <- read.table('data/exonCoord_BulkAndCT_status', header = T)
percSig <- read.table('data/percSig_byCelltype',header = T)
cu <-  read.table('data/consUseExons_pseudoBulk', header = T)
cu_ct <-  read.table('data/consUseExons_cellType', header = T)
disAs_dist_df <- read.table('data/diseaseAssoc', header = T)
PbyC <- read.table('data/phastCons_vs_inclusion_byCelltype',header = T)

# Functions ------------------
plotStatusDist <- function(statusMat){
  statusMatSubset <- as.matrix(statusMat[,5:9])
  
  statusMatSubset[statusMatSubset == "Same"] <- 3
  statusMatSubset[statusMatSubset == "Different"] <- 4
  statusMatSubset[statusMatSubset == "NonSig"] <- -2
  statusMatSubset[statusMatSubset == "cons"] <- 1
  statusMatSubset[statusMatSubset == "lowCounts"] <- 0
  
  statusMatSubset <- as.matrix(apply(statusMatSubset, 2, as.numeric))
  colnames(statusMatSubset) <- c("Astro","EN","IN","Oligo","OPCs")
  
  colors = structure(c("lightgrey","#1BB6AF","white","#B25D91","#172869"), names = c(0,1,-2,3,4)) 
  lgd = list(labels = c("Low Counts", "Constitutive", "Non Sig",
                        "Sig - Same Dir"), 
             title = "Exon \ncoordination", at =  c(0,1,-2,3),
             nrow = 2, title_position = "leftcenter")
  
  H0 = Heatmap(statusMatSubset,show_row_dend = FALSE, show_column_dend = FALSE,
               col = colors, name = "Direction of \nchange wrt Bulk", heatmap_legend_param = lgd)
  
  H = draw(H0, heatmap_legend_side = "bottom")
  
  return(H)
}

## Panel a ------------------

p_df <- as.data.frame(table(sapply(1:nrow(oneCT), function(ix) 
  length(which(oneCT[ix,5:9] =="Same")))))

p_df <- p_df %>% 
  arrange(desc(Var1)) %>%
  mutate(prop = Freq*100 / sum(p_df$Freq)) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

p_df$Var1 <- as.factor(p_df$Var1)


pie = ggplot(p_df, aes(x="", y=prop, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = paste0(Var1,"\n",round(prop,1),"%")), color = "white", size=4) +
  scale_fill_manual(values = c("#EFC7E6","#e3a7c0","#CB87B4","#B25D91","#7A0177","#49006A"))

### Output ------------------
pie

## Panel b ------------------

gBar <- ggplot(percSig, aes(x = CellType, y = PercSig)) + 
  geom_bar(stat = "identity",fill = "#B25D91", width = 0.7) + 
  geom_errorbar(aes(x=CellType, ymin=PercSig-SE, 
                    ymax=PercSig+SE),
                width = 0.05, color = "grey50") + 
  geom_hline(yintercept=25, linetype="dashed", color = "black") + 
  geom_hline(yintercept=50, linetype="dashed", color = "black") + 
  geom_hline(yintercept=75, linetype="dashed", color = "black") + 
  scale_y_continuous(breaks = seq(0, 75, by = 25)) +
  theme(legend.position = "none") + 
  theme_classic(base_size = 12) +
  labs(y = expression(italic(paste(" | ", log[2](oddsRatio), " | ")>= 1)),
       x = "Cell type") +
  theme(legend.position = "none")

### Output ------------------
gBar

## Panel c ------------------
untestability <- ggplot(cu, aes(x = exonType, y = percentCons, fill = exonType)) + 
  geom_bar(stat = "identity", width = 0.7) + 
  geom_errorbar(aes(x=exonType, ymin=percentCons-SE, 
                    ymax=percentCons+SE),
                width = 0.05, color = "grey60") + 
  geom_hline(yintercept=25, linetype="dashed", color = "black") + 
  geom_hline(yintercept=50, linetype="dashed", color = "black") + 
  geom_hline(yintercept=75, linetype="dashed", color = "black") + 
  theme(legend.position = "none") + 
  scale_y_continuous(breaks = seq(0, 75, by = 25)) +
  theme_classic(base_size = 12) +
  labs(y = "constitutive use of >=1 \nexon in >= 1 celltype (%)",
       x = "Exon type") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#6C6C9D","#172869"))

### Output ------------------
untestability

## Panel d ------------------

untestability_byCT <- ggplot(cu_ct %>% filter(value == "cons"), 
                             aes(x = name, y = percent, fill = exonType)) + 
  geom_bar(stat = "identity",position = "dodge", width = 0.7) + 
  geom_errorbar(aes(x=name, ymin=percent-SE, 
                    ymax=percent+SE),
                width = 0.05, color = "grey60",position = position_dodge(.9)) + 
  geom_hline(yintercept=25, linetype="dashed", color = "black") + 
  geom_hline(yintercept=50, linetype="dashed", color = "black") + 
  theme_classic(base_size = 12) +
  scale_y_continuous(breaks = seq(0, 50, by = 25)) +
  labs(y = "Untestable genes \nwith constitutive exons (%)", x = "") +
  theme(axis.text.x = element_text(angle = 90)) + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#6C6C9D","#172869"))

### Output ------------------
untestability_byCT

## Panel e ------------------
matCounts <- as.matrix(cbind(disAs_dist_df[,2]-disAs_dist_df[,3],disAs_dist_df[,3]))
sigLevel <- format(fisher.test(matCounts)$p.value,digits = 3)

da <- ggplot(disAs_dist_df, aes(x = Association, y = Coordinated *100 / Tested, fill = Association)) +
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(x=Association, ymin=Perc-SE, 
                    ymax=Perc+SE),
                width = 0.05, color = "grey50") +
  geom_signif(annotations = sigLevel,
              y_position = 43, xmin = 1, xmax = 2) +
  scale_fill_manual(values = c("#4D9A97","#1BB6AF")) + 
  labs(y = "Coordinated (%)", x = "Disease association") +
  theme_classic(base_size = 14)

### Output ------------------
da

## Panel g ------------------
H.oneCT <- plotStatusDist(oneCT)

### Output ------------------
H.oneCT


## Panel h and i -------------

PC_absLOR <- ggplot(PbyC %>% filter(cellType == "ExcitatoryNeurons"), 
                    aes(x = MinPhastcons, y = logOddsRatio, fill = as.factor(group2))) + 
  geom_boxplot(aes(factor(group2)), color = c("black","grey60")) +
  geom_signif(aes(factor(group2)),
              comparisons = list(c(1,2)),
              map_signif_level = F) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = c("#CB87B4","#172869")) +
  theme(legend.position = "none") +
  labs(x = "minimum Phastcons", 
       y = expression(italic( paste(" | ", log[2](oddsRatio), " | "))))

PC_doubleIn <- ggplot(PbyC %>% filter(cellType == "ExcitatoryNeurons"), 
                      aes(x = MinPhastcons, y = doublein, fill = as.factor(group2))) + 
  geom_boxplot(aes(factor(group2)),color = c("black","grey60")) +
  geom_signif(aes(factor(group2)),
              comparisons = list(c(1,2)),
              map_signif_level = F) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = c("#CB87B4","#172869")) +
  theme(legend.position = "none") +
  labs(x = "minimum Phastcons", 
       y = expression(italic( inclusion / conserved)))

### Output ------------------
PC_absLOR | PC_doubleIn
