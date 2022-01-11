# Intro ----------------------
# /bin/R
# By Anoushka Joglekar 10.2021
# Figures for snisor-seq paper
# Figure 6

# Setup ----------------------
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)


# Read input -----------------
cdf <- read.table('data/percSig_endSite_exon_pseudobulk', header = T)
tss_df <- read.table('data/percSig_tss_exon', header = T)
tss_superMat_sB <- as.matrix(read.table('data/tss_exon_tertiaryRep', header = T))
polya_df <- read.table('data/percSig_polya_exon', header = T)
polya_superMat_sB <- as.matrix(read.table('data/polya_exon_tertiaryRep', header = T))


## Panel a -------------------
tss_df$Var1 <- as.factor(tss_df$Var1)

pie_tss = ggplot(tss_df, aes(x="", y=prop, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = paste0(Var1,"\n",round(prop,1),"%")), color = "white", size=4) +
  scale_fill_manual(values = c("#EFC7E6","#e3a7c0","#CB87B4","#B25D91","#7A0177","#49006A"))

### Output -------------------
pie_tss


## Panel b -------------------
g_dPI_t = ggplot(cdf %>% filter(EndSite == "TSS"), 
                 aes(x=dPI,y=PercentageSig,fill=EndSite))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(x=dPI, ymin=PercentageSig-StandardError, 
                    ymax=PercentageSig+StandardError),
                width = 0.05, color = "grey50") + 
  geom_hline(yintercept=5, linetype="dashed", color = "black") + 
  geom_hline(yintercept=10, linetype="dashed", color = "black") + 
  scale_fill_manual(values = c("#172869")) + 
  scale_x_continuous(breaks = seq(0, 0.7, by = 0.1)) + 
  ylim(c(0,10)) + 
  theme(legend.position = "none") + 
  theme_classic(base_size = 12) +
  labs(x = expression(italic( paste(" | ",Delta * Pi ," | ")))) +
  theme(legend.position = "none",
        axis.title.y = element_blank())

### Output -------------------
g_dPI_t

## Panel c -------------------
colors = structure(c("lightgrey","#1BB6AF","white","#B25D91"), names = c(0,1,-2,3))

lgd = list(labels = c("lowCounts","Constitutive", "Non Sig",
                      "Sig"), 
           title = "TSS - Exon \ncoordination", at =  c(0, 1,-2,3),
           nrow = 2, title_position = "leftcenter")

H_tmp = Heatmap(tss_superMat_sB, show_row_names = FALSE, col = colors, 
                cluster_rows = T, show_row_dend = F, cluster_columns = F,
                name = "Sig in Bulk by celltype",heatmap_legend_param = lgd,
                show_column_names = FALSE)

H_tss = draw(H_tmp, heatmap_legend_side = "bottom")

### Output -------------------
H_tss

## Panel d -------------------
polya_df$Var1 <- as.factor(polya_df$Var1)

pie_polya = ggplot(polya_df, aes(x="", y=prop, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = paste0(Var1,"\n",round(prop,1),"%")), color = "white", size=4) +
  #scale_fill_brewer(palette="GnBu",type = "div", direction = -1) +
  scale_fill_manual(values = c("#EFC7E6","#e3a7c0","#CB87B4","#B25D91","#7A0177","#49006A"))


### Output -------------------
pie_polya

## Panel e -------------------
g_dPI_p = ggplot(cdf %>% filter(EndSite == "PolyA"), 
                 aes(x=dPI,y=PercentageSig,fill=EndSite))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(x=dPI, ymin=PercentageSig-StandardError, 
                    ymax=PercentageSig+StandardError),
                width = 0.05, color = "grey50") + 
  geom_hline(yintercept=5, linetype="dashed", color = "black") + 
  geom_hline(yintercept=10, linetype="dashed", color = "black") + 
  scale_fill_manual(values = c("#6C6C9D")) + 
  scale_x_continuous(breaks = seq(0, 0.7, by = 0.1)) + 
  ylim(c(0,10)) + 
  theme(legend.position = "none") + 
  theme_classic(base_size = 12) +
  labs(x = expression(italic( paste(" | ",Delta * Pi ," | ")))) +
  theme(legend.position = "none",
        axis.title.y = element_blank())

### Output -------------------
g_dPI_p

## Panel f -------------------
lgd = list(labels = c("lowCounts","Constitutive", "Non Sig",
                      "Sig"), 
           title = "PolyA - Exon \ncoordination", at =  c(0, 1,-2,3),
           nrow = 2, title_position = "leftcenter")

H_tmp2 = Heatmap(polya_superMat_sB, show_row_names = FALSE, col = colors, 
                 cluster_rows = T, show_row_dend = F, cluster_columns = F,
                 name = "Sig in Bulk by celltype",heatmap_legend_param = lgd)

H_polyA = draw(H_tmp2, heatmap_legend_side = "bottom")

### Output -------------------
H_polyA

### Panel g ------------------
### Run independently using ScisorWiz - code not included here 
