# Intro ----------------------
# /bin/R
# By Anoushka Joglekar 10.2021
# Figures for snisor-seq paper
# Figure 4

# Setup ----------------------
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(ggsignif)

# Read input -----------------
bulk_exon_dist <- read.table('data/DIST.maxabsLOR.perGene.tab',header = T)
bulk_exon_adj <- read.table('data/ADJ.maxabsLOR.perGene.tab',header = T)
mc <- read.table('data/adjExons_coordStatus_lengthOfIntron', header = T)
mm1 <- read.table('data/adjExons_spliceSiteStrength_geneID',header = T)
mm2 <- read.table('data/adjExons_spliceSiteStrength_maxEnt',header = T)
tss_exon_pC <- read.table('data/TSS_exon_mindpi_phastcons_17way.tab',header = T)
polya_exon_pC <- read.table('data/TSS_exon_mindpi_phastcons_17way.tab',header = T)


# Functions ------------------
populateMatrix_oR <- function(inMat,colID){
  tot = nrow(inMat)
  L = list()
  for(cutoff in 1:7){
    sP <- nrow(inMat %>% filter(FDR <= 0.05 & abs(logOddsRatio) >= as.numeric(cutoff)))
    perc <- sP*100/tot
    se <- sqrt(perc*(100-perc)/tot)
    L[[cutoff]] <- list(cutoff,sP*100/tot,se,colID)
  }
  return(L)
}

## Panels b and c ----------------
d <- populateMatrix_oR(bulk_exon_dist,"Distant")
a <- populateMatrix_oR(bulk_exon_adj,"Adjacent")

cdf_oR <- as.data.frame(rbind(do.call("rbind",a),do.call("rbind",d)))
colnames(cdf_oR) <- c("LOR","PercentageSig","StandardError","ExonPairType")
cdf_oR$PercentageSig <- as.numeric(as.vector(cdf_oR$PercentageSig))
cdf_oR$StandardError <- as.numeric(as.vector(cdf_oR$StandardError))
cdf_oR$LOR <- as.numeric(as.vector(cdf_oR$LOR))
cdf_oR$ExonPairType <- factor(cdf_oR$ExonPairType, levels = c("Adjacent","Distant"))

g_oR_a = 
  ggplot(cdf_oR %>% filter(ExonPairType == "Adjacent"), 
         aes(x=LOR,y=PercentageSig,fill=ExonPairType))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(x=LOR, ymin=PercentageSig-StandardError, 
                    ymax=PercentageSig+StandardError),
                width = 0.4, color = "grey50") + 
  geom_hline(yintercept=25, linetype="dashed", color = "black") + 
  geom_hline(yintercept=50, linetype="dashed", color = "black") + 
  scale_x_continuous(breaks = seq(0, 7, by = 1)) +
  scale_y_continuous(breaks = seq(0, 75, by = 25)) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = c("#F1BB7B")) + 
  labs(x = expression(italic( paste(" | ", log[2](oddsRatio), " | "))),
       y = expression(italic("Tested genes with \ncoordination >= threshold (%)"))) +
  theme(legend.position = "none")

g_oR_d = 
  ggplot(cdf_oR %>% filter(ExonPairType == "Distant"), 
         aes(x=LOR,y=PercentageSig,fill=ExonPairType))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(x=LOR, ymin=PercentageSig-StandardError, 
                    ymax=PercentageSig+StandardError),
                width = 0.4, color = "grey50") + 
  geom_hline(yintercept=25, linetype="dashed", color = "black") + 
  geom_hline(yintercept=50, linetype="dashed", color = "black") + 
  scale_x_continuous(breaks = seq(0, 7, by = 1)) +
  scale_y_continuous(breaks = seq(0, 75, by = 25)) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = c("#FA796C")) + 
  labs(x = expression(italic( paste(" | ", log[2](oddsRatio), " | "))),
       y = "") +
  theme(legend.position = "none")

### Output -----------------------
g_oR_a | g_oR_d

## Panels e and f ----------------

a <- bulk_exon_adj %>% filter(FDR <= 0.05)
d <- bulk_exon_dist %>% filter(FDR <= 0.05)

aLOR <- as.data.frame(rbind(cbind(a$logOddsRatio,"Adjacent"),
                            cbind(d$logOddsRatio,"Distant")))
colnames(aLOR) <- c("LOR","ExonType")
aLOR$LOR <- as.numeric(as.vector(aLOR$LOR))
aLOR$ExonType <- factor(aLOR$ExonType, levels = c("Adjacent","Distant"))

ylim1 = boxplot.stats(abs(aLOR$LOR))$stats[c(1, 5)]

bp = ggplot(aLOR,aes(y = abs(LOR), x = ExonType, fill = ExonType, color = ExonType)) +
  geom_boxplot(outlier.shape = NA, alpha =0.5) + 
  theme_classic(base_size = 12) +
  scale_fill_manual(values = c("#F1BB7B", "#FA796C")) +
  scale_color_manual(values = c("#F1BB7B", "#FA796C")) +
  geom_signif(comparisons = list(c("Adjacent","Distant")),
              y_position = 18.5,color = "grey50") + 
  theme(legend.position = "bottom") + 
  labs(y = expression(italic( paste(" | ", log[2](oddsRatio), " | "))),
       x = "Exon Type") + 
  coord_cartesian(ylim = ylim1*1.1)

dp = ggplot(aLOR, aes(x = LOR,fill = ExonType, color = ExonType)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0,linetype="dashed", color = "grey50") +
  scale_fill_manual(values = c("#F1BB7B", "#FA796C")) +
  scale_color_manual(values = c("#F1BB7B", "#FA796C")) +
  theme_classic(base_size = 12) + theme(legend.position = "bottom") + 
  labs(y = "Frequency",
       x = expression(italic( log[2](oddsRatio))))

### Output -----------------------
bp | dp 

## Panel g  ----------------------

anno <- sapply(as.vector(unique(mc$variable)), function(var)
  wilcox.test(as.vector(mc %>% filter(Type == "Coordinated" & variable == var) %>% select(value))$value,
            as.vector(mc %>% filter(Type == "Non-Coordinated" & variable == var) %>% select(value))$value)$p.value)

g1 = ggplot(mc, aes(x = variable, y = value, fill = Type)) +
  geom_boxplot(position = "dodge", outlier.alpha = 0.5, outlier.size = 0.2) + 
  theme_classic(base_size = 15) + 
  ylim(0,8e04) +
  geom_signif(annotation=formatC(anno, digits=2),
              y_position=7.5e4, xmin=seq(0.75,2.75, 1), 
              xmax=seq(1.25,3.25, 1), 
              tip_length = 0) +
  scale_fill_manual(values = c("#5B1A18","#FC8D62")) + 
  xlab("Position of intron") +
  ylab("Length of intron")

### Output -----------------------
g1

## Panel h -----------------------

### Gene ID
anno_m1 <- sapply(c("Acc1","Don1","Acc2","Don2"), function(ss) 
  wilcox.test(as.vector(mm1 %>% filter(Type == "Non-Coordinated" & variable == ss) %>% select(value))$value,
              as.vector(mm1 %>% filter(Type == "Coordinated" & variable == ss) %>% select(value))$value)$p.value)

gm1 <- ggplot(mm1, aes(x = variable, y = value, color = Type)) + 
  geom_boxplot(outlier.alpha = 0.5, lwd = 0.7) +
  theme_classic(base_size = 15) +
  geom_signif(annotation=formatC(anno_m1, digits=2),
              y_position=14, xmin=seq(0.8,3.8, 1), 
              xmax=seq(1.2,4.2, 1), 
              color = "black") +
  scale_color_manual(values = c("#5B1A18","#FC8D62")) +
  xlab("Splice Site") +
  ylab("Gene ID score") +
  theme(legend.position = "bottom")

### Max Ent

anno_m2 <- sapply(c("Acc1","Don1","Acc2","Don2"), function(ss) 
  wilcox.test(as.vector(mm2 %>% filter(Type == "Non-Coordinated" & variable == ss) %>% select(value))$value,
              as.vector(mm2 %>% filter(Type == "Coordinated" & variable == ss) %>% select(value))$value)$p.value)

gm2 <- ggplot(mm2, aes(x = variable, y = value, color = Type)) + 
  geom_boxplot(outlier.alpha = 0.5, lwd = 0.7) +
  theme_classic(base_size = 15) +
  geom_signif(annotation=formatC(anno_m2, digits=2),
              y_position=14, xmin=seq(0.8,3.8, 1), 
              xmax=seq(1.2,4.2, 1), 
              color = "black") +
  scale_color_manual(values = c("#5B1A18","#FC8D62")) +
  xlab("Splice Site") +
  ylab("Max Ent score") +
  theme(legend.position = "bottom")

### Output -----------------------
gm1 | gm2 


## Panel i -----------------------
exon_pc <- bulk_exon_adj
exon_pc <- exon_pc %>% mutate(Status = case_when(logOddsRatio >= 0 ~ "Pos", TRUE ~ "Neg"))
exonPos = exon_pc %>% filter(Status == "Pos") %>% select(phastcons_min,logOddsRatio)
cor_exon_pos <- cor.test(exonPos$phastcons_min,exonPos$logOddsRatio)

pC1 <- ggplot(data = exonPos,
              aes(phastcons_min, abs(logOddsRatio))) +
  geom_point(col = "#F8866F") + 
  geom_smooth(method = 'lm',color = "brown4") +
  labs(x = "Primate Phastcons Score",
       y = expression(italic( log[2](oddsRatio)))) + 
  theme_classic(base_size = 12) +
  annotate(geom = "text",x = 0.25,y = 15,
           label = paste(expression(r^2),"=",round(cor_exon_pos$estimate^2,2),"\np-value < 0.05")) +
  theme(legend.position = "none")

### Output --------------------------
pC1

## Panel j --------------------------
### TSS

cor_tss <- cor.test(tss_exon_pC$phastcons,abs(tss_exon_pC$dPI))

pT <- ggplot(data = tss_exon_pC,
             aes(phastcons, abs(dPI))) +
  geom_point(col = "#D67236") + 
  geom_smooth(method = 'lm',color = "#7B2827") +
  labs(x = "Primate Phastcons Score",
       #y = expression(italic( paste(" | ", Log[2](oddsRatio), " | ")))) + 
       y = expression(italic( paste(" | ", Delta * Pi, " | ")))) + 
  theme_classic(base_size = 12) + 
  annotate(geom = "text",x = 0.75,y = 0.95,
           label = c(paste(expression(r^2)," =",round(cor_tss$estimate^2,3), "\np-value > 0.05"))) +
  theme(legend.position = "none")

### PolyA
cor_polya <- cor.test(polya_exon_pC$phastcons,abs(polya_exon_pC$dPI))


pP <- ggplot(data = polya_exon_pC,
             aes(phastcons, abs(dPI))) +
  geom_point(col = "#F1BB7B") + 
  geom_smooth(method = 'lm',color = "#7B2827") +
  labs(x = "Primate Phastcons Score",
       #y = expression(italic( paste(" | ", log[2](oddsRatio), " | ")))) + 
       y = expression(italic( paste(" | ", Delta * Pi, " | ")))) + 
  theme_classic(base_size = 12) + 
  annotate(geom = "text",x = 0.75,y = 0.95,
           label = c(paste(expression(r2)," =",round(cor_polya$estimate^2,3), "\np-value < 0.05"))) +
  theme(legend.position = "none")

### Output --------------------------
pT | pP


