#/bin/R
# By Simon Hardwick, 2021
# Figure 3 for SnISOr-Seq paper

# Setup --------------------

library("dplyr")
library("ggplot2")
library("tidyr")
library("tidyverse")
library("pheatmap")
library("reshape2")

# Figure 3a generated with Adobe Illustrator using toy example

# Figure 3b --------------------

ONT_altExons <- read.csv('data/all.altExons.matrix_mod.tab', header = TRUE, sep = "\t")

ONT_altExons_atLeast2cellTypes <- ONT_altExons %>%
  rename(exon_id_noStrand = X.exonID__) %>%
  melt(id.var = "exon_id_noStrand", variable.name = "cellType", value.name = "psi") %>%
  filter(psi != "NA") %>%
  group_by(exon_id_noStrand) %>%
  summarise(num_types = n(), max_psi = max(psi), min_psi = min(psi)) %>%
  filter(num_types >= 2) %>%
  mutate(maxMinusMin = max_psi - min_psi) %>%
  separate(col = exon_id_noStrand, into = c("chr", "start", "end"), sep = "_", remove = FALSE) %>%
  mutate(length = as.numeric(end) - as.numeric(start) + 1)

ONT_altExons_atLeast2cellTypes$category <- ifelse(ONT_altExons_atLeast2cellTypes$maxMinusMin<=0.25, "maxMinusMin<=0.25",
                                                  ifelse(ONT_altExons_atLeast2cellTypes$maxMinusMin>0.25 & ONT_altExons_atLeast2cellTypes$maxMinusMin<=0.75, "0.25<maxMinusMin<=0.75", 
                                                         "maxMinusMin>0.75"))

ggplot(ONT_altExons_atLeast2cellTypes, aes(x=length, fill=category)) + geom_density(alpha=0.4) + xlim(c(0,250)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Figure 3c --------------------

ASD_asExons <- read.csv('data/ASD_asExons.txt', header = FALSE)

schizo_asExons <- read.csv('data/Schizophrenia_asExons.txt', header = FALSE)

ALS_asExons <- read.csv('data/ALS_asExons.txt', header = FALSE)

ASD_yes <- ONT_altExons %>%
  rename(exon_id_noStrand = X.exonID__) %>%
  filter(exon_id_noStrand %in% ASD_asExons$V1) %>%
  select(exon_id_noStrand,ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes) %>%
  filter(ExcitatoryNeurons != "NA") %>%
  filter(InhibitoryNeurons != "NA") %>%
  filter(Astrocytes != "NA") %>%
  filter(Oligodendrocytes != "NA") %>%
  mutate(max = pmax(ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes)) %>%
  mutate(min = pmin(ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes)) %>%
  mutate(maxMinusMin = max - min) %>%
  mutate(disease = "ASD") %>%
  mutate(type = "yes")

ASD_no <- ONT_altExons %>%
  rename(exon_id_noStrand = X.exonID__) %>%
  filter(!exon_id_noStrand %in% ASD_asExons$V1) %>%
  select(exon_id_noStrand,ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes) %>%
  filter(ExcitatoryNeurons != "NA") %>%
  filter(InhibitoryNeurons != "NA") %>%
  filter(Astrocytes != "NA") %>%
  filter(Oligodendrocytes != "NA") %>%
  mutate(max = pmax(ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes)) %>%
  mutate(min = pmin(ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes)) %>%
  mutate(maxMinusMin = max - min) %>%
  mutate(disease = "ASD") %>%
  mutate(type = "no")

Schizo_yes <- ONT_altExons %>%
  rename(exon_id_noStrand = X.exonID__) %>%
  filter(exon_id_noStrand %in% schizo_asExons$V1) %>%
  select(exon_id_noStrand,ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes) %>%
  filter(ExcitatoryNeurons != "NA") %>%
  filter(InhibitoryNeurons != "NA") %>%
  filter(Astrocytes != "NA") %>%
  filter(Oligodendrocytes != "NA") %>%
  mutate(max = pmax(ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes)) %>%
  mutate(min = pmin(ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes)) %>%
  mutate(maxMinusMin = max - min) %>%
  mutate(disease = "Schizophrenia") %>%
  mutate(type = "yes")

Schizo_no <- ONT_altExons %>%
  rename(exon_id_noStrand = X.exonID__) %>%
  filter(!exon_id_noStrand %in% schizo_asExons$V1) %>%
  select(exon_id_noStrand,ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes) %>%
  filter(ExcitatoryNeurons != "NA") %>%
  filter(InhibitoryNeurons != "NA") %>%
  filter(Astrocytes != "NA") %>%
  filter(Oligodendrocytes != "NA") %>%
  mutate(max = pmax(ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes)) %>%
  mutate(min = pmin(ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes)) %>%
  mutate(maxMinusMin = max - min) %>%
  mutate(disease = "Schizophrenia") %>%
  mutate(type = "no")

ALS_yes <- ONT_altExons %>%
  rename(exon_id_noStrand = X.exonID__) %>%
  filter(exon_id_noStrand %in% ALS_asExons$V1) %>%
  select(exon_id_noStrand,ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes) %>%
  filter(ExcitatoryNeurons != "NA") %>%
  filter(InhibitoryNeurons != "NA") %>%
  filter(Astrocytes != "NA") %>%
  filter(Oligodendrocytes != "NA") %>%
  mutate(max = pmax(ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes)) %>%
  mutate(min = pmin(ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes)) %>%
  mutate(maxMinusMin = max - min) %>%
  mutate(disease = "ALS") %>%
  mutate(type = "yes")

ALS_no <- ONT_altExons %>%
  rename(exon_id_noStrand = X.exonID__) %>%
  filter(!exon_id_noStrand %in% ALS_asExons$V1) %>%
  select(exon_id_noStrand,ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes) %>%
  filter(ExcitatoryNeurons != "NA") %>%
  filter(InhibitoryNeurons != "NA") %>%
  filter(Astrocytes != "NA") %>%
  filter(Oligodendrocytes != "NA") %>%
  mutate(max = pmax(ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes)) %>%
  mutate(min = pmin(ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes)) %>%
  mutate(maxMinusMin = max - min) %>%
  mutate(disease = "ALS") %>%
  mutate(type = "no")

neuroDisease_combined <- rbind(ASD_yes,ASD_no,ALS_yes,ALS_no,Schizo_yes,Schizo_no)

ggplot(neuroDisease_combined, aes(x=factor(disease, levels = c("Schizophrenia","ALS","ASD")), y=maxMinusMin, 
                                  fill=factor(type, levels = c("yes","no")))) + geom_boxplot(notch = TRUE, outlier.shape = NA) + 
  ggtitle("NeuroDisease exon variability") + 
  theme(aspect.ratio=1, axis.line = element_line(colour = "black"))

# Figure 3d --------------------

ONT_altExons_ASD_4mainTypes <- ONT_altExons %>%
  rename(exon_id_noStrand = X.exonID__) %>%
  filter(exon_id_noStrand %in% ASD_asExons$V1) %>%
  select(exon_id_noStrand,ExcitatoryNeurons,InhibitoryNeurons,Astrocytes,Oligodendrocytes) %>%
  filter(ExcitatoryNeurons != "NA") %>%
  filter(InhibitoryNeurons != "NA") %>%
  filter(Astrocytes != "NA") %>%
  filter(Oligodendrocytes != "NA") %>%
  column_to_rownames(var = "exon_id_noStrand") %>%
  as.matrix

ONT_altExons_micro_4mainTypes_ids <- ONT_altExons %>%
  rename(exon_id_noStrand = X.exonID__) %>%
  filter(exon_id_noStrand %in% ASD_asExons$V1) %>%
  filter(ExcitatoryNeurons != "NA") %>%
  filter(InhibitoryNeurons != "NA") %>%
  filter(Astrocytes != "NA") %>%
  filter(Oligodendrocytes != "NA") %>%
  separate(col = exon_id_noStrand, into = c("chr", "start", "end"), sep = "_", remove = FALSE) %>%
  mutate(length = as.numeric(end) - as.numeric(start) + 1) %>%
  filter(length <=27) %>%
  select(exon_id_noStrand) %>%
  mutate(exonType = "microExon")

ONT_altExons_nonMicro_4mainTypes_ids <- ONT_altExons %>%
  rename(exon_id_noStrand = X.exonID__) %>%
  filter(exon_id_noStrand %in% ASD_asExons$V1) %>%
  filter(ExcitatoryNeurons != "NA") %>%
  filter(InhibitoryNeurons != "NA") %>%
  filter(Astrocytes != "NA") %>%
  filter(Oligodendrocytes != "NA") %>%
  separate(col = exon_id_noStrand, into = c("chr", "start", "end"), sep = "_", remove = FALSE) %>%
  mutate(length = as.numeric(end) - as.numeric(start) + 1) %>%
  filter(length >27) %>%
  select(exon_id_noStrand) %>%
  mutate(exonType = "nonMicroExon")

ONT_altExons_ASD_microNonMicro_combined <- rbind(ONT_altExons_micro_4mainTypes_ids,ONT_altExons_nonMicro_4mainTypes_ids)

ONT_altExons_ASD_microNonMicro_combined <- ONT_altExons_ASD_microNonMicro_combined %>%
  column_to_rownames(var = "exon_id_noStrand")

pheatmap(ONT_altExons_ASD_4mainTypes, show_rownames = FALSE, annotation_row = ONT_altExons_ASD_microNonMicro_combined)

# Figure 3e --------------------

gencode_exons <- read.csv('data/gencode.v34.annotation_exonsUnique_noStrand.txt', header = FALSE)

ONT_altExons_bulkPSI <- read.csv('data/all.altExons_withGeneID_bulkPSI.tab', header = FALSE, sep = "\t")

ONT_altExons_bulkPSI <- ONT_altExons_bulkPSI %>%
  rename(exon_id_noStrand = V1, gene_id = V2, bulkPSI = V3)

ONT_altExons_novel_atLeast2cellTypes_maxMinusMin_bulkPSI <- ONT_altExons %>%
  rename(exon_id_noStrand = X.exonID__) %>%
  filter(!exon_id_noStrand %in% gencode_exons$V1) %>%
  melt(id.var = "exon_id_noStrand", variable.name = "cellType", value.name = "psi") %>%
  filter(psi != "NA") %>%
  group_by(exon_id_noStrand) %>%
  summarise(num_types = n(), max_psi = max(psi), min_psi = min(psi)) %>%
  mutate(maxMinusMin = max_psi - min_psi) %>%
  filter(num_types >= 2) %>%
  left_join(ONT_altExons_bulkPSI, by = "exon_id_noStrand") %>%
  select(exon_id_noStrand, maxMinusMin, bulkPSI)

ggplot(ONT_altExons_novel_atLeast2cellTypes_maxMinusMin_bulkPSI, aes(x=bulkPSI, y=maxMinusMin)) + geom_point() + geom_smooth(method = "loess") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Figure 3f generated using genome browser
