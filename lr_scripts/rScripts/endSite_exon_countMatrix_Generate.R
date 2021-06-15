#!/bin/R
# By Anoushka Joglekar 02.2021
# Modified 04.2021 to add reads per gene argument
## Generate inclusion-exclusion count matrix of exons associated with end-sites
## This is in order to evaluate coordination between alternative exons and transcription
## start and end sites

# Setup -----
library(dplyr)
library(tidyr)
library(data.table)
library(parallel)
library(pbmcapply)

## Arguments -----
# 1. All-Info file
# 2. EndSite
# 3. Grouping - Celltype or Bulk
# 4. Number of threads
# 5. Endsite assigned reads per gene

args <- commandArgs(trailing=TRUE)

print("Reading in file")

allInfo <- fread(args[1])
colnames(allInfo) <- c("Read","Gene","Celltype","Barcode","UMI",
		"Stretch","TSS","PolyA","Exon","Status","NumIntrons")

endSite <- args[2]

if(args[3] == "Celltype"){
	groupVar = c("Read","Gene","Celltype",endSite,"Exon")
} else {groupVar = c("Read","Gene",endSite,"Exon")}

nThreads <- as.integer(args[4])
rpg <- as.integer(args[5])

print("Filtering lowly expressed genes and organizing data")

complete_ES = allInfo %>% filter(.[[endSite]] != paste0("No",endSite)) %>% select(all_of(groupVar)) %>%
	group_by(Gene) %>% add_count() %>% filter(n >= rpg) %>% select(-n) %>% as.data.frame()

start_end = complete_ES %>% rowwise() %>% 
	mutate(start = unlist(strsplit(Exon,"_"))[2], end = rev(unlist(strsplit(Exon,"_")))[2]) %>% 
	as.data.frame()

endSite_exon <- start_end %>%
        separate_rows(Exon,sep = ";%;") %>% filter(Exon != "") %>% 
	group_by(Read) %>% add_count() %>% filter(n>=3) %>% slice(2:(n()-1)) %>%  ## internal exons only
	as.data.frame() %>% 
	select(all_of(endSite),Exon) %>% distinct() %>%
	separate(endSite, into = c("chr","es_s","es_e","strand"), sep ="_", remove=FALSE) %>%
	separate(Exon, into = c("chr_e","s","e","strand_e"), sep ="_", remove=FALSE) %>%
	mutate_at(c("es_s","es_e","s","e"),as.numeric) %>%
	mutate(status = case_when(endSite == "TSS" & strand == "-" & es_e >= e ~ "correct",	## CAN CHANGE TO OVERLAP WINDOW
			endSite == "TSS" & strand == "+" & es_s <= s ~ "correct",	## 5-->3 "+" TSS - Exon - PolyA
			endSite == "PolyA" & strand == "-" & es_s >= s ~ "correct",	## 3-->5 "-" PolyA - Exon - TSS
			endSite	== "PolyA" & strand == "+" & e <= es_e ~ "correct",
			TRUE ~ "incorrect")) %>%
	filter(status == "correct") %>%
	select(all_of(endSite),Exon)

getCounts <- function(es, start_end, endSite_exon){
	data_subset = start_end %>% filter(.[[endSite]] == es) %>% 
		select(-Read) %>% group_by_all() %>% summarize(n=n(),.groups='keep') %>% 
		mutate_at(c("start","end"),as.numeric) %>% as.data.frame()
	exonList = as.vector(endSite_exon %>% filter(.[[endSite]] == es) %>% 
		select(Exon) %>% as.matrix())
	dd = list()
	for(exon in exonList){
		d = list()
		ix = 1
		for(r in 1:nrow(data_subset)){
			l = NULL
			if(grepl(exon,data_subset[r,]$Exon)){	## If exon present add count to inclusion
				l = c(as.matrix(data_subset[r,groupVar[-c(1,length(groupVar))]]),exon,data_subset$n[r],'0')
			} else if (!grepl(exon,data_subset[r,]$Exon)){		## Else check if exon chain spans exon
				s = as.integer(unlist(strsplit(exon,"_"))[2])
				e = as.integer(unlist(strsplit(exon,"_"))[3])
				if(data_subset[r,]$start <= s & data_subset[r,]$end >= e){	## If it does, add count to exclusion
					l = c(as.matrix(data_subset[r,groupVar[-c(1,length(groupVar))]]),exon,'0',data_subset$n[r])} 
			}
			if(!is.null(l)){
	                	d[[ix]] = l
				ix = ix + 1}
		}
	dd[[exon]] = do.call("rbind",d)
	}
	df = plyr::ldply(unname(dd), data.frame)
	colnames(df) = c(groupVar[-1],"Inclusion","Exclusion")
	df <- df %>% mutate_at(c("Inclusion","Exclusion"),as.vector) %>% 
		mutate_at(c("Inclusion","Exclusion"),as.numeric) %>%
		group_by_at(groupVar[-1]) %>% 
		summarize(Inclusion = sum(Inclusion), Exclusion = sum(Exclusion),.groups='keep')
        return(df)
}

print("Extracting counts, might take a while")
bigList <- pbmclapply(unique(endSite_exon[,endSite]),function(es) getCounts(es,start_end,endSite_exon),
		mc.cores = nThreads,ignore.interactive = TRUE)

print("Converting to dataframe and writing output")
bigDF <- as.data.frame(data.table::rbindlist(bigList))
data.table::fwrite(bigDF, file = paste0("inclusionExclusionCounts_",endSite,"by_",args[3],"_",rpg,"readsPerGene",".gz"),
	sep = "\t", quote = FALSE, row.names = FALSE, compress="gzip")
