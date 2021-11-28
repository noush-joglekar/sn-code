#!/bin/R
# By Anoushka Joglekar 02.2021
# Modified 04.2021 to add reads per gene argument
## Generate inclusion-exclusion count matrix of exons associated with end-sites
## This is in order to evaluate coordination between alternative exons and transcription
## start and end sites

## Modified 09.2021 - Keeping track of preceeding and succeeding exon boundaries
## Modified again 09.2021 - Making sure read ends after the acceptor

# Setup -----
library(dplyr)
library(tidyr)
library(data.table)
library(parallel)
library(pbmcapply)
library(stringr)

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

precursorDF <- start_end %>%
	separate_rows(Exon,sep = ";%;") %>% filter(Exon != "") %>%
	group_by(Read) %>% mutate(PrevExon = lag(Exon), NextExon = lead(Exon)) %>% 
	na.omit() %>% ungroup() %>% rowwise() %>%		## automatically selects internal exons only
	mutate(donor=unlist(strsplit(PrevExon,"_"))[3],
		acceptor=unlist(strsplit(NextExon,"_"))[2] ) %>% 
	select(-c(PrevExon,NextExon)) %>% 
	as.data.frame()

longEnoughReads <- precursorDF %>% group_by(.dots = groupVar[3], Gene, Exon) %>%
	mutate(md=min(donor), ma = max(acceptor)) %>% filter(start < md & end > ma) %>% ungroup()

endSite_exon <- longEnoughReads %>% select(all_of(endSite),Exon,donor,acceptor) %>% distinct() %>%
	separate(endSite, into = c("chr","es_s","es_e","strand"), sep ="_", remove=FALSE) %>%
	mutate_at(c("es_s","es_e","donor","acceptor"),as.numeric) %>%
	mutate(status = case_when(endSite == "TSS" & strand == "-" & es_s >= acceptor ~ "correct",	
			endSite == "TSS" & strand == "+" & es_e <= donor ~ "correct",	## 5-->3 "+" TSS - Donor - Exon - Acceptor - PolyA
			endSite == "PolyA" & strand == "-" & donor >= es_e ~ "correct",	## 3-->5 "-" PolyA - Donor - Exon - Acceptor - TSS
			endSite	== "PolyA" & strand == "+" & acceptor <= es_s ~ "correct",
			TRUE ~ "incorrect")) %>%
	filter(status == "correct") %>%
	select(all_of(endSite),Exon, donor, acceptor,es_s,es_e,strand)

## have to add something in here for exclusion counts where the read has to start before the min donor and end after 
## the min acceptor

#X = dim(a %>% group_by(Exon) %>% filter(Exon == exon) %>% 
#mutate_at(c("start","end","donor","acceptor"),as.numeric) %>% 
#mutate(md = min(donor), ma = max(acceptor)) %>% filter(start <= md & end >= ma) %>% as.data.frame())

filterIneligiblePairs <- function(eS_e,df){
	if((eS_e$strand[1] == "+" & endSite == "TSS") || (eS_e$strand[1] == "-" & endSite == "PolyA") ){
		conservativePairs <- left_join(df, eS_e, by = c(endSite,'Exon')) %>% 
			group_by(Exon) %>% slice_min(n = 1, order_by = donor) %>%
			filter(es_e <= donor) %>% as.data.frame()
	} else if ((eS_e$strand[1] == "-" & endSite == "TSS") || (eS_e$strand[1] == "+" & endSite == "PolyA")){
		conservativePairs <- left_join(df, eS_e, by = c(endSite,'Exon')) %>%
                       	group_by(Exon) %>% slice_max(n = 1, order_by = acceptor) %>%
                       	filter(es_s >= acceptor) %>% as.data.frame()
	}
	return(conservativePairs)
}


getCounts <- function(es, start_end, endSite_exon, longEnoughReads, args){
	data_subset <- start_end %>% filter(.[[endSite]] == es & Read %in% longEnoughReads$Read) %>% 
		select(-Read) %>% group_by_all() %>% summarize(n=n(),.groups='keep') %>% 
		mutate_at(c("start","end"),as.numeric) %>% as.data.frame()
	exonList <- as.vector(endSite_exon %>% filter(.[[endSite]] == es) %>% 
		select(Exon) %>% as.matrix())
	eS_e <- endSite_exon %>% filter(.[[endSite]] == es)
	dd = list()
	for(exon in exonList){
		inCounts <- data_subset %>% filter(grepl(exon,Exon)) %>%
			group_by(.dots = groupVar[3]) %>% summarize(Inclusion = sum(n),.groups='keep') %>%
			mutate(Exon = exon)
		exCounts <- data_subset %>% filter(!grepl(exon,Exon)) %>%
			group_by(.dots = groupVar[3]) %>% mutate(exon = exon) %>%
			separate(exon, into = c("chr_e","s","e","strand_e"), sep ="_") %>%
			mutate_at(c("start","end","s","e"),as.numeric) %>%
			filter(start <= s & end >= e) %>% summarize(Exclusion = sum(n),.groups='keep') %>%
			mutate(Exon = exon)
		if(args[3] == "Celltype"){
			dd[[exon]] <- full_join(inCounts,exCounts,by = c("Celltype", "Exon")) %>% replace(is.na(.),0)
		} else {
			dd[[exon]] <- full_join(inCounts,exCounts,by = groupVar[-c(1:2)]) %>% replace(is.na(.),0)
		}
	}
	df = plyr::ldply(unname(dd), data.frame)
	df$Gene <- data_subset$Gene[1]
	df[,endSite] <- data_subset[1,endSite]
	eligPairs <- filterIneligiblePairs(eS_e,df)
	return(eligPairs)
}


bigList <- pbmclapply(unique(unname(unlist(endSite_exon[,endSite]))),function(es) getCounts(es,start_end,endSite_exon,longEnoughReads,args),
		mc.cores = nThreads,ignore.interactive = TRUE)

print("Converting to dataframe")
bigDF <- as.data.frame(data.table::rbindlist(bigList))


print("Writing output to file")
data.table::fwrite(bigDF, file = paste0("inclusionExclusionCounts_",endSite,"by_",args[3],"_",rpg,"readsPerGene",".gz"),
	sep = "\t", quote = FALSE, row.names = FALSE, compress="gzip")
 
