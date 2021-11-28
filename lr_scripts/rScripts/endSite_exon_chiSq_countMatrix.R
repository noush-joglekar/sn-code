#!/bin/R
# By Anoushka Joglekar 02.2021
# Modified 04.2021
# ChiSq test of association of end-site and exon 
# Modified 09.2021
# Making sure that only the necessary things are being tested

library(dplyr)
library(tidyr)
library(parallel)
library(data.table)
library(pbmcapply)

args <- commandArgs(trailing = TRUE)
# 1. matrix of counts
# 2. endSite name
# 3. cellType (or bulk)
# 4. numThreads
# 5. outDir
## config and or threshold

perCT <- function(countMatrix,args,outCT){
	print(paste0("Starting analysis in ",outCT))
	print(nrow(countMatrix %>% select(Exon) %>% distinct()))
	endSiteType <- args[2]
	nThreads <- as.integer(args[4])
	outDir <- args[5]

	eligible <- countMatrix %>% group_by(Exon) %>% 
		add_count() %>% filter(n>1) %>% 
		summarize(tot = sum(Inclusion) + sum(Exclusion))

	countDF <- countMatrix %>% filter(Exon %in% eligible$Exon)
	
	allIter = pbmclapply(eligible$Exon, function(e) getChiSq(e, countDF, endSiteType),
			mc.cores = nThreads,ignore.interactive = TRUE)

	counts = rbindlist(do.call( rbind, allIter)[,1]) %>% select(-c(donor,acceptor,es_s,es_e,strand))
	a = do.call( rbind, allIter)[,2]
	results = rbindlist(purrr::map(a, ~ purrr::compact(.)) %>% purrr::keep(~length(.) != 0))
	eV = rbindlist(do.call( rbind, allIter)[,3])
	write.table(eV, file.path(outDir,paste0(args[2],"_ExonCoord_",outCT,"_expValues")),
                        quote = FALSE, row.names = FALSE, sep = "\t")
	if(nrow(results) > 1){
		results$FDR = p.adjust(results$pval, method = "BY")
	
		write.table(results, file.path(outDir,paste0(args[2],"_ExonCoord_",outCT,"_results")),
			quote = FALSE, row.names = FALSE, sep = "\t")
	
		write.table(counts, file.path(outDir,paste0(args[2],"_ExonCoord_",outCT,"_counts")),
			quote = FALSE, row.names = FALSE, sep = "\t")
	
		sink(file.path(outDir,paste0("REPORT_",args[2],"_",outCT,".txt")), append=TRUE, split=TRUE)
		Sys.time()
		cat("Exon coordination of",args[2],"usage by",outCT,"\n")		
		cat("Total exon-endSite pairs tested:",length(results$FDR),"\n")
		cat("Number of", args[2], "with |deltaPSI| >= 0.1:",nrow(results %>% filter(abs(dPSI) >= 0.1)),"\n")
		cat("Number of significant exons by BY:",nrow(results %>% filter(FDR <= 0.05)),"\n")
		cat("Number of BY sig exon-endSite pairs with deltaPSI >= 0.1:",nrow(results %>% filter(FDR <= 0.05 & abs(dPSI) >= 0.1)),"\n")
		cat("Number of BY sig genes with deltaPSI >= 0.1:",nrow(results %>% filter(FDR <= 0.05 & abs(dPSI) >= 0.1) %>% select(Gene) %>% distinct()),"\n")
		cat('\n')
		Sys.time()
		sink()
		return()
	}
}

filterIneligibleES <- function(mat,endSiteType){
	if((endSiteType == "TSS" & mat$strand[1] == "+") || (endSiteType == "PolyA" & mat$strand[1] == "-")){
		mat = mat %>% filter(es_e <= min(donor)) %>% group_by_at(endSiteType) %>%
			slice_min(n = 1, order_by = "donor") %>% as.data.frame()
	} else if ((endSiteType == "TSS" & mat$strand[1] == "-") || (endSiteType == "PolyA" & mat$strand[1] == "+")){
		mat = mat %>% filter(es_s >= max(acceptor)) %>% group_by_at(endSiteType) %>%
			slice_max(n = 1, order_by = "acceptor") %>% as.data.frame()
	}
	if(nrow(mat) >= 2){
		return(mat)
	} else {
		return(NULL)}
}

getChiSq <- function(exon, countDF, endSiteType){
	mat <- countDF %>% filter(Exon == exon) 
	newMat <- filterIneligibleES(mat,endSiteType)
	if(!is.null(newMat)){
		mat = newMat %>% select(Inclusion,Exclusion)
		criterionStatus = chisq_criterion(mat)
		ev = criterionStatus[[2]]
		cons = criterionStatus[[3]]
		if (criterionStatus[[1]]=="True"){
			pVal <- chisq.test(mat)$p.value
			loR <- logOddsCalc(mat)
			pi = mat$Inclusion/sum(mat$Inclusion)
			pe = mat$Exclusion/sum(mat$Exclusion)
			mat$dPI = pi-pe
			if (sum(as.numeric(as.character(mat$dPI))[1:2]) == 0){
				ix = sample(c(1,2),1,0.5)
			} else {ix = which(abs(mat$dPI) == max(abs(mat$dPI)) )}
			dPSI = mat$dPI[ix]	
			subDF = newMat
			subDF$dPI = mat$dPI
			subDF$pval = pVal
			res = data.frame("Gene" = subDF$Gene[1],
					"Exon" = exon,
					"dPSI"= dPSI,
					"pval"= pVal,
					"logOddsRatio" = loR,
					"nrow" = nrow(mat))
			expVal = data.frame("Gene" = subDF$Gene[1],
                                	"Exon" = exon,
                                	"EV" = ev,
					"cons" = cons,
                               		"nrow" = nrow(mat))
			return(list(subDF,res,expVal))
		} else {
			res = NULL
			subDF = NULL
			m = countDF %>% filter(Exon == exon)
			expVal = data.frame("Gene" = m$Gene[1],
					"Exon" = exon,
					"EV" = ev,
					"cons" = cons,
					"nrow" = nrow(mat))
			return(list(subDF,res,expVal))
		}
	} else {return(NULL)}
}


chisq_criterion <- function(mat){	## Expected value in 80% (rounded to nearest integer) cells is atleast 5, and of all cells atleast 1
	ev <- min(rowSums(mat))*min(colSums(mat))/sum(mat)
	evMat <- sapply(colSums(mat), 
		function(c) rowSums(mat)*c/sum(mat))
	indEV <- length(which(evMat > 5))
	if(sum(mat)>=25 & any(apply(evMat,2,median) <5)){
		cs = "cons"
	} else {cs = "lowCounts"}
	threshold <- round(0.8*nrow(mat)*ncol(mat),0)
	if (ev >= 1 & indEV >= threshold){
		cs = "testable"
		return(list("True",ev, cs))
	 } else {return(list("False",ev, cs))}
}

logOddsCalc <- function(mat){	## Set 0's to 0.5 and explicitly calculate log odds ratio for 2x2 matrices
	if(nrow(mat)==2){
		mat[mat==0]=0.5
		lor = log2((mat[1,1]*mat[2,2])/(mat[2,1]*mat[1,2]))
		return(lor)
	} else {
		return("NA")}
}

outDir <- args[5]

if(!dir.exists(outDir)){dir.create(outDir)}

if(tolower(args[3]) == "celltype"){
	countMat <- as.data.frame(fread(args[1]))
	ct_split = split.data.frame(countMat,countMat$Celltype)
	l = lapply(names(ct_split),function(name) perCT(ct_split[[name]],args,name))
} else {
	countMat <- as.data.frame(fread(args[1]))
	perCT(countMat,args,"Bulk")
}
