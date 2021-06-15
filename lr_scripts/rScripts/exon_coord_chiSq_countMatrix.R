#!/bin/R
# By Anoushka Joglekar 04.2021
# ChiSq test for exon coordination

library(dplyr)
library(tidyr)
library(parallel)
library(data.table)
library(pbmcapply)

args <- commandArgs(trailing = TRUE)
# 1. matrix of counts
# 2. cellType (or pseudoBulk)
# 3. numThreads
# 4. outDir
## config and or threshold

perCT <- function(countMatrix,args,outCT){
	print(paste0("Starting analysis in ",outCT))
	print(nrow(countMatrix))
	nThreads <- as.integer(args[3])
	outDir <- args[4]

	allIter = pbmclapply(1:nrow(countMatrix), function(l) getChiSq(l, countMatrix),
			mc.cores = nThreads,ignore.interactive = TRUE)

	results = rbindlist(do.call( rbind, allIter)[,1])
	eV = rbindlist(do.call( rbind, allIter)[,2])

	write.table(eV, file.path(outDir,paste0("ExonCoord_",outCT,"_expValues")),
                        quote = FALSE, row.names = FALSE, sep = "\t")
	if(nrow(results) > 1){
		results$FDR = p.adjust(results$pval, method = "BY")
	
		write.table(results, file.path(outDir,paste0("ExonCoord_",outCT,"_results")),
			quote = FALSE, row.names = FALSE, sep = "\t")
	
		sink(file.path(outDir,paste0("REPORT_",outCT,".txt")), append=TRUE, split=TRUE)
		Sys.time()
		cat("Exon coordination of",outCT,"\n")		
		cat("Total exon pairs tested:",length(results$FDR),"\n")
		cat("Number of exons with |LOR| >= 1:",nrow(results %>% filter(abs(logOddsRatio) >= 1)),"\n")
		cat("Number of significant exons by BY:",nrow(results %>% filter(FDR <= 0.05)),"\n")
		cat("Number of BY sig exon pairs with |LOR| >= 1:",nrow(results %>% filter(FDR <= 0.05 & abs(logOddsRatio) >= 1)),"\n")
		cat("Number of BY sig genes with |LOR| >= 1:",
			nrow(results %>% filter(FDR <= 0.05 & abs(logOddsRatio) >= 1) %>% select(gene) %>% distinct()),"\n")
		cat('\n')
		Sys.time()
		sink()
		return()
	}
}


getChiSq <- function(l, countDF){
	line <- countDF[l,]
	mat <- matrix(as.numeric(line[,5:8]), 2, byrow = T)
	criterionStatus = chisq_criterion(mat)
	ev = criterionStatus[[2]]
	cons = criterionStatus[[3]]
	if (criterionStatus[[1]]=="True"){
		expVal <- line[,-c(5:8)]
		pVal <- chisq.test(mat)$p.value
		loR <- logOddsCalc(mat)
		line$pval <- pVal
		line$logOddsRatio <- loR
		expVal$EV <- ev
		expVal$cons <- cons
		return(list(line,expVal))
	} else {
		expVal <- line[,-c(5:8)]
		expVal$EV <- ev
		expVal$cons <- cons
		line = NULL
		return(list(line,expVal))
	}
}


chisq_criterion <- function(mat){	## Expected value in 80% (rounded to nearest integer) cells is atleast 5, and of all cells atleast 1
	ev <- min(rowSums(mat))*min(colSums(mat))/sum(mat)
	a <- which(sapply(colSums(mat),
                function(c) rowSums(mat)*c/sum(mat)) > 5)
	indEV <- length(a)
	a <- paste(a,collapse="_")
	allowedCB <- c("1_2","1_3","2_4","3_4")
	if (indEV==2 & a %in% allowedCB){	## simplistic: if row or column too small then becomes constitutive
		cs = "cons"
	} else {cs = "lowCounts"}
	threshold <- round(0.8*nrow(mat)*ncol(mat),0)
	if (ev >= 1 & indEV >= threshold){
		cs = "testable"
		return(list("True",ev,cs))
	 } else {return(list("False",ev,cs))}
}

logOddsCalc <- function(mat){	## Set 0's to 0.5 and explicitly calculate log odds ratio for 2x2 matrices
	if(nrow(mat)==2){
		mat[mat==0]=0.5
		lor = log2((mat[1,1]*mat[2,2])/(mat[2,1]*mat[1,2]))
		return(lor)
	} else {
		return("NA")}
}

outDir <- args[4]

if(!dir.exists(outDir)){dir.create(outDir)}

if(tolower(args[2]) == "celltype"){
	countMat <- as.data.frame(fread(args[1]))
	ct_split = split.data.frame(countMat,countMat$cellType)
	l = lapply(names(ct_split),function(name) perCT(ct_split[[name]],args,name))
} else {
	countMat <- as.data.frame(fread(args[1]))
	perCT(countMat,args,"pseudoBulk")
}
