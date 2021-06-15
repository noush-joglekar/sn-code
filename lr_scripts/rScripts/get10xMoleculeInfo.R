#! /bin/R

# By Anoushka Joglekar 01.2021 to extract Gene-Barcode-UMI pairs from 10x data
## This is a precursor to script that will look at ONT read assignments and gauge
## whether there has been molecule inflation because of indels in barcodes and/or 
## UMIs

# Requires R version >= 4

library(DropletUtils)
library(dplyr)
library(R.utils)
library(plyr)

args <- commandArgs(trailing=TRUE)
bcFile <- args[1]
h5File <- args[2]
outName <- args[3]

bc <- read.table(args[1])

## read in 10x HDF5 file
mol.info <- read10xMolInfo(h5File, get.umi=TRUE, get.gene = TRUE, 
		get.cell = TRUE, barcode.length = 16)

molecule_df <- mol.info$data %>% as.data.frame() %>% 
	filter(cell %in% bc$V1)

molecule_df <- molecule_df %>% 
	mutate(Gene = mol.info$genes[gene])

## convert from 2-bit encoding to nucleotide sequences
twoBit_toSeq <- function(deci,uLen){
	bin_str <- intToBin(deci)
	len <- length(strsplit(bin_str,"")[[1]])
	leading <- 2*uLen - len
	if(leading > 0){
		padded <- paste(strrep(0,leading), bin_str, sep = "")
	} else {padded = bin_str}
	split_bin <- sapply(seq(1, 2*uLen, 2),
		function(ix) substr(padded,ix,ix+1))
	mapped_toSeq <- mapvalues(split_bin, from = c("00","01","10","11"), 
				to = c("A","C","G","T"),warn_missing = FALSE)
	mapped_str <- paste(mapped_toSeq, collapse ="")
	return(mapped_str)
}


par_umis <- unlist(mclapply(molecule_df$umi, function(umi) twoBit_toSeq(umi,12),mc.cores=12))

molecule_df$UMI <- par_umis

## preserving structure reqd for downstream processing
molecule_df <- molecule_df[,c(4, 6,2, 1, 7)]

write.table(molecule_df, paste0(outName,".csv"), 
	quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
