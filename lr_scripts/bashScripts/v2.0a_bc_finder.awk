## quick and dirty script
## by Hagen Tilgner 2021.02

function comp(x){
    x=toupper(x) ;
    if(x=="A"){return("T");}
    if(x=="T"){return("A");} 
    if(x=="C"){return("G");} 
    if(x=="G"){return("C");} 
    if(x=="N"){return("N");}   
} 
function revcomp(x){
    y=""; 
    for(i=length(x);i>=1;i--){y=y""comp(substr(x,i,1));} 
    return(y);
}
function findElements(localSeq,localExpecPos,localK,localGene,kmerGenepairsLocal,barcodeUMIPairLength){
   
    # A. initializing
    foundKmers="NONE";
    foundBCs="NONE"
    readRegions2Check="NONE"
    foundPos="NONE";
    nfoundKmers=0;

    # B. polyT-localization and barcode localization
    from=localExpecPos-50;
    if(from<1){from=1;}
    to=localExpecPos+50;
    if(to>length(localSeq)-(localK-1)){to=length(localSeq)-(localK-1);}
    
    for(j=from;j<=to;j++){
	substring=substr(localSeq,j,localK);
	if(substring"\t"localGene in kmerGenepairsLocal){
	    nfoundKmers++;
	    if(foundKmers=="NONE"){
		foundPos=j;
		foundKmers=substring;
		for(k=1;k<=kmerGenepairsLocal[substring"\t"localGene];k++){
		    foundBCs=foundBCs","kmerGenepair2barcodeUMIPair[substring"\t"localGene"\t"k];
		    print substring"\t"localGene"\t"k"\t"kmerGenepair2barcodeUMIPairPosition[substring"\t"localGene"\t"k] > "xyz.1"
		    substringInRead=substr(localSeq,j + 1 - kmerGenepair2barcodeUMIPairPosition[substring"\t"localGene"\t"k], barcodeUMIPairLength);
		    readRegions2Check=readRegions2Check","substringInRead;
		    
		}
		foundBCs=substr(foundBCs,6,length(foundBCs)-5);
		readRegions2Check=substr(readRegions2Check,6,length(readRegions2Check)-5);
		if(length(substringInRead)!=28){
		    print length(substringInRead),substringInRead,from,to,localSeq > "abc.1"
		}
	    }
	    else{
		foundPos=foundPos","j
		foundKmers=foundKmers","substring;
		for(k=1;k<=kmerGenepairsLocal[substring"\t"localGene];k++){
		    foundBCs=foundBCs","kmerGenepair2barcodeUMIPair[substring"\t"localGene"\t"k];
		    print substring"\t"localGene"\t"k"\t"kmerGenepair2barcodeUMIPairPosition[substring"\t"localGene"\t"k] > "xyz.2"
		    substringInRead=substr(localSeq,j + 1 - kmerGenepair2barcodeUMIPairPosition[substring"\t"localGene"\t"k], barcodeUMIPairLength);
		    readRegions2Check=readRegions2Check","substringInRead;
		    if(length(substringInRead)!=28){
			print length(substringInRead),substringInRead,from,to,localSeq > "abc.1"
		    }
		}
	    }
	}
    }
    uniqueCounter=0; 
    foundBCsUniq="NONE";
    foundRegions2CheckUniq="NONE"
    if(foundBCs!="NONE"){
	n=split(foundBCs,localArr,",");
	n2=split(readRegions2Check,localArr2,",");
	if(n!=n2){
	    print "ERROR:n="n";n2="n2" forfoundBCs="foundBCs" and readRegions2Check="readRegions2Check > "/dev/stderr";
	    exit(0);
	}
	
	foundBCsUniq=localArr[1];
	foundRegions2CheckUniq=localArr2[1];
	uniqueCounter=1;
	for(k in seen){delete seen[k];}
	seen[localArr[l]"\t"localArr2[l]]=1;
	
	for(l=2;l<=n;l++){
	    if(!(localArr[l]"\t"localArr2[l] in seen)){
		seen[localArr[l]"\t"localArr2[l]]=1;
		foundBCsUniq=foundBCsUniq","localArr[l];
		foundRegions2CheckUniq=foundRegions2CheckUniq","localArr2[l]
		uniqueCounter++;
	    }
	}
    }
    
   
    return(uniqueCounter"\t"foundBCsUniq"\t"foundRegions2CheckUniq);
}

BEGIN{
    
    print "### executing v2.0a_bc_finder.awk" > "/dev/stderr";
    if(!unzipCommand){print "ERROR: no value for unzipCommand" > "/dev/stderr";exit(0);}
    if(!geneSetFile){print "ERROR: no value for geneSetFile" > "/dev/stderr";exit(0);}
    if(!read2GeneFile){print "ERROR: no value for read2GeneFile" > "/dev/stderr";exit(0);}
    if(!fastqGZ){print "ERROR: no value for fastqGZ" > "/dev/stderr";exit(0);}
    if(!bcUmiGeneFile){print "ERROR: no value for bcUmiGeneFile" > "/dev/stderr";exit(0);}
    #if(!outputFastqGZ){print "ERROR: no value for outputFastqGZ" > "/dev/stderr";exit(0);}
    if(!expectedBarcodePosFWD){print "ERROR: no value for expectedBarcodePosFWD" > "/dev/stderr";exit(0);}
    if(!expectedBarcodePosREV){print "ERROR: no value for expectedBarcodePosREV" > "/dev/stderr";exit(0);}
    #if(!outputGZString){print "ERROR: no value for outputGZString" > "/dev/stderr";exit(0);}

    outputGZString="gzip > "outputFastqGZ;

    print "# reading geneSetFile="geneSetFile > "/dev/stderr";
    comm="cat "geneSetFile;
    while(comm|getline){
	split($1,a,".");
	legalGenes[a[1]]=1;
    }

    print "# reading bcUmiGeneFile="bcUmiGeneFile > "/dev/stderr";
    comm="cat "bcUmiGeneFile;
    while(comm|getline){

	gene=$2;
	if(!(gene in legalGenes)){continue;}
	#print "got here" > "x2"
	gene2nMols[gene]++;
	barcodeUMIPairTrunk=toupper($4)""substr(toupper($5),1,6);
	barcodeUMIPair=toupper($4)""toupper($5);
	barcode=toupper($4)
	
	# treating each kmer (usually k=11)in the 22mer ()
	for(i=1;i<=length(barcodeUMIPairTrunk)-10;i++){
	    kmer=substr(barcodeUMIPairTrunk,i,11);
	    kmerGenepairs[kmer"\t"gene]++;
	    kmerGenepair2barcodeUMIPair[kmer"\t"gene"\t"kmerGenepairs[kmer"\t"gene]]=barcodeUMIPair;
	    kmerGenepair2barcodeUMIPairPosition[kmer"\t"gene"\t"kmerGenepairs[kmer"\t"gene]]=i;
	    
	}
    }   

#    for(k in kmerGenepairs){
#	print k > "y";
#    }
    
    print "# reading read2GeneFile="read2GeneFile > "/dev/stderr";
    comm="zcat "read2GeneFile;
    while(comm|getline){
	split($2,a,".");
	if(!(a[1] in legalGenes)){continue;}
	read2Gene[$1]=a[1];
    }

    print "# reading fastqGZ="fastqGZ > "/dev/stderr";
    comm=unzipCommand" "fastqGZ;
    while(comm|getline){
	lineNumber++;
	#if(lineNumber>4000000){
	#    break;
	#}
	if(lineNumber % 10000000 ==0){print "lineNumber="lineNumber > "/dev/stderr";}
	if(lineNumber % 4 == 1){readID=substr($1,2,length($1)-1);continue;}
	if(lineNumber % 4 == 2){seq=toupper($0);continue;}
	#print "got here with readID="readID > "z";
	if(lineNumber>=4 && lineNumber % 4 == 0 && length(seq)>=400 && readID in read2Gene){

	    
	    qualString=$0;
	    # I.finding sequences: the 125 value needs to be adjusted, so that we have enough sequence  
	    fwd=findElements(substr(seq,1,expectedBarcodePosFWD+100),expectedBarcodePosFWD,11,read2Gene[readID],kmerGenepairs,28);
	    #split(fwd,fwdArr,"___");
	    extend=expectedBarcodePosREV+100
	    rev=findElements(revcomp(substr(seq,length($0) - extend + 1,extend)),expectedBarcodePosREV,11,read2Gene[readID],kmerGenepairs,28);
	    #split(rev,revArr,"___");
	    
	    # II. printing
	    print readID"\t"read2Gene[readID]"\t"fwd"\t"rev
	}
    }
}
