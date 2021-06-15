#! /bin/bash
# By Anoushka Joglekar 02.2021
### Sort gene list and chunk into 1000 genes per file
### Break up the 10X Gene-BC-UMI and the ONT AllInfo files according to this chunking metric
### Run the python script multi-threaded on all the chunks in a for loop
### In the same for loop, run the postprocessing distance script

tx_gbu=$1;
allInfoFile=$2;
outName=$3;
pyDir=".";
shDir="../bashScripts";
awkScript=$shDir"/v2.0a_bc_finder.awk"
pyScript=$pyDir"/minEditDistance_post11merMatching.py"

echo "Preparing files for chunking";

cat $tx_gbu | awk '{print $2}' | sort -u > geneNames
split -l 1000 -d geneNames geneSet
rm geneNames

d=$((`ls geneSet* | wc -l` -1 ))

echo "Chunking files into 1000 gene bits";

for i in `eval echo {00..$d}` ; do awk -v i=$i -v aI=$allInfoFile \
BEGIN{comm="cat geneSet"i; while(comm|getline) {a[$1]=$1;} \
comm="zcat "aI; while(comm|getline) {split($2,gene,"."); \
if(gene[1] in a) {print }} }' > lrTech_${i} ; gzip lrTech_${i} & done

for i in `eval echo {00..$d}` ; do awk 'NR==FNR {a[$1]=$1;next} \
$2 in a {print }' geneSet${i} $tx_gbu > 10x_${i} & done

echo "Running python script in parallel"

for i in `eval echo {00..$d}`
do
echo "Processing chunk $i"

time awk -v unzipCommand=zcat -v geneSetFile=geneSet${i} \
-v bcUmiGeneFile=10x_${i} -v read2GeneFile=lrTech_${i}.gz \
-v fastqGZ=file_${i}.fastq.gz -v expectedBarcodePosFWD=50 \
-v expectedBarcodePosREV=35 -f $awkScript | gzip -c > combined_${i}.gz

sleep 5

python $pyScript --inFile combined_${i}.gz --outName summaryOut_${i} ## per 1000 gene block

done


echo "Concatenating and cleaning up"
zcat combined_*.gz | gzip > $outName"_potentialMatches.gz"
awk 'FNR==1 && NR!=1{next;}{print}' summaryOut_*_full.csv > $outName"_full.csv"
awk 'FNR==1 && NR!=1{next;}{print}' summaryOut_*_min.csv > $outName"_min.csv"
rm lrTech_* 10x_* geneSet* file_* summaryOut_* combined_* xyz.*

echo "Done!"
