#! /bin/bash
# By Anoushka Joglekar 02.2021
## Figured out that the efficiency of the edit distance code for 3
## separate packages (Levenshtein, ssw, and edit.distance) is abysmal
## when one tries to parallelize in python using joblib

## Instead, here I will be implementing a bash parallelization in a for-loop
### Sort gene list and chunk into number of threads
### Break up the 10X and the PB/ONT allInfo files according to this chunking metric
### Run the python script 1-threaded on all the chunks

tx_gbu=$1;
allInfoFile=$2;
numThreads=$3;
outName=$4;
pyScript="../pythonScripts/moleculeInflation_10x_forParallel.py";

echo "Preparing files for chunking";

cat $tx_gbu | awk '{print $2}' | sort -u > geneNames
split -n $numThreads -d geneNames geneSet
rm geneNames

d=`echo $numThreads-1 | bc -l`

echo "Chunking files across $numThreads threads";

for i in `eval echo {00..$d}` ; do awk 'NR==FNR {a[$1]=$1;next} \
{split($2,gene,"."); if(gene[1] in a) {print }}' geneSet${i} \
$allInfoFile > lrTech_${i} & done

for i in `eval echo {00..$d}` ; do awk 'NR==FNR {a[$1]=$1;next} \
$2 in a {print }' geneSet${i} $tx_gbu > 10x_${i} & done

sleep 1
echo "Running python script in parallel"

for i in `eval echo {00..$d}` ; do \
python $pyScript --ont lrTech_${i} --tx 10x_${i} \
--outName outFile_${i} --threads 1 & done

wait

iter=0;
for i in `eval echo {00..$d}`
do
if [ -f outFile_${i}.csv ]; then
((iter=iter+1))
fi
done

if [ $iter == $numThreads ]; then
echo "Concatenating and cleaning up"
awk 'FNR==1 && NR!=1{next;}{print}' outFile_* > $outName
rm lrTech_* 10x_* geneSet* outFile_*
fi

echo "Done!"
