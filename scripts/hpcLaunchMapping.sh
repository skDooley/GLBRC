#!/bin/bash

cd /mnt/research/ShadeLab/GLBRC 
batchn=0
for fname in scripts/hpc_scripts/filteredBamMaker/1*;
do
    ((batchn++))
    baseFileName=$(basename "$fname" .sb)
    #echo $baseFileName
    statsFile="stats/filtered/$baseFileName.tsv"
    if [ ! -f $statsFile ]; then
    	echo $batchn $statsFile $fname
    	sbatch --job-name=bowtie$batchn -A shadeash-colej --cpus-per-task=16 --mem=30G --output=logs/bowtie_4_2_$batchn.slurm.log --ntasks=1 --time=0:45:00 $fname
    fi
done