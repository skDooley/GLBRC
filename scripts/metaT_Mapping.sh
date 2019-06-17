#!/bin/bash -login
#PBS -l walltime=50:00:00,nodes=01:ppn=8,mem=80gb
#PBS -q main
#PBS -M adina.chuang@gmail.com
#PBS -m abe
#PBS -A ged

THREADS=10
MY_PATH=/mnt/research/ShadeLab/GLBRC/
META_T=/mnt/research/ShadeLab/GLBRC/mapping/metaT
CONTIGS_FILE=$MY_PATH/annotatedContigs.fa

#Get the sample files
cd $META_T/unpaired
files=(*.gz)

#Go to the base dir
cd $MY_PATH

#for each sample
for f in "${files[@]}"; do 
	#Step 1. Make a batch script file for each sample
	cat scripts/hpc_scripts/header.sb >scripts/hpc_scripts/filteredBamMaker/$f.sb
		
	#Step 2. Separate the combined paired-end reads file into 2 separate files (PE1, PE2)
	echo "split-paired-reads.py --gzip -1 $META_T/paired/$f.pe1 -2 $META_T/paired/$f.pe2 $META_T/unpaired/$f" >>scripts/hpc_scripts/filteredBamMaker/$f.sb

	#Step 3. Trim adapters and QC reads
	echo "java -jar ~/bin/trimAdapters/trimmomatic-0.38.jar PE -phred33 -summary $META_T/logs/trimStats/$f.log -threads 10 $META_T/paired/$f.pe1 $META_T/paired/$f.pe2 $META_T/trimmed/$f.pe1 $META_T/trimmed/$f.se1 $META_T/trimmed/$f.pe2 $META_T/trimmed/$f.se2 ILLUMINACLIP:/mnt/home/dooleys1/bin/trimAdapters/NexteraPE-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" >>scripts/hpc_scripts/filteredBamMaker/$f.sb

	#Step 4. Align the trimmed reads to the metagenomic assembly
	echo "bowtie2 -p 10 -x $CONTIGS_FILE -1 $META_T/trimmed/$f.pe1 -2 $META_T/trimmed/$f.pe2 -S $META_T/sams/$f.sam >$META_T/flagstats/$f.stat 2>&1 " >> scripts/hpc_scripts/filteredBamMaker/$f.sb
	
	#Step 5. Compress and filter the sam file to make a bam (filter = remove improperly paired reads)
	echo "samtools view --threads 10 -f 0x2 -b -S $META_T/sams/$f.sam -t $$CONTIGS_FILE >$META_T/bams/$f.bam" >> scripts/hpc_scripts/filteredBamMaker/$f.sb
	
	#Step 6. Remove the sam file to save space
	echo "rm $META_T/sams/$f.sam">> scripts/hpc_scripts/filteredBamMaker/$f.sb

	#Step 7. Sort the reads
	echo "samtools sort -o $META_T/bams/$f.bam $META_T/bams/$f.sorted.bam" >> scripts/hpc_scripts/filteredBamMaker/$f.sb
	echo "rm $META_T/bams/$f.bam" >> scripts/hpc_scripts/filteredBamMaker/$f.sb

	#Step 8. Index the reads
	echo "samtools index -@ 10 $META_T/bams/$f.sorted.bam" >> scripts/hpc_scripts/filteredBamMaker/$f.sb

	#Step 9. Get the read counts along the contigs
	echo "samtools idxstats $META_T/bams/$f.sorted.bam >$META_T/stats/$f.tsv" >> scripts/hpc_scripts/filteredBamMaker/$f.sb

	#Step 10. Launch the script on the hpc
	sbatch scripts/hpc_scripts/filteredBamMaker/$f.sb
	
done


# CONT_PATH=/mnt/research/ShadeLab/GLBRC/assemblies
# CONTIGS_FILE=assemblyGeneSeqs.fa
# CONTIGS_FILE_DB_DIR=/mnt/research/ShadeLab/GLBRC/assemblies/bowtieDB
# PAIRED_DIR=/mnt/research/ShadeLab/GLBRC/jgi_transfer/Raw_Data/metag-pipeline/paired-reads
# #declare -a files=("11815.1.220223.CTTGT.fastq.gz" "11815.1.220223.AGTCA.fastq.gz" "11814.1.220215.TAGCT.fastq.gz" "11814.1.220215.GGCTA.fastq.gz" "11861.3.224524.CGTACG.fastq.gz")
# cd $PAIRED_DIR
#module load bedtools
#python coverage-bed-reference.py $CONTIGS_FILE > $CONTIGS_FILE.bed

#for x in mapping-data/*sorted; do echo "bamToBed -i $f > $f.bed"; done > bamtobed.sh
#cat bamtobed.sh | $PAR_PATH/parallel
#for x in mapping-data/*bed; do echo "coverageBed -a $CONTIGS_FILE.bed -b $f -d > $f.bed2"; done > coveragebed.sh
#cat coveragebed.sh | $PAR_PATH/parallel
#for x in mapping-data/*bed2; do echo "python bedcoverage-to-coverage.py $f > $f.counts"; done > bedcoveragefinal.sh
#cat bedcoveragefinal.sh | $PAR_PATH/parallel
