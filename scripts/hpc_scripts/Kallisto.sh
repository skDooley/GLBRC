#!/bin/bash
GLBRC=/mnt/research/ShadeLab/GLBRC
READ_TYPE=metaT
META_DIR=$GLBRC/mapping/$READ_TYPE
THREADS=20
MEM="40G"

# ASSEMBLY=$SCRATCH/bowtieDB/GLBRC_Metagenome/final.contigs.fa 
ASSEMBLY=$SCRATCH/bowtieDB/MAG_Assembly.fa #Assembly to be indexed
COUNTS=$META_DIR/kallistoCounts
INDEX_DIR=$SCRATCH/bowtieDB #Path of index file kallisto will create
INDEX_NAME="MAGsIndex.kai"
FASTQ_DIR=$SCRATCH/$READ_TYPE/cleaned_fastqs
SCRATCH_OUT=$SCRATCH/$READ_TYPE/kallistoOut
SCRIPTS_DIR=$GLBRC/scripts/hpc_scripts/kallisto
HEADER=$GLBRC/scripts/hpc_scripts/MappingHeader.sb
R1_EXT="R1.fastq" #The extension type for reads 1 files. Eg Sample1.R1.fastq
R2_EXT="R2.fastq" #The extension type for reads 2 files. Eg Sample1.R2.fastq  


cd $FASTQ_DIR
R1S=($(ls *.$R1_EXT)) #Get a list of the pe1/r1 reads
R2S=($(ls *.$R2_EXT)) #Get a list of the pe2/r2 reads
nsamples=${#R1S[@]} #Number of samples to process

#Double check that we have the same number of R1/R2 files
if [ $nsamples -ne ${#R2S[@]} ]; then
    echo "Different Number of R1s ($nsamples) and R2s (${#R2S[@]})"
    exit -1
fi

#Create the sbatch script
SAMPLE_STRING=""
for ((i=0; i<nsamples; i++)); do
	sample=${R1S[$i]}
	sample=${sample:0:-9}

	if [ -f "$COUNTS/$sample.mags.krona.kegg.minpath.tab" ]; then
		continue;
	fi

	mkdir -p $SCRATCH_OUT/$sample
	cat $HEADER >$SCRIPTS_DIR/$sample.sb
	echo "" >>$SCRIPTS_DIR/$sample.sb
	echo "cd $FASTQ_DIR" >>$SCRIPTS_DIR/$sample.sb
	echo "mkdir -p $SCRATCH_OUT/$sample" >>$SCRIPTS_DIR/$sample.sb
	echo "kallisto quant --threads $THREADS --pseudobam -i $INDEX_DIR/$INDEX_NAME -o $SCRATCH_OUT/$sample ${R1S[$i]} ${R2S[$i]} 2>&1" >>$SCRIPTS_DIR/$sample.sb
	echo "cd $SCRATCH_OUT/$sample" >>$SCRIPTS_DIR/$sample.sb
	echo -en "echo -en \"Sorting\"\n\n" >>$SCRIPTS_DIR/$sample.sb
	echo "samtools sort --threads $THREADS -m 2G -t -n pseudoalignments.bam -o pseudoalignments_sorted.bam" >>$SCRIPTS_DIR/$sample.sb
	echo "samtools index -@ $THREADS pseudoalignments_sorted.bam" >>$SCRIPTS_DIR/$sample.sb
		echo -en "
retval=\$?
if [ \$retval -ne 0 ]; then
    echo \"Return code for index was not zero but  \$retval\"
fi\n" >>$SCRIPTS_DIR/$sample.sb
	echo -en "\necho -en \"Getting Coverage\"" >>$SCRIPTS_DIR/$sample.sb
	echo "bedtools coverage -g $SCRATCH/bowtieDB/MAG_AssemblyChrSizes.txt -sorted -hist -a $SCRATCH/bowtieDB/MAG_Assembly.map.bed -b pseudoalignments_sorted.bam > pseudoalignments.mag.map.hist" >>$SCRIPTS_DIR/$sample.sb
	echo -en "
retval=\$?
if [ \$retval -ne 0 ]; then
    echo \"Return code for bedtools was not zero but  \$retval\"
fi\n" >>$SCRIPTS_DIR/$sample.sb
	echo "python $GLBRC/scripts/get_coverage_for_genes.py -i <(echo pseudoalignments.mag.map.hist) > pseudoalignments.mags.coverage" >>$SCRIPTS_DIR/$sample.sb
	echo "python $GLBRC/scripts/genes.to.kronaTable.py -i $SCRATCH/bowtieDB/mags/PROKKA.MAG_Assembly.ec -m $GLBRC/annotations/metagenomics-workshop/reference_db/kegg/ec.to.pwy -H $GLBRC/annotations/metagenomics-workshop/reference_db/kegg/pwy.hierarchy -n $sample -l <(grep \"minpath 1\" $SCRATCH/bowtieDB/mags/PROKKA.MAG_Assembly.kegg.minpath) -c pseudoalignments.mags.coverage -o $COUNTS/$sample.mags.krona.kegg.minpath.tab" >>$SCRIPTS_DIR/$sample.sb
	
	echo "echo \"Run Completed\"" >>$SCRIPTS_DIR/$sample.sb
	echo -en "$i $sample Job#: "

	short=${sample:0:-4}
	if [ $i = $((nsamples - 1)) ]; then
		sbatch --account shadeash-colej --time=0:30:00 --mem=$MEM --cpus-per-task=$THREADS --job-name=$short --error=$GLBRC/logs/kallisto/$sample.err --output=$GLBRC/logs/kallisto/$sample.out $SCRIPTS_DIR/$sample.sb --mail-user=dooley.shanek@gmail.com --mail-type=ALL 
	else
		sbatch --account shadeash-colej --time=0:30:00 --mem=$MEM --cpus-per-task=$THREADS --job-name=$short --error=$GLBRC/logs/kallisto/$sample.err --output=$GLBRC/logs/kallisto/$sample.out $SCRIPTS_DIR/$sample.sb 
	fi

done

# echo $SCRIPTS_DIR/$sample.sb










# OUTPUT_DIR=$GLBRC/mapping/$READ_TYPE/kallistoOut #Output directory for kallisto counts


############### You shouldn't have to change anything below here ###############
# Step 1. Create the Kallisto index for the metagenomic assembly

# TODAY=`date +%d-%m-%y`
# echo -en "$(date) \t Building Kallisto Index \n\t kallisto index -i $INDEX_NAME $ASSEMBLY\n\n" >kallisto_progress.txt
# cd $INDEX_DIR
# kallisto index -i $INDEX_NAME $ASSEMBLY
# echo -en "$(date) \t Kallisto Index completed\n" >>$OUTPUT_DIR/kallisto_progress.txt

# Step 2. kallisto defaults paired-end and you need to specify pairA_1.fastq pairA_2.fastq pairB_1.fastq pairB_2.fastq, etc... To do that, I am going to do a for loop and create a string to specify all the pairs.


# Step 3. Run Kallisto
# echo -en "$(date) \t Running Kallisto with $nsamples samples ($((nsamples*2)) fastq files) \n\t kallisto quant -t $SLURM_CPUS_PER_TASK -i $INDEX_DIR/$INDEX_NAME -o $SCRATCH_OUT $SAMPLE_STRING\n" #>>$OUTPUT_DIR/kallisto_progressG.txt
# kallisto bus -t $SLURM_CPUS_PER_TASK -i $INDEX_DIR/$INDEX_NAME -o $SCRATCH_OUT $SAMPLE_STRING
# echo -en "$(date) \t Kallisto Run complete\n\nResources:\n" #>>$OUTPUT_DIR/kallisto_progressG.txt
# cp $SCRATCH_OUT/* $OUTPUT_DIR/
#js -j $SLURM_JOB_ID >>$OUTPUT_DIR/kallisto_progress.txt
