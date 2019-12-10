filename=$1


for f in mapped/*.sam; 
do 
	bamname=${f/mapped/bams}.bam
	bamname=${bamname/.sam/}
	sortedbam=${bamname/.bam/.sorted.bam}
	statsFile=${sortedbam/bams/stats}
	statsFile=${statsFile/.sorted.bam/.tsv}

	if [ ! -f $statsFile ]; then
		# echo "samtools view  --threads 12 -f 0x2 -q 20 -o $bamname $f"
		# echo "samtools sort  --threads 12 -o $sortedbam $bamname"
		# echo "rm $bamname"
		# echo "samtools index -@ 12 $sortedbam"
		echo "samtools idxstats $sortedbam > $statsFile"
	fi
done