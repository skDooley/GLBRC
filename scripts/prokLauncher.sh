#!/bin/bash --login
#SBATCH --time=24:00:00            # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=6                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=16         # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=50G                  # memory required per node - amount of memory (in bytes)
#SBATCH --job-name bammaker        # you can give your job a name for easier identification (same as -J)
#SBATCH --account shadeash-colej 

cd /mnt/research/germs/shane/glbrc  
batchn=0
for fname in scripts/hpc_scripts/annotation/*.sb;
do
    ((batchn++))
    echo $batchn $fname
    #endFileName=`expr index "$fname" \.`-9
    #jname="Blast_${fname:8:endFileName}"
    sbatch --job-name=Prokka$batchn -A shadeash-colej --cpus-per-task=16 --mem=50G --output=logs/prokka_3_19_$batchn.slurm.log --ntasks=1 --time=48:00:00 $fname
  	
done





