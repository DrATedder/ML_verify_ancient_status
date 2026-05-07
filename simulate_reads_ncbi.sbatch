#!/bin/sh
#SBATCH --job-name=InSilicoSeq
#SBATCH --account=43070_203920
#SBATCH --partition=compute
#SBATCH --cpus-per-task=20
#SBATCH --time=18:00:00

module load bioconda
conda activate inSilicoSeq

list=(5 6 7 8 9 10 11 12 13 15 15 16 17 18 19 20)

for i in {1..100};
   do
   echo UoBsim${i}
   iss generate --cpus ${SLURM_CPUS_PER_TASK} --ncbi bacteria --n_genomes_ncbi ${list[$RANDOM % ${#list[@]} ]} --model hiseq --output /storage02/or-microbio/simulated_modern_reads/new_sims/UoBsim${i}
done
