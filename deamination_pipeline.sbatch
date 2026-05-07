!/bin/sh
#SBATCH --job-name=deamination_pipeline
#SBATCH --account=43070_203920
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00

module purge
module load bioconda
conda activate deamination

INPUT=/storage02/or-microbio/meta_sim_test_25
OUTPUT=/storage02/or-microbio/meta_sim_test_25

cd ${INPUT}

for r1 in "${INPUT}"/UoB*_R1.fastq
do
    filename=$(basename "$r1")
    name="${filename%%_*}"
    echo "$name"
    simulate_damage --r1 ${INPUT}/${name}_R1.fastq --r2 ${INPUT}/${name}_R2.fastq --out_r1 ${OUTPUT}/a${name}_R1.fastq --out_r2 ${OUTPUT}/a${name}_R2.fastq --ct5 0.03 --ctint 0.4 --ga3 0.01 --gaint 0.7 --overhang 10
done
