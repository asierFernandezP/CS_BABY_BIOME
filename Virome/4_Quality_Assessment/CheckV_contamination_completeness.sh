#!/bin/bash
#SBATCH --job-name=CheckV_cont_comp
#SBATCH --output=CheckV_cont_comp.out
#SBATCH --mem=20gb
#SBATCH --time=01:59:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

contig_file=$1 #path to FASTA file with the predicted viral contigs
output=$2 #path to output directory

# Clean environment and load modules
module purge; ml CheckV/1.0.1-foss-2021b-DIAMOND-2.1.8; ml list 

# Run CheckV
checkv end_to_end \
    $contig_file \
    $output \
    -t ${SLURM_CPUS_PER_TASK} \
    -d /scratch/hb-llnext/databases/checkv-db-v1.5
