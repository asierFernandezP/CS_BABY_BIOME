#!/bin/bash
#SBATCH --job-name=geNomad_ME
#SBATCH --output=geNomad_ME.out
#SBATCH --mem=64gb
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=16
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

contig_file=$1 #path to FASTA file with contigs
output_dir=$2 #path to output directory

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /home2/p304845/Conda_envs/geNomad_conda/; conda list

# Run geNomad
genomad end-to-end \
        --conservative \
        --cleanup \
        $contig_file \
        $output_dir \
        /scratch/hb-llnext/databases/geNomad_db/

conda deactivate
