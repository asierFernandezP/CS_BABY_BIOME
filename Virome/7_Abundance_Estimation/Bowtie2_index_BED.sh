#!/bin/bash
#SBATCH --job-name=Bowtie2_index
#SBATCH --output=Bowtie2_index.out
#SBATCH --mem=50gb
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

viral_rep_seqs=$1 # FASTA file with all the vOTUs representative viral sequences

# Build bowtie2 index
module purge; ml Bowtie2; module list

bowtie2-build \
    $viral_rep_seqs \
    Viral_rep_seqs_DB \
    --large-index \
    --threads ${SLURM_CPUS_PER_TASK}

# Clean environment and load modules
module purge; ml bioawk/1.0-GCC-11.2.0; module list

# Generate BED file
bioawk -c fastx '{print $name"\t0\t"length($seq)}' $viral_rep_seqs > Viral_rep_seqs_DB.bed

# Change permissions
chmod 440 Viral_rep_seqs_DB.bed
