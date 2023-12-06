#!/bin/bash
#SBATCH --job-name=vOTU_clustering_processing
#SBATCH --output=vOTU_clustering_proc.out
#SBATCH --mem=4gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

clusters=$1 #CheckV aniclust.py STEP5_combined_sequences_clusters.tsv output file
viral_genomes_file=$2 #file with all the viral genomes used as input for the dereplication

#Load modules
ml Python/3.10.8-GCCcore-12.2.0 seqtk; ml list

# Generate a txt and a FASTA file with the representative sequences of each vOTU 
cat $clusters | cut -f1 > rep_seqs_vOTUs.txt
seqtk subseq $viral_genomes_file rep_seqs_vOTUs.txt > rep_seqs_vOTUs.fa

# Set permissions 
chmod 440 rep_seqs_vOTUs.txt rep_seqs_vOTUs.fa
