#!/bin/bash
#SBATCH --job-name=CheckV_der_all_viruses
#SBATCH --output=CheckV_der_all_viruses.out
#SBATCH --mem=20gb
#SBATCH --time=08:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

contig_file=$1 #path to FASTA file all the viral sequences (own viral database + public DBs)
blast_db=$2 #path to blast+ database generated from all viral sequences

# Clean environment, load modules
module purge; ml Python/3.10.8-GCCcore-12.2.0 CheckV; module list

# Run dereplication for contigs with CheckV scripts

#Using the blast+ database as input, caalculate pairwise ANI by combining local alignments between sequence pairs:
python /scratch/p304845/CS_Baby_Biome/5_DEREPLICATION_vOTUs/anicalc.py \
    -i $blast_db  \
    -o viruses_ani.tsv

#Finally, perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
python /scratch/p304845/CS_Baby_Biome/5_DEREPLICATION_vOTUs/aniclust.py \
    --fna $contig_file \
    --ani viruses_ani.tsv \
    --out viral_clusters.tsv \
    --min_ani 95 \
    --min_tcov 85 \
    --min_qcov 0
