#!/bin/bash
#SBATCH --job-name=rRNA_DB
#SBATCH --output=rRNA_DB.out
#SBATCH --mem=8gb
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

viral_fragments_file=$1 #path to output FASTA file with viral fragments (STEP3)
viral_file_name="$(basename "${viral_fragments_file}")" #extract filename
viral_file_name="${viral_file_name%.*}" #extract filename without the extension

# Clean environment, load BLAST+
module purge; ml BLAST+; module list

# First, create a blast+ database of the viral fragments:
makeblastdb \
    -dbtype 'nucl' \
    -in $viral_fragments_file \
    -out ${viral_file_name}_db

# Next, use megablast from blast+ package to perform all-vs-all blastn of sequences:
# Using the SILVA_138.1_SSURef_NR99 DB
blastn \
    -task 'blastn' \
    -evalue 0.001 \
    -query /scratch/p304845/CS_Baby_Biome/4_rRNA_FILTERING/SILVA_138.1_SSURef_NR99_tax_silva.fasta \
    -db ${viral_file_name}_db \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen slen' \
    -out ${viral_file_name}_SILVA_SSU_blast.tsv \
    -num_threads ${SLURM_CPUS_PER_TASK}

# Using the SILVA_138.1_LSURef_NR99 DB
blastn \
    -task 'blastn' \
    -evalue 0.001 \
    -query /scratch/p304845/CS_Baby_Biome/4_rRNA_FILTERING/SILVA_138.1_LSURef_NR99_tax_silva.fasta \
    -db ${viral_file_name}_db \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen slen' \
    -out ${viral_file_name}_SILVA_LSU_blast.tsv \
    -num_threads ${SLURM_CPUS_PER_TASK}

# Set file permissions
chmod 440 *_blast.tsv
