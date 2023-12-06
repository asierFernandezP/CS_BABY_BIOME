#!/bin/bash
#SBATCH --job-name=CheckV_der_contigs 
#SBATCH --output=CheckV_der_contigs.out
#SBATCH --mem=8gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=16
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

contig_file=$1 #path to FASTA file with predicted viral contigs
contig_file_name="$(basename "${contig_file}")" #extract filename
contig_file_name="${contig_file_name%.*}" #extract filename without the extension

echo -e '\n---- RUNNING BLAST+ ----'
mkdir -p Dereplication/
cd Dereplication/

# Clean environment, load modules 
module purge; ml CheckV BLAST+; module list

# Run dereplication for contigs with CheckV scripts

#First, create a blast+ database:
makeblastdb \
	-in $contig_file \
	-dbtype nucl \
	-out ${contig_file_name}_db.fa

#Next, use megablast from blast+ package to perform all-vs-all blastn of sequences:
blastn \
	-query $contig_file \
	-db ${contig_file_name}_db.fa \
	-outfmt '6 std qlen slen' \
	-max_target_seqs 10000 \
	-out ${contig_file_name}_blast.tsv \
	-num_threads ${SLURM_CPUS_PER_TASK}

echo -e '\n---- RUNNING anicalc.py and aniclust.py scripts ----'

#Next, calculate pairwise ANI by combining local alignments between sequence pairs:
/scratch/p304845/CS_Baby_Biome/2_ME_REMOVAL/anicalc.py \
	-i ${contig_file_name}_blast.tsv  \
	-o ${contig_file_name}_ani.tsv

#Finally, perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
/scratch/p304845/CS_Baby_Biome/2_ME_REMOVAL/aniclust.py \
	--fna $contig_file \
	--ani ${contig_file_name}_ani.tsv \
	--out ${contig_file_name}_clusters.tsv \
	--min_ani 95 \
	--min_tcov 85 \
	--min_qcov 0

echo -e '\n---- Dereplication DONE ----'
