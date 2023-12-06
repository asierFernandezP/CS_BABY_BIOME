#!/bin/bash
#SBATCH --job-name=rRNA_cont_processing
#SBATCH --output=rRNA_cont_processing.out
#SBATCH --mem=1gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

# Clean environment and load modules
module purge; ml seqtk; module list

viral_fragments_file=$1 #path to output FASTA file with viral fragments (STEP3)
viral_file_name="$(basename "${viral_fragments_file}")" #extract filename
viral_file_name="${viral_file_name%.*}" #extract filename without the extension
viral_SILVA_SSU_blast=$2 #Blast results file of viral fragments vs SILVA SSU DB
viral_SILVA_LSU_blast=$3 #Blast results file of viral fragments vs SILVA LSU DB

# Get the contigs IDs with rRNA sequences aligning to them (hit covering > 50% of the gene length)
awk '$8 - $7 + 1 > $14*0.5' $viral_SILVA_SSU_blast | cut -d$'\t' -f2 | sort | uniq > Viral_SSU_contigs.txt
awk '$8 - $7 + 1 > $14*0.5' $viral_SILVA_LSU_blast | cut -d$'\t' -f2 | sort | uniq > Viral_LSU_contigs.txt
cat Viral_LSU_contigs.txt > Viral_rRNA_contigs.txt
comm -13 Viral_LSU_contigs.txt Viral_SSU_contigs.txt >> Viral_rRNA_contigs.txt

# Exclude rRNA contigs from FASTA file with predicted viral regions
grep '>' $viral_fragments_file | sed 's/^>//' | sort > Viral_all_contigs.txt
sort Viral_rRNA_contigs.txt | comm -13 - Viral_all_contigs.txt > ${viral_file_name}_no_rRNA_contigs.txt

seqtk subseq \
    -l 80 \
    $viral_fragments_file \
    ${viral_file_name}_no_rRNA_contigs.txt > ${viral_file_name}_no_rRNA.fa

# Remove intermediate files and set permissions
rm Viral_all_contigs.txt
chmod 440 ${viral_file_name}_no_rRNA_contigs.txt ${viral_file_name}_no_rRNA.fa Viral_rRNA_contigs.txt
