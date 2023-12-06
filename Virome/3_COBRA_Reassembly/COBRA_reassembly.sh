#!/bin/bash
#SBATCH --job-name=COBRA
#SBATCH --output=COBRA_%A_%a.out
#SBATCH --mem=8gb
#SBATCH --time=01:29:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

FASTQ_dir=$1 #directory with the FASTQ files with clean reads
sample_list=${FASTQ_dir}/$2 #file with the list of all samples in the directory
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})
contig_dir=$3 # directory with the metaSPAdes contigs/scaffolds 
viral_contig_file=$4 # FASTA file with the predicted viral contigs

echo -e '\n-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'

echo -e '\n---- Preparing COBRA input files for '${SAMPLE_ID}' SAMPLE ----'

# Clean environment, load environment 
module purge; module load SeqKit; module list

mkdir -p 3_COBRA/
mkdir -p 3_COBRA/SAM_files 3_COBRA/BAM_files 3_COBRA/COVERAGE_files 3_COBRA/RESULT_files 3_COBRA/TMP_files 3_COBRA/LOG_files 3_COBRA/OUTPUT_files

# Generate a FASTA file with the viral sequences identified from the sample
seqkit grep -r -p $SAMPLE_ID $viral_contig_file | sed 's/^[^_]*_/>/' > 3_COBRA/TMP_files/${SAMPLE_ID}_viruses.fa

echo -e '\n---- Preparing mapping file for '${SAMPLE_ID}' SAMPLE ----'

# Clean environment, load environment 
module purge; module load Bowtie2 SAMtools; module list


bowtie2-build ${contig_dir}/${SAMPLE_ID}_metaspades_scaffolds.fa ${SAMPLE_ID}

bowtie2 \
	--very-sensitive \
	-x ${SAMPLE_ID} \
	-1 ${FASTQ_dir}/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
	-2 ${FASTQ_dir}/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
	--no-unal \
	--threads ${SLURM_CPUS_PER_TASK} \
	-S 3_COBRA/SAM_files/${SAMPLE_ID}.sam 

samtools view -bS 3_COBRA/SAM_files/${SAMPLE_ID}.sam  > 3_COBRA/BAM_files/${SAMPLE_ID}.bam
samtools sort -o 3_COBRA/BAM_files/${SAMPLE_ID}_sorted.bam 3_COBRA/BAM_files/${SAMPLE_ID}.bam
samtools index -@ $((${SLURM_CPUS_PER_TASK}-1)) 3_COBRA/BAM_files/${SAMPLE_ID}_sorted.bam

echo -e '\n---- Preparing coverage file for '${SAMPLE_ID}' SAMPLE ----'

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/CoverM_env; conda list

coverm contig \
    --coupled ${FASTQ_dir}/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz ${FASTQ_dir}/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
    --reference ${contig_dir}/${SAMPLE_ID}_metaspades_scaffolds.fa \
    -o 3_COBRA/COVERAGE_files/${SAMPLE_ID}_coverage.txt

sed -i '1d' 3_COBRA/COVERAGE_files/${SAMPLE_ID}_coverage.txt # remove header from CoverM output coverage file

echo -e '\n---- Running COBRA for '${SAMPLE_ID}' SAMPLE ----'

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3 BLAST+; module list
source activate /scratch/hb-llnext/conda_envs/COBRA_env; conda list

python /scratch/hb-llnext/conda_envs/COBRA_env/cobra/cobra.py \
    -f ${contig_dir}/${SAMPLE_ID}_metaspades_scaffolds.fa \
    -q 3_COBRA/TMP_files/${SAMPLE_ID}_viruses.fa \
    -c 3_COBRA/COVERAGE_files/${SAMPLE_ID}_coverage.txt \
    -m 3_COBRA/BAM_files/${SAMPLE_ID}_sorted.bam \
    -a metaspades \
    -mink 21 \
    -maxk 55 \
    -o 3_COBRA/RESULT_files/${SAMPLE_ID}

mv 3_COBRA/RESULT_files/${SAMPLE_ID}/log 3_COBRA/LOG_files/${SAMPLE_ID}_log.txt 
sed "s/^>/>${SAMPLE_ID}_/" 3_COBRA/RESULT_files/${SAMPLE_ID}/*.fasta >> 3_COBRA/RESULT_files/all_predicted_viral_sequences_extended.fa

# Remove intermediate files and directories
rm ${SAMPLE_ID}*bt2 3_COBRA/SAM_files/${SAMPLE_ID}.sam 3_COBRA/BAM_files/${SAMPLE_ID}*bam 3_COBRA/BAM_files/${SAMPLE_ID}*bai
rm -r 3_COBRA/RESULT_files/${SAMPLE_ID}/intermediate.files 3_COBRA/RESULT_files/${SAMPLE_ID}/*fasta 3_COBRA/RESULT_files/${SAMPLE_ID}/debug.txt  3_COBRA/RESULT_files/${SAMPLE_ID}/*new.fa 3_COBRA/TMP_files/${SAMPLE_ID}_viruses.fa
find 3_COBRA/SAM_files 3_COBRA/BAM_files 3_COBRA/TMP_files -type d -empty -delete
