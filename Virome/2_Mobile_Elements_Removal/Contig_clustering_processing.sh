#!/bin/bash
#SBATCH --job-name=Contig_clustering_plasmid_processing
#SBATCH --output=Clustering_plasmid_proc.out
#SBATCH --mem=4gb
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

clusters=$1 #CheckV aniclust.py all_predicted_viral_contigs_clusters.tsv output file
plasmids=$2 #geNomad all_predicted_viral_contigs_plasmid_summary.tsv output file
viruses=$3 #geNomad all_predicted_viral_contigs_virus_summary.tsv output file
contig_file=$4 #File with all the predicted viral contigs from STEP1

contig_file_path="$(dirname "${contig_file}")" #extract path
contig_file_name="$(basename "${contig_file}")" #extract filename
contig_file_name="${contig_file_name%.*}" #extract filename without the extension

#Load modules
ml Python/3.10.8-GCCcore-12.2.0 seqtk; ml list

#Process geNomad output file to get a file with the names of identified plasmid sequences (1 per row)
awk 'NR>1' $plasmids | cut -f1 | sort > Plasmid_sequences.txt
#Execute the python script that outputs Viral_sequences_no_plasmids.txt file with only the viral contigs
python Contig_clustering_processing.py $clusters Plasmid_sequences.txt

#Extract from the FASTA file with the contigs those that are predicted to be viral (after dereplication)
seqtk subseq $contig_file Viral_sequences_no_plasmids.txt > all_predicted_viral_sequences_no_plasmids.fa #full contigs

#Extract from the FASTA file with the contigs those that are predicted to be plasmids (according to geNomad results)
seqtk subseq $contig_file Plasmid_sequences.txt > all_predicted_plasmids.fa

# Set permissions
chmod 440 all_predicted_plasmids.fa all_predicted_viral_sequences_no_plasmids.fa Plasmid_sequences.txt

exec > Clustering_geNomad_summary_info.txt
echo -e "\n"

echo -e "################################### CLUSTERING/PLASMID REMOVAL SUMMARY STATS ###################################\n"
#1. Get summary stats of geNomad results
n_plasmids_geNomad=$(cat Plasmid_sequences.txt | wc -l)
n_viruses_geNomad=$(cat $viruses | cut -f1 | cut -f1 -d "|" | sort | uniq | wc -l)
n_original_contigs=$(cat $contig_file |grep ">"| wc -l)
echo "The number of contigs classified as plasmids by geNomad is: $n_plasmids_geNomad (Plasmid_sequences.txt)"
echo "The number of contigs classified as viruses by geNomad is: $n_viruses_geNomad"
n_unclassified=$(( $n_original_contigs - $n_plasmids_geNomad - $n_viruses_geNomad ))
echo "The total number of contigs unclassified by geNomad is: $n_unclassified"
n_no_plasmids=$(( $n_unclassified + $n_viruses_geNomad ))
echo -e "The number of non-plasmid contigs according to geNomad (before dereplication) is: $n_no_plasmids\n"

mean_length_plasmids=$(awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' all_predicted_plasmids.fa | grep -v ">" |  awk '{ total += $1; count++ } END { print total/count }')
median_length_plasmids=$(awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' all_predicted_plasmids.fa | grep -v ">" |  sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')
echo -e "The mean length of contigs classified as plasmids by geNomad is: $mean_length_plasmids (median = $median_length_plasmids)\n"

#2. Get summary stats of results after dereplication
n_clusters=$(cat $clusters | wc -l)
echo "The total number of clusters found after dereplication is: $n_clusters"
n_no_plasmids_derep=$(cat Viral_sequences_no_plasmids.txt | wc -l)
echo "The total number of contigs that DO NOT cluster with plasmids identified by geNomad is: $n_no_plasmids_derep (Viral_sequences_no_plasmids.txt)"
n_plasmids_derep=$(( $n_original_contigs -  $n_no_plasmids_derep ))
echo -e "The total number of contigs that cluster with plasmids identified by geNomad is: $n_plasmids_derep\n"

mean_length_no_plasmids=$(awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' all_predicted_viral_sequences_no_plasmids.fa | grep -v ">" |  awk '{ total += $1; count++ } END { print total/count }')
median_length_no_plasmids=$(awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' all_predicted_viral_sequences_no_plasmids.fa | grep -v ">" |  sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')
echo -e "The mean length of predicted viral contigs that DO NOT cluster with plasmids identified by geNomad is: $mean_length_no_plasmids (median = $median_length_no_plasmids)\n"

echo "The file all_predicted_plasmids.fa contains the sequences of all contigs that have been predicted to be plasmids by geNomad."
echo -e "The file all_predicted_viral_sequences_no_plasmids.fa contains the sequences of all contigs that do not cluster with predicted plasmids. It can be used as input for CheckV.\n"
echo -e "###################################################### END ######################################################\n"
