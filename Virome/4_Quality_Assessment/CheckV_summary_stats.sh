#!/bin/bash
#SBATCH --job-name=CheckV_summary_stats
#SBATCH --output=CheckV_summary.out
#SBATCH --mem=2gb
#SBATCH --time=00:19:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

CheckV_dir=$1 #directory with CheckV results
cd $CheckV_dir
mkdir CheckV_summary_stats_files
exec > CheckV_summary_stats_files/CheckV_summary_info.txt
echo -e "\n" 

echo -e "################################### CHECKV SUMMARY STATS ###################################\n"
# 1. Get files and summary stats of CheckV results
n_provirus=$(cat proviruses.fna | grep ">" | wc -l)
n_virus=$(cat viruses.fna | grep ">" | wc -l)
echo "The number of viral sequences detected is: $n_virus (viruses.fna file)"
echo "The number of proviral sequences detected is: $n_provirus (proviruses.fna file)"
sum=$(( $n_provirus + $n_virus ))
echo -e "The total number of viral/proviral sequences identified by CheckV is: $sum\n"

# Get the summary stats of the total number of contigs to be kept (completeness > 50%) or discarded. 
n_selected_contigs=$(cat selected_CheckV_contigs.txt| wc -l)
n_filtered_contigs=$(cat filtered_CheckV_contigs.txt| wc -l)
awk 'NR>1' quality_summary.tsv | grep "Not-determined" | cut -f1 | sort > not_determined_comp_CheckV_contigs.txt
n_not_det_contigs=$(cat not_determined_comp_CheckV_contigs.txt| wc -l)
echo "The total number of contigs with completeness >50% is: $n_selected_contigs"
echo -e "The total number of contigs with not-determined completeness is: $n_not_det_contigs\n"
echo "The total number of contigs with completeness =<50% (including Not-determined) is: $n_filtered_contigs"

# 2. Get the file and summary stats of viral contigs selected (completeness > 50%) but with 0 viral genes
awk 'NR>1' quality_summary.tsv | grep "no viral genes detected" | grep "Complete\|High-quality\|Medium-quality" | cut -f1 | sort > selected_nonviralgenes_CheckV_contigs.txt
n_selected_contigs_no_viral_genes=$(cat selected_nonviralgenes_CheckV_contigs.txt| wc -l)
echo -e "The total number of contigs with completeness >50% but 0 viral genes is: $n_selected_contigs_no_viral_genes\n"

# 3. Get files and summary stats of viral contigs with mutiple viral regions

# Get a list of contigs with more than 1 viral sequence detected
# Print the number of contigs
cat quality_summary.tsv  | grep ">1 viral region detected"| cut -f1 | sort > multiple_viral_region_contigs.txt
n_multiple_viral_region_contigs=$(cat multiple_viral_region_contigs.txt| wc -l)
n_multiple_viral_regions=$(sed 's/.*/&_/' multiple_viral_region_contigs.txt | grep -f - proviruses.fna | wc -l) #sed adds an underscore at the end of each line
echo "The number of contigs with more than 1 viral sequence detected is: $n_multiple_viral_region_contigs"
echo -e "The number of viral sequences from the contigs with more than 1 viral sequence detected is: $n_multiple_viral_regions\n"

# Get a list of contigs with more than 1 viral region detected with low quality/undetermined (will be filtered)
# Print the number of contigs and the number of viral regions
awk 'NR>1' quality_summary.tsv | grep ">1 viral region detected" | awk '$8 == "Low-quality" || $8 == "Not-determined"' | cut -f1 | sort > filtered_multiple_viral_region_contigs.txt
n_multiple_viral_region_contigs_filtered=$(cat filtered_multiple_viral_region_contigs.txt| wc -l)
n_multiple_viral_regions_filtered=$(sed 's/.*/&_/' filtered_multiple_viral_region_contigs.txt | grep -f - proviruses.fna | wc -l) #sed adds an underscore at the end of each line
echo "The number of contigs with more than 1 viral sequence detected and with low/undetermined quality is: $n_multiple_viral_region_contigs_filtered"
echo -e "The number of viral sequences from the contigs with more than 1 viral sequence detected and with low/undetermined quality is: $n_multiple_viral_regions_filtered\n"

# Get a list of contigs with more than 1 viral region detected with medium/high-quality (or complete) (will be kept)
# Print the number of contigs and the number of viral regions
awk 'NR>1' quality_summary.tsv | grep ">1 viral region detected" | awk '$8 != "Low-quality" && $8 != "Not-determined"' | cut -f1 | sort > selected_multiple_viral_region_contigs.txt
n_multiple_viral_region_contigs_high_q=$(cat selected_multiple_viral_region_contigs.txt| wc -l)
n_multiple_viral_regions_high_q=$(sed 's/.*/&_/' selected_multiple_viral_region_contigs.txt | grep -f - proviruses.fna| wc -l) #sed adds an underscore at the end of each line
echo "The number of contigs with more than 1 viral sequence detected and with medium/high/complete quality is: $n_multiple_viral_region_contigs_high_q"
echo -e "The number of viral sequences from the contigs with more than 1 viral sequence detected and with medium/high/complete quality is: $n_multiple_viral_regions_high_q\n"

# Double check that all the contigs with multiple viral regions are classified as proviruses
n_multiple_viral_regions_as_viruses=$(grep -wf multiple_viral_region_contigs.txt viruses.fna| wc -l) 
echo -e "The number of viral sequences from the contigs with more than 1 viral sequence detected and NOT classified as proviruses is: $n_multiple_viral_regions_as_viruses\n"

# 4.Get the summary stats of the final number of sequences/contigs selected after running CheckV
# CheckV_sequences.fna will contain full contigs (CheckV_full_contigs.sh) or trimmed viral regions (CheckV_viral_regions.sh).
n_final_sequences=$(cat CheckV_sequences.fna| grep ">" | wc -l)
echo "The final number of viral sequences/contigs with completeness >50% is: $n_final_sequences"
echo -e "The final sequences/contigs are available in CheckV_sequences.fna file\n"

# Set permissions
chmod 440 selected_CheckV_contigs.txt filtered_CheckV_contigs.txt

# Move all the generated summary stats files to a folder
mv *multiple_viral_region_contigs* *CheckV*.txt CheckV_summary_stats_files
echo -e "########################################### END ###########################################\n"
