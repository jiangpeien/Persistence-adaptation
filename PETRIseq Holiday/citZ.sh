#!/bin/bash
#SBATCH --job-name=citZ_select
#SBATCH --nodes=2
#SBATCH --cpus-per-task=12
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=8G

export PATH=/gpfs/home/pj2062/miniconda/envs/PETRI/bin:/gpfs/home/pj2062/miniconda2/envs/PETRI/bin:/gpfs/home/pj2062/miniconda/bin:/gpfs/home/pj2062/miniconda2/bin:/cm/shared/apps/slurm/current/sbin:/cm/shared/apps/slurm/current/bin:/cm/local/apps/environment-modules/4.4.1/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/sbin:/cm/local/apps/environment-modules/current/bin:/opt/ibutils/bin:/usr/lpp/mmfs/bin:/gpfs/home/pj2062/.local/bin:/gpfs/home/pj2062/bin

module load samtools
mkdir PETRI_holiday_cluster_variant/citZ_selected 
cd PETRI_holiday_cluster_variant/citZ_selected 

threads=12       # Adjust this based on your CPU resources

# Directory containing all subdirectories with fastq files
BASE_PATH="PETRI_JE2_D2D4D23_20250117"

# File containing barcodes
BARCODE_FILE="PETRI_holiday_cluster_variant/citZ_D2D4D23_barcodes.csv"

# Specify the chromosome and position
chromosome="NC_007793.1"
position=1799669

# Output file
output_file="PETRI_holiday_cluster_variant/citZ_D2D4D23_1799669base.csv"

# Path to reference genome
REF_GENOME="PETRI_JE2_D2D4D23_20250117/scripts/USA300_FPR3757.fa"

# Initialize the output file with headers
echo -e "Sample\tG\tA\tT\tC" > $output_file

# Function to count bases from mpileup output
count_bases() {
    local mpileup_output=$1
    local sample_name=$2

    # Initialize counts
    local g_count=0
    local a_count=0
    local t_count=0
    local c_count=0

    if [[ -n "$mpileup_output" ]]; then
        # Iterate over each base in the pileup column
        for base in $(echo "$mpileup_output" | awk '{print $5}' | fold -w1); do
            case "$base" in
                G|g) ((g_count++)) ;;
                A|a) ((a_count++)) ;;
                T|t) ((t_count++)) ;;
                C|c) ((c_count++)) ;;
            esac
        done
    fi

    # Write the counts to the output file
    echo -e "${sample_name}\t${g_count}\t${a_count}\t${t_count}\t${c_count}" >> $output_file
}

# Read each barcode from the barcode file
while IFS= read -r barcode
do
    echo "Processing barcode: $barcode"
    # Construct the filename from the barcode
    # Insert "_R2" before the first underscore and after the prefix
    prefix=$(echo "$barcode" | cut -d'_' -f1)  # This isolates the part before the first underscore (BS22735A)
    suffix=$(echo "$barcode" | cut -d'_' -f2-) # This captures everything after the first underscore (bc1_11_bc2_83_bc3_73)
    fastq_file="${BASE_PATH}/${prefix}_R2_trimmed/${prefix}_R2_${suffix}_R2_trimmed.fastq.gz"

    # Alignment to reference genome
    ALIGN_OUTPUT="${barcode}.bam" #sorted indexed bam
    bwa aln -n 1 -t "$threads" "$REF_GENOME" "$fastq_file" > temp.sai 
    bwa samse -n 14 "$REF_GENOME" temp.sai "$fastq_file" > temp.sam 
    samtools view  -@ "$threads" -Sb temp.sam | samtools sort -@ "$threads" -o "$ALIGN_OUTPUT"
    samtools index "$ALIGN_OUTPUT"
    
    # Deduplication using umi_tools
    DEDUP_OUTPUT="${barcode}.dedup.bam"
    umi_tools dedup -I "$ALIGN_OUTPUT" -S "$DEDUP_OUTPUT"
    samtools index "$DEDUP_OUTPUT"
    
    rm temp.sai temp.sam "$ALIGN_OUTPUT" "${barcode}.bam.bai"
    
    # Get the mpileup output for the specific position
    mpileup_output=$(samtools mpileup -r ${chromosome}:${position}-${position} "$DEDUP_OUTPUT" | grep "^${chromosome}")

    # Count bases and write to the output file
    count_bases "$mpileup_output" "$barcode"

    echo "Finished processing barcode: $barcode"
done < "$BARCODE_FILE"


echo "All barcodes processed. Mpileup completed for all samples."

