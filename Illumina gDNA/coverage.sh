#!/bin/bash
#SBATCH --partition=cpu_medium
#SBATCH --job-name=cov
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=14:00:00

set -euo pipefail
cd /gpfs/home/pj2062/PROJECTS/gDNA
REF="/gpfs/scratch/pj2062/reference/USA300_FPR3757/USA300_FPR3757.fasta"           # <-- set this
DIR="nodup_rg"
OUTDIR="picard_wgs_metrics"
#mkdir -p "$OUTDIR"

module load picard
# One-time reference indexes (safe to rerun)
#java -jar /gpfs/share/apps/picard/2.26.10/raw/picard/build/libs/picard.jar CreateSequenceDictionary R="$REF" O="${REF%.*}.dict"

# 1) Run CollectWgsMetrics on each BAM
for i in {4..20}; do
  BAM="$DIR/$i.bam"
  [ -f "$BAM" ] || { echo "Missing $BAM" >&2; continue; }
  java -Xmx6g -jar /gpfs/share/apps/picard/2.26.10/raw/picard/build/libs/picard.jar CollectWgsMetrics \
    I="$BAM" O="$OUTDIR/$i.wgs.txt" R="$REF" \
    MINIMUM_BASE_QUALITY=20 MINIMUM_MAPPING_QUALITY=20 COVERAGE_CAP=1000000
done

# 2) Extract MEDIAN_COVERAGE into one table
echo -e "sample\tmedian_coverage" > wgs_median_coverage.tsv
for MET in "$OUTDIR"/*.wgs.txt; do
  s=$(basename "$MET" .wgs.txt)
  med=$(awk -F'\t' '
    /^#/ {next}
    !h && NF { for(i=1;i<=NF;i++) if($i=="MEDIAN_COVERAGE"){c=i; h=1; next} }
    h && NF { print $c; exit }
  ' "$MET")
  echo -e "${s}\t${med:-NA}"
done >> wgs_median_coverage.tsv

