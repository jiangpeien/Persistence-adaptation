#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --job-name=pipeline
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=1-20

module load bwa
module load samtools
module load freebayes
module load picard

cd /gpfs/home/pj2062/PROJECTS/gDNA
REF=/gpfs/scratch/pj2062/reference/USA300_FPR3757/USA300_FPR3757.fasta
BAMDIR=bam
VCFDIR=vcf_noduplicates
nodup=nodup
#mkdir -p ${BAMDIR} ${VCFDIR}

# pick the sample prefix for this array task
samp=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_1_20.txt)

echo ">>> Processing $samp (task $SLURM_ARRAY_TASK_ID)"

# 1) Alignment ( for >100bp reads)
bwa mem -t $SLURM_CPUS_PER_TASK $REF \
    trimmed/${samp}_R1_001.trim.fastq.gz \
    trimmed/${samp}_R2_001.trim.fastq.gz \
  | samtools view -bSh - \
  | samtools sort -@ $SLURM_CPUS_PER_TASK -o ${BAMDIR}/${samp}.sorted.bam

samtools index ${BAMDIR}/${samp}.sorted.bam

# 2) Mark (and optionally remove) duplicates
java -jar /gpfs/share/apps/picard/2.26.10/raw/picard/build/libs/picard.jar MarkDuplicates \
  I="bam/${samp}.sorted.bam" \
  O="nodup/${samp}.dedup.bam" \
  M="metrics/${samp}.dup_metrics.txt" \
  REMOVE_DUPLICATES=true \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
  TMP_DIR=tmp 

# 4) Index the deduplicated BAM
samtools index ${nodup}/${samp}.dedup.bam

# 2) Variant calling with freebayes (haploid)

freebayes \
  --fasta-reference $REF \
  --ploidy 1 \
  --min-mapping-quality 20 \
  --min-base-quality 20 \
  --min-alternate-count 5 \
  ${nodup}/${samp}.dedup.bam \
> ${VCFDIR}/${samp}.dedup.freebayes.raw.vcf


echo ">>> Done $samp"
