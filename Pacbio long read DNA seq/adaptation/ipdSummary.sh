#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --job-name=D1A
#SBATCH --nodes=2
#SBATCH --cpus-per-task=12
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=4G

cd /gpfs/home/pj2062/SCRATCH/pacbio_20230731
#This script is an example to align the D1A reads to D0C and call the modification. We can do the same for all the samples separately
/gpfs/share/apps/smrtlink/5.1.0/smrtcmds/bin/pbmm2 align D0C_polished_assembly.fasta \
/gpfs/data/sequence/results/yanailab/PacBio/r64466e_20230731_182046/Microbial_Genome_Analysis/D1A/ccs_kinetics_bystrandify.subreads.bystrand.bam \
D1A_alignD0C.bam --preset HIFI --sort

/gpfs/share/apps/smrtlink/5.1.0/smrtcmds/bin/pbindex D1A_alignD0C.bam

/gpfs/share/apps/smrtlink/5.1.0/smrtcmds/bin/ipdSummary D1A_alignD0C.bam --methylFraction \
--reference D0C_polished_assembly.fasta --gff D1A_alignD0C_methy.gff
