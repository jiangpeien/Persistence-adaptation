#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --job-name=121liftoff
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=1G

module load condaenvs/gpu/liftoff
module load minimap2/2.15
cd /pamA_pacbio
liftoff -g /3017gene_NC_007793.gff -o RU121_polished_assembly_RC_3017gene.gff -p 12 RU121_polished_assembly_RC.fasta USA300_FPR3757.fasta 
