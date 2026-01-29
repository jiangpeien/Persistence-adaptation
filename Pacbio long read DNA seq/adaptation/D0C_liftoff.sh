#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --job-name=D0Cliftoff
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=1G

module load condaenvs/gpu/liftoff
module load minimap2/2.15
cd /pacbio_20230731
liftoff -g /3017gene_NC_007793.gff -o D0C_3017gene.gff -p 12 D0C_polished_assembly.fasta USA300_FPR3757.fasta 
