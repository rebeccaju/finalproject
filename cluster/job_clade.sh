#!/bin/bash

#SBATCH --partition=eeb354
#SBATCH --job-name=clade_iqtree
#SBATCH --time=30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

module load IQ-TREE/1.6.12

iqtree -s por.co1.alignment.fasta -nt AUTO
