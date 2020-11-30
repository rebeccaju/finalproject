#!/bin/bash

#SBATCH --partition=eeb354
#SBATCH --job-name=scler_iqtree
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

module load IQ-TREE/1.6.12

iqtree -s alignment.cat.fasta -bb 1000 -nt AUTO
