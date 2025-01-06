#!/bin/bash

#<------------------------Request for Resources----------------------->
#SBATCH -p parallel
#SBATCH --mem=20G
#SBATCH -t 160:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#<------------------------Set environment variables------------------->

module load R/4.4.1

#<------------------------Run R script with arguments----------------->
Rscript semi-p12-n100-c50-r05.R $1 $2
