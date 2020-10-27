#!/bin/sh 
#SBATCH --job-name=2020_06_02_run_cognac_for_genus_Escherichia 
#SBATCH --mail-user=rcrawfo@umich.edu 
#SBATCH --mail-type=BEGIN,END 
#SBATCH --cpus-per-task=12 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --mem-per-cpu=100g 
#SBATCH --time=3-00:00:00 
#SBATCH --account=esnitkin1 
#SBATCH --partition=largemem 
#SBATCH --output=/scratch/esnitkin_root/esnitkin/rcrawfo/cognac_paper/lib/%x-%j.log 
 
cd /scratch/esnitkin_root/esnitkin/rcrawfo/cognac_paper/lib 
Rscript 2020_06_02_run_cognac_for_genus.R Escherichia 
