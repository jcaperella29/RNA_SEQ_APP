#!/bin/bash
#SBATCH --job-name=rnaseq_shiny
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12000
#SBATCH --time=02:00:00

SIF_PATH="/mnt/c/Users/jcape/Downloads/RNA_SEQ_APP/rnaseq_shiny.sif"
APP_PATH="/mnt/c/Users/jcape/Downloads/RNA_SEQ_APP"

singularity run --bind $APP_PATH:/srv/shiny-server/app $SIF_PATH
