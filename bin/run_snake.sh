#!/bin/bash
#SBATCH --export=NONE
#SBATCH -J verify_crispr_ont
#SBATCH -o logs/%x.o%j
#SBATCH -e logs/%x.e%j
#SBATCH --ntasks 1
#SBATCH --time 120:00:00
#SBATCH --mem=8G
#SBATCH --partition=bbc

cd $SLURM_SUBMIT_DIR

snakemake_module="bbc2/snakemake/snakemake-8.25.2"

module load $snakemake_module

# make logs dir if it does not exist already. 
logs_dir="logs/"
[[ -d $logs_dir ]] || mkdir -p $logs_dir


echo "Start snakemake workflow." >&1                   
echo "Start snakemake workflow." >&2     

snakemake \
-p \
--latency-wait 120 \
--use-envmodules \
--use-conda \
--jobs 100 \
--executor cluster-generic --cluster-generic-submit-cmd "mkdir -p logs/{rule}; sbatch \
-p bbc,long \
--export=ALL \
--nodes 1 \
--ntasks-per-node {threads} \
--mem={resources.mem_gb}G \
-t 120:00:00 \
-o logs/{rule}/{resources.log_prefix}.o \
-e logs/{rule}/{resources.log_prefix}.e" # SLURM hangs if output dir does not exist, so we create it before running sbatch on the snakemake jobs.
#--slurm \
#--default-resources slurm_account=${SLURM_JOB_USER} slurm_partition=${SLURM_JOB_PARTITION}

echo "snakemake workflow done." >&1                   
echo "snakemake workflow done." >&2                
