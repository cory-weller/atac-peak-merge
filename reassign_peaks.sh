#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 50G
#SBATCH --time 24:00:00
#SBATCH --partition norm

module purge

module load R/4.3
module load samtools

job_start=${1}
job_end=${2}

if [ -z "${job_end}" ]; then
    job_end=${job_start}
fi

job_list=$(seq $job_start $job_end)
max_jobs=4

parallel -j ${max_jobs} Rscript reassign_peaks.R ::: ${job_list}

