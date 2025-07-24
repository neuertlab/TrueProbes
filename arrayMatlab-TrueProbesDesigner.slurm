mod#!/bin/bash
#SBATCH --mail-user=hughes.j.hughes@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --nice=10000
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=300G
#SBATCH --array=49
#SBATCH --output=Out_%A_%a.out

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load blast+/2.14.1
module load matlab/2023b
matlab -nodisplay -nosplash -r  "A0_BKJH_ProbeDesign_Wrapper_cluster_V5(${SLURM_ARRAY_TASK_ID},1);"
