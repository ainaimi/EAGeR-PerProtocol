#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=7-00:00:00
#SBATCH --job-name=ipw-gComp-sim
#SBATCH --mem=120g
#SBATCH --partition=week-long-cpu
#SBATCH --array=1-4

module purge
module load R

# Natural Course
Rscript --no-save --no-restore --verbose ./R/gF.R "compl1" 0.7142857 200 10000 0 $SLURM_ARRAY_TASK_ID  > gF_NaturalCourse.Rout 2>&1

# Exposed
Rscript --no-save --no-restore --verbose ./R/gF.R "compl1" 0.7142857 200 10000 0 $SLURM_ARRAY_TASK_ID 1 1 > gF_Exposed.Rout 2>&1

# Unxposed
Rscript --no-save --no-restore --verbose ./R/gF.R "compl1" 0.7142857 200 10000 0 $SLURM_ARRAY_TASK_ID 0 0 > gF_Unxposed.Rout 2>&1

# Interaction 6
Rscript --no-save --no-restore --verbose ./R/gF.R "compl1" 0.7142857 200 10000 0 $SLURM_ARRAY_TASK_ID 1 1 0 6 > gF_Interaction6.Rout 2>&1

# Interaction 8
Rscript --no-save --no-restore --verbose ./R/gF.R "compl1" 0.7142857 200 10000 0 $SLURM_ARRAY_TASK_ID 1 1 0 8 > gF_Interaction8.Rout 2>&1

# Interaction 12
Rscript --no-save --no-restore --verbose ./R/gF.R "compl1" 0.7142857 200 10000 0 $SLURM_ARRAY_TASK_ID 1 1 0 12 > gF_Interaction12.Rout 2>&1

# Interaction 20
Rscript --no-save --no-restore --verbose ./R/gF.R "compl1" 0.7142857 200 10000 0 $SLURM_ARRAY_TASK_ID 1 1 0 20 > gF_Interaction20.Rout 2>&1
