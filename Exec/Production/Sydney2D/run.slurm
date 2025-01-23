#!/bin/bash
#SBATCH --nodes=2                 # Number of nodes
#SBATCH --time=01:00:00             # Job should run for up to 6 hours
#SBATCH --partition=debug
###SBATCH --gres=gpu:1
###SBATCH --mem=700GB
###SBATCH --qos=high
#SBATCH --account=co2snow  # Where to charge NREL Hours
#SBATCH --output=job.out.%j
#SBATCH --error=job.err.%j

export OMP_NUM_THREADS=4

# Compute number of MPI ranks
let "NUM_RANKS = ($SLURM_JOB_NUM_NODES * $SLURM_CPUS_ON_NODE) / $OMP_NUM_THREADS"
# DEBUG
echo $NUM_RANKS

# Run the code
srun -n $NUM_RANKS ./PeleLMeX3d.gnu.x86-spr.MPI.OMP.ex  inputs.3d-regt $RESTART_TEXT amr.max_wall_time=0.9

