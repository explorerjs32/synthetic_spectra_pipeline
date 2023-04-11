#!/bin/bash
#SBATCH --job-name=MCMC                                # Job name
#SBATCH --mail-type=END,FAIL                           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fmendez@ufl.edu                    # Where to send mail
#SBATCH --nodes=1                                      # Run all processes on a single node
#SBATCH --ntasks=1                                     # Run a single task
#SBATCH --cpus-per-task=4                              # Number of CPU cores per task
#SBATCH --mem=10gb                                     # total memory limit
#SBATCH --time=240:00:00                               # Time limit hrs:min:sec
#SBATCH --output=../log_files/mcmc_out_%A_%a.out       # Standard output log
#SBATCH --error ../log_files/mcmc_error_%A_%a.error    # Standard error log

date;hostname;pwd

# Run the MCMC code
python3 /blue/rezzeddine/share/fmendez/pipeline/mcmc.py

date
