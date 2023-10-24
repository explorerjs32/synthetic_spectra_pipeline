#!/bin/bash
#SBATCH --job-name=line_by_line_abundances              # Job name
#SBATCH --mail-type=END,FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fmendez@ufl.edu                     # Where to send mail
#SBATCH --account=rezzeddine                            # Account name
#SBATCH --qos=rezzeddine-b                              # Request normal or burst cores
#SBATCH --nodes=1                                       # Run all processes on a single node
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=4                               # Number of CPU cores per task
#SBATCH --mem-per-cpu=2gb                               # total memory limit
#SBATCH --time=24:00:00                                 # Time limit hrs:min:sec
#SBATCH --output=../log_files/abund_out_%A_%a.out       # Standard output log
#SBATCH --error ../log_files/abund_error_%A_%a.error    # Standard error log
#SBATCH --array=0-10



date;hostname

# Define the input index for the elements
input_idx=${SLURM_ARRAY_TASK_ID}

scp /blue/rezzeddine/share/fmendez/pipeline/line_by_line_abundances.py /blue/rezzeddine/share/fmendez/pipeline/line_by_line_abundances_${SLURM_ARRAY_TASK_ID}.py

# Define the star and elements you want to get abundances for
star="$1"
elements=("Mg" "Si" "Ca" "Sc" "Ti" "V" "Cr" "Mn" "Co" "Ni" "Fe")

element=${elements[${input_idx}]}

# Calculate the line by line abundances
echo Calculating the ${element} abundance of ${star}

python3 /blue/rezzeddine/share/fmendez/pipeline/line_by_line_abundances_${SLURM_ARRAY_TASK_ID}.py -s ${star} -e ${element}

# Remove the extra files
rm /blue/rezzeddine/share/fmendez/pipeline/line_by_line_abundances_${SLURM_ARRAY_TASK_ID}.py

date