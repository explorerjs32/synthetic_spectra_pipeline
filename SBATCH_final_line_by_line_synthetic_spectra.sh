#!/bin/bash
#SBATCH --job-name=synthspec_abundances                           # Job name
#SBATCH --mail-type=END,FAIL                                      # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fmendez@ufl.edu                               # Where to send mail
#SBATCH --account=rezzeddine                                      # Account name
#SBATCH --qos=rezzeddine-b                                        # Request normal or burst cores
#SBATCH --nodes=1                                                 # Run all processes on a single node
#SBATCH --ntasks=1                                                # Run a single task
#SBATCH --cpus-per-task=2                                         # Number of CPU cores per task
#SBATCH --mem-per-cpu=2gb                                         # total memory limit
#SBATCH --time=24:00:00                                           # Time limit hrs:min:sec
#SBATCH --output=../log_files/abund_synthspec_out_%A_%a.out       # Standard output log
#SBATCH --error ../log_files/abund_synthspec_error_%A_%a.error    # Standard error log
#SBATCH --array=0-10                                              # Array range



date;hostname

# Define the input index for the elements
input_idx=${SLURM_ARRAY_TASK_ID}

# Create the copies for all the necesary files
scp /blue/rezzeddine/share/fmendez/pipeline/final_line_by_line_abundance_synthetic_spectra.py /blue/rezzeddine/share/fmendez/pipeline/final_line_by_line_abundance_synthetic_spectra_${input_idx}.py
scp /blue/rezzeddine/share/fmendez/pipeline/abund_spec_create.sh /blue/rezzeddine/share/fmendez/pipeline/abund_spec_create_${input_idx}.sh
scp /blue/rezzeddine/share/fmendez/Turbospectrum2019-master/COM-v19.1/ts_script_abu.com /blue/rezzeddine/share/fmendez/Turbospectrum2019-master/COM-v19.1/ts_script_abu_${input_idx}.com
scp /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth.sh /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth_${input_idx}.sh
scp /blue/rezzeddine/share/fmendez/pipeline/inst_broadening.f90 /blue/rezzeddine/share/fmendez/pipeline/inst_broadening_${input_idx}.f90
scp /blue/rezzeddine/share/fmendez/pipeline/rot_vel.f90 /blue/rezzeddine/share/fmendez/pipeline/rot_vel_${input_idx}.f90

# Activate the copies of the bash scripts
chmod a+x /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth_${input_idx}.sh

# Define the star and elements you want to get abundances for
star="$1"
elements=("Fe" "Mg" "Si" "Ca" "Sc" "Ti" "V" "Cr" "Mn" "Co" "Ni")

element=${elements[${input_idx}]}

# Define the variables to input to the abundance code
ts_bash=abund_spec_create_${input_idx}.sh
ts_script=ts_script_abu_${input_idx}.com
broadening_bash=spectra_smooth_${input_idx}.sh
inst_broad_script=inst_broadening_${input_idx}.f90
inst_broad_script_out=inst_broadening_${input_idx}
rot_vel_broad_script=rot_vel_${input_idx}.f90
rot_vel_broad_script_out=rot_vel_${input_idx}

# Create the abundance synthetic spectra
python3 /blue/rezzeddine/share/fmendez/pipeline/final_line_by_line_abundance_synthetic_spectra_${input_idx}.py -tsbs ${ts_bash} -tss ${ts_script} -bbs ${broadening_bash} -ibs ${inst_broad_script_out} -rbs ${rot_vel_broad_script_out} -s ${star} -e ${element} &

bg_pid=$!
wait $bg_pid
wait


# Delete the the extra files
rm /blue/rezzeddine/share/fmendez/pipeline/final_line_by_line_abundance_synthetic_spectra_${input_idx}.py
rm /blue/rezzeddine/share/fmendez/pipeline/abund_spec_create_${input_idx}.sh
rm /blue/rezzeddine/share/fmendez/Turbospectrum2019-master/COM-v19.1/ts_script_abu_${input_idx}.com
rm /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth_${input_idx}.sh
rm /blue/rezzeddine/share/fmendez/pipeline/inst_broadening_${input_idx}.f90
rm /blue/rezzeddine/share/fmendez/pipeline/rot_vel_${input_idx}.f90
rm /blue/rezzeddine/share/fmendez/pipeline/inst_broadening_${input_idx}
rm /blue/rezzeddine/share/fmendez/pipeline/rot_vel_${input_idx}


date
