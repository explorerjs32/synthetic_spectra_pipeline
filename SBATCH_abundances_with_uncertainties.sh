#!/bin/bash
#SBATCH --job-name=abundances                           # Job name
#SBATCH --mail-type=END,FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fmendez@ufl.edu                     # Where to send mail
#SBATCH --qos=rezzeddine-b                           # Request normal or burst cores
#SBATCH --nodes=1                                       # Run all processes on a single node
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=2                               # Number of CPU cores per task
#SBATCH --mem-per-cpu=2gb                                      # total memory limit
#SBATCH --time=24:00:00                                 # Time limit hrs:min:sec
#SBATCH --output=../log_files/abund_out_%A_%a.out        # Standard output log
#SBATCH --error ../log_files/abund_error_%A_%a.error     # Standard error log
#SBATCH --array=1-21                                    # Array range

date;hostname

# Define the input index for the elements
input_idx=$1

# Create the copies for all the necesary files
scp /blue/rezzeddine/share/fmendez/pipeline/abundances_with_uncertainties.py /blue/rezzeddine/share/fmendez/pipeline/abundances_with_uncertainties_${input_idx}_${SLURM_ARRAY_TASK_ID}.py
scp /blue/rezzeddine/share/fmendez/pipeline/abund_spec_create.sh /blue/rezzeddine/share/fmendez/pipeline/abund_spec_create_${input_idx}_${SLURM_ARRAY_TASK_ID}.sh
scp /blue/rezzeddine/share/fmendez/Turbospectrum2019-master/COM-v19.1/ts_script_abu.com /blue/rezzeddine/share/fmendez/Turbospectrum2019-master/COM-v19.1/ts_script_abu_${input_idx}_${SLURM_ARRAY_TASK_ID}.com
scp /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth.sh /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth_${input_idx}_${SLURM_ARRAY_TASK_ID}.sh
scp /blue/rezzeddine/share/fmendez/pipeline/inst_broadening.f90 /blue/rezzeddine/share/fmendez/pipeline/inst_broadening_${input_idx}_${SLURM_ARRAY_TASK_ID}.f90
scp /blue/rezzeddine/share/fmendez/pipeline/rot_vel.f90 /blue/rezzeddine/share/fmendez/pipeline/rot_vel_${input_idx}_${SLURM_ARRAY_TASK_ID}.f90

# Define the star and elements you want to get abundances for
star="$2"
elements=("Mg" "Si" "Ca" "Sc" "Ti" "V" "Cr" "Mn" "Co" "Ni")

element=${elements[${input_idx}]}

# Get the solar abundance from Axplund 2020
sol_abundance=`grep -H -m1 $element /blue/rezzeddine/share/fmendez/literature_tables/asplund_solar_abund.csv | cut -d' ' -f3`

# Define the abundances range
step=1
min=$(awk '{print $1-$2}' <<<"${sol_abundance} ${step}")
max=$(awk '{print $1+$2}' <<<"${sol_abundance} ${step}")

# Create an array of abundances
abund_arr=( $(seq -f "%f" ${min} 0.1 ${max}) )

# Assign each of the abundance to the code and create the synthetic spectra
idx_dif=1
index=$(awk '{print $1-$2}' <<<"${SLURM_ARRAY_TASK_ID} ${idx_dif}")

echo ${star} ${element} ${sol_abundace} ${abund_arr[$index]}

# Define the variables to input to the abundance code
ts_bash=abund_spec_create_${input_idx}_${SLURM_ARRAY_TASK_ID}.sh
ts_script=ts_script_abu_${input_idx}_${SLURM_ARRAY_TASK_ID}.com
abundance=${abund_arr[$index]}
broadening_bash=spectra_smooth_${input_idx}_${SLURM_ARRAY_TASK_ID}.sh
inst_broad_script=inst_broadening_${input_idx}_${SLURM_ARRAY_TASK_ID}.f90
inst_broad_script_out=inst_broadening_${input_idx}_${SLURM_ARRAY_TASK_ID}
rot_vel_broad_script=rot_vel_${input_idx}_${SLURM_ARRAY_TASK_ID}.f90
rot_vel_broad_script_out=rot_vel_${input_idx}_${SLURM_ARRAY_TASK_ID}

# Derive the abundances
echo Running abundances_with_uncertainties_${input_idx}_${SLURM_ARRAY_TASK_ID}.py for ${star} ${element} using a ${abund_arr[$index]} abundance value
python3 /blue/rezzeddine/share/fmendez/pipeline/abundances_with_uncertainties_${input_idx}_${SLURM_ARRAY_TASK_ID}.py model_create.sh scratch ${ts_bash} ${ts_script} ${abundance} ${broadening_bash} ${inst_broad_script} ${inst_broad_script_out} ${rot_vel_broad_script} ${rot_vel_broad_script_out} ${star} ${element} &

bg_pid=$!
wait $bg_pid
wait


# Delete the the extra files
rm /blue/rezzeddine/share/fmendez/pipeline/abundances_with_uncertainties_${input_idx}_${SLURM_ARRAY_TASK_ID}.py
rm /blue/rezzeddine/share/fmendez/pipeline/abund_spec_create_${input_idx}_${SLURM_ARRAY_TASK_ID}.sh
rm /blue/rezzeddine/share/fmendez/Turbospectrum2019-master/COM-v19.1/ts_script_abu_${input_idx}_${SLURM_ARRAY_TASK_ID}.com
rm /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth_${input_idx}_${SLURM_ARRAY_TASK_ID}.sh
rm /blue/rezzeddine/share/fmendez/pipeline/inst_broadening_${input_idx}_${SLURM_ARRAY_TASK_ID}.f90
rm /blue/rezzeddine/share/fmendez/pipeline/rot_vel_${input_idx}_${SLURM_ARRAY_TASK_ID}.f90
rm /blue/rezzeddine/share/fmendez/pipeline/inst_broadening_${input_idx}_${SLURM_ARRAY_TASK_ID}
rm /blue/rezzeddine/share/fmendez/pipeline/rot_vel_${input_idx}_${SLURM_ARRAY_TASK_ID}


date
