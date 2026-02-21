#!/bin/bash
#PBS -N wave_scaling
#PBS -q cpu
#PBS -l select=4:ncpus=8:mpiprocs=16
#PBS -l walltime=02:00:00
#PBS -e error_scaling.log
#PBS -o output_scaling.log

# Navigate to the working directory
cd /work/u10905938/Wave_Equation

# Execute the python script directly through the container
apptainer exec --bind /work /software/containers/mk/amsc_mk_2025.sif bash -c "
    source /u/sw/etc/bash.bashrc
    module load gcc-glibc dealii
    python3 scripts/weak_scaling.py
"
