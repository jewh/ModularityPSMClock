#!/bin/bash
# set the number of nodes.
#SBATCH --nodes=44
# set the number of CPUs required.
#SBATCH --ntasks-per-node=48
# set the amount of memory needed for each CPU.
#SBATCH --mem-per-cpu=8000
# set max wallclock time (hh:mm:ss).
#SBATCH --time=30:00:00
# set the time partition for the job. 
#SBATCH --partition=medium
# set name of job (AND DATE!)
#SBATCH --job-name=20230829_changing_freq_and_TM
# mail alert at start, end, and abortion of execution
#SBATCH --mail-type=ALL
# send mail to this address
#SBATCH --mail-user=james.hammond@merton.ox.ac.uk
# run the application

module load Julia/1.8.2-linux-x86_64

julia requirements.jl
julia map_parameter_space.jl > map_parameter_space.log
