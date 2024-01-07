#!/bin/bash
# set the number of nodes.
#SBATCH --nodes=1
# set the number of CPUs required.
#SBATCH --ntasks-per-node=48
# set the amount of memory needed for each CPU.
#SBATCH --mem-per-cpu=8000
# set max wallclock time (hh:mm:ss).
#SBATCH --time=00:10:00
# set the time partition for the job. 
#SBATCH --partition=devel
# set name of job (AND DATE!)
#SBATCH --job-name=debug_param_sweep
# mail alert at start, end, and abortion of execution
#SBATCH --mail-type=ALL
# send mail to this address
#SBATCH --mail-user=james.hammond@merton.ox.ac.uk
# run the application

module load Julia/1.8.2-linux-x86_64

julia requirements.jl
julia debug_param_sweep.jl > debug_param_sweep.log
