#!/bin/bash
## job name 
#SBATCH --job-name=interpolate

#SBATCH --mail-user=bperfect@uw.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


#SBATCH --account=stf
#SBATCH --partition=stf
## Resources 
## Nodes 
#SBATCH --nodes=1       
## Tasks per node (Slurm assumes you want to run 28 tasks per node unless explicitly told otherwise)
#SBATCH --ntasks-per-node=28 
## Walltime (hh:mm:ss) 
#SBATCH --time=3:00:00 
## Memory per node 
#SBATCH --mem=120G 
## Specify the working directory for this job 
#SBATCH --workdir=/gscratch/stf/bperfect/

## Modules needed to run
module load matlab_2017a

matlab -nodisplay -nosplash -nodesktop -r "run /gscratch/stf/bperfect/hyak_interpolate_simname.m"; exit; > log_interp_simname
