#!/bin/bash
## job name 
#SBATCH --job-name=simname_datarun

#SBATCH --mail-user=bperfect@uw.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


#SBATCH --account=stf
#SBATCH --partition=stf
## Resources 
## Nodes 
#SBATCH --nodes=2       
## Tasks per node (Slurm assumes you want to run 28 tasks per node unless explicitly told otherwise)
#SBATCH --ntasks-per-node=28 
## Walltime (hh:mm:ss) 
#SBATCH --time=12:00:00 
## Memory per node 
#SBATCH --mem=120G 
## Specify the working directory for this job 
#SBATCH --workdir=/gscratch/stf/bperfect/

## Modules needed to run
module load icc_17-impi_2017
module load netcdf_fortran+c_4.4.1.1-icc_17 
##old=f2e5n2e4
##new=simname
##find /gscratch/stf/bperfect -exec touch {} \;
##cp seamount.in seamount_${new}.in
##sed -i -e 's/'$old'/'$new'/g' seamount_${new}.in

mpirun -np 56 /gscratch/stf/bperfect/oceanM /gscratch/stf/bperfect/seamount_simname.in > log_simname
