#!/bin/bash
#PBS -N R.gcta
#PBS -S /bin/bash
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=4gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err

cd $PBS_O_WORKDIR
module load R
module load gcta/1.24.4

time R --no-save < 09_perms_est_cross-tissue_h2.r --args 15
