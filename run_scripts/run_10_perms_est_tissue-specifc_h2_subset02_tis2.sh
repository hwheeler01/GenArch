#!/bin/bash
#PBS -N R.gcta.ts
#PBS -S /bin/bash
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=4gb
#PBS -o $HOME/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e $HOME/${PBS_JOBNAME}.e${PBS_JOBID}.err

cd $PBS_O_WORKDIR
module load R
module load gcta/1.24.4

time R --no-save < 10_perms_est_tissue-specifc_h2.r --args 02 2
