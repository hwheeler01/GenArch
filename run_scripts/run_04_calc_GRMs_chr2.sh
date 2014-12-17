#!/bin/bash

#PBS -N R.gcta
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb
#PBS -o $HOME/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e $HOME/${PBS_JOBNAME}.e${PBS_JOBID}.err

cd $PBS_O_WORKDIR

module load R
module load gcta/1.24.4

time R --vanilla < 04_calc_GRMs.r --args 2
