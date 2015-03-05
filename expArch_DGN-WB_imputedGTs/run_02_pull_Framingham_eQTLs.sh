#!/bin/bash

#PBS -N parse.framing
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err

cd $PBS_O_WORKDIR
#module load python
module load perl

##time python 02_pull_Framingham_eQTLs_faster.py ##WAY TOO SLOW, perl >150X faster
time perl 02_pull_Framingham_eQTLs.pl
