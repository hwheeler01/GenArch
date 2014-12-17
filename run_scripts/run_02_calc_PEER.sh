#!/bin/bash

#PBS -N R.PEER
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=10
#PBS -l mem=16gb
#PBS -o $HOME/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e $HOME/${PBS_JOBNAME}.e${PBS_JOBID}.err

cd $PBS_O_WORKDIR
module load R

R CMD INSTALL /home/hwheeler/R_peer_source_1.3.tgz
R --vanilla < 02_calc_PEER.r
