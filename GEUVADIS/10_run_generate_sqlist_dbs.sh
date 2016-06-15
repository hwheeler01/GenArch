#!/bin/bash
#PBS -N make_db
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load python/2.7.9

./generate_sqlite_dbs.py --input_folder /group/im-lab/nas40t2/hwheeler/cross-tissue/ --results_sub_folder GEUVADIS --betas_sub_folder GEUVADIS/geu_weights/ --output_folder geu_db_files/ --betas_include_clause elasticNet_alpha1 --results_include_clause elasticNet_alpha1

