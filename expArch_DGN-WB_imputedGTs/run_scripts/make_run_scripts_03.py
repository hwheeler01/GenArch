#!/usr/bin/env python

'''make a run script for each chr and output a qsub file'''

qsubfile = open('../qsub.txt','w')
prescript = '03_calc_GRMs'

for i in range(1,23):
    outfilename = 'run_' + prescript + '_chr' + str(i) + '.sh'
    outfile = open(outfilename,'w')
    output = '''#!/bin/bash
#PBS -N R.gcta.''' + str(i) +'''\n#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load R
module load gcta/1.24.4

time R --vanilla < ''' + prescript + '.r --args ' + str(i)
    outfile.write(output)
    qsubfile.write('qsub run_scripts/' + outfilename + '\nsleep 3\n')


