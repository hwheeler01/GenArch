#!/usr/bin/env python

'''make a run script for each chr and output a qsub file'''

qsubfile = open('../qsub.txt','w')
prescript = '06_calc_knownTrans_GRMs'

for i in range(0,20):
    newi = str(i).zfill(2) ##zfill adds 0 in front of single digits
    outfilename = 'run_' + prescript + '_' + newi + '.sh'
    outfile = open(outfilename,'w')
    output = '''#!/bin/bash
#PBS -N R.gcta.''' + newi +'''\n#PBS -S /bin/bash
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=10gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load R
module load gcta/1.24.4

time R --no-save < ''' + prescript + '.r --args ' + newi + '\n'
    outfile.write(output)
    qsubfile.write('qsub run_scripts/' + outfilename + '\nsleep 3\n')


