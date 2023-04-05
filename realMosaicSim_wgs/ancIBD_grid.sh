#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N grid #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=50G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:3:1
#$ -l 'h=!bionode0[1-9]'


i=$SGE_TASK_ID
i=$(($i-1))

# a total of 21 lengths
covs=(cov1 cov1over2 cov1over4)
cov=${covs[$i]}

echo grid search for IBD segments of length 12cM with coverage $cov

python3 ancIBD_grid.py -l 12 --cov $cov
python3 ancIBD_grid.py -l 0 --cov $cov