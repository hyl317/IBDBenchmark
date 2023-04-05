#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N simGeno #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=50G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -l 'h=!bionode0[1-9]'



python3 simGeno_noIBDcopied.py -n 100