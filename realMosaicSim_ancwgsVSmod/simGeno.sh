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
#$ -t 1:5:1

i=$SGE_TASK_ID
i=$(($i-1))

# a total of 21 lengths
lens=(4 8 12 16 20)
l=${lens[$i]}

echo simulating IBD segments of length $l

python3 simGeno.py -l $l -n 500