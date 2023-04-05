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
#$ -t 1:6:1

i=$SGE_TASK_ID
i=$(($i-1))

# a total of 27 lengths
lens=(2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 11 12 13 14 15 16 17 18 19 20)
l=${lens[$i]}

echo simulating IBD segments of length $l

python3 simGeno.py -l $l -n 250