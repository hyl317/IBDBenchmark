#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N copy #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=10G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID

#git add --all -v
#git commit -v -m "all relevant files, excluding large files"
#git push -v --set-upstream origin master

find . -iregex '.*\.\(sh\|py\|ipynb\)$' -exec cp --parents {} /mnt/archgen/users/yilei/script_IBDbenchmark \;