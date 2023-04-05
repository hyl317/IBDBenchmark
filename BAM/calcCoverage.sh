#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N samtools_coverage #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=20G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:52:1

i=$SGE_TASK_ID
path2bam="$(sed "${i}q;d" bamfile.list)"
echo calculating coverage for $path2bam
samtools coverage -q 30 -Q 30 -o $path2bam.coverage $path2bam