#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N damage #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=45G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:52:1

ch=3
i=$SGE_TASK_ID
id="$(sed "${i}q;d" /mnt/archgen/users/yilei/IBDsim/genotypes/iid.list)"

path2BAM="/mnt/archgen/users/yilei/IBDsim/BAM/$id.bam"

echo damage profiling for $id

ref="/mnt/archgen/users/yilei/Data/1000G_release/raw_data/hs37d5.fa"
/projects1/tools/mapdamage_2.0.6/mapDamage -i $path2BAM -r $ref -Q 30 -v
