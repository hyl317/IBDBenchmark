#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N gridsearch #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -pe smp 4 #needs 8 CPU cores
#$ -l h_vmem=30G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 101:150:1
#$ -tc 125


i=$SGE_TASK_ID
i=$(($i-1))

coverages=(cov5 cov2 cov1 cov3over4 cov1over2)
cov_index=$(($i/50))
cov=${coverages[$cov_index]}
b=$((1+$i-50*$cov_index))
echo "doing grid search for batch$b, $cov"


cd ./callIBD/$cov/batch$b

#######################################
# clean up previous runs
# cd ./grid
# ls -1 | xargs rm -rf
# cd ../
# rm -r grid 
#######################################

python3 /mnt/archgen/users/yilei/IBDsim/downsample/hapBLOCK_grid.py --min_cm1 5
rm -r ./grid