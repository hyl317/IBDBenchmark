#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N germline #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=40G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:1:1


b=$SGE_TASK_ID
echo "running hapIBD for batch$b"

##### For germline, use the same input as hapIBD (aka 'filter 1')


base=/mnt/archgen/users/yilei/IBDsim/downsample/wgs
cd callIBD/wgs

for cov in cov5 cov2 cov1 cov3over4 cov1over2 cov1over4 cov1over10;
do
    cd ./$cov
    if [ ! -d batch$b ]; then
        mkdir batch$b
    fi
    cd batch$b

    rm -r ./GERMLINE
    if [ ! -d GERMLINE ]; then
        mkdir GERMLINE
    fi
    cd GERMLINE

    ## Convert vcf to IMPUTE/SHAPEIT format
    for ch in {1..22};
    do
        bcftools view -Oz -o tmp.batch$b.ch$ch.vcf.gz ../hapIBD/batch$b.merged.hapIBD.vcf.gz $ch
        bcftools index tmp.batch$b.ch$ch.vcf.gz
        bcftools convert --hapsample tmp.ch$ch.batch$b tmp.batch$b.ch$ch.vcf.gz

        ## Running hapIBD
        /mnt/archgen/users/yilei/bin/germline2/g2 tmp.ch$ch.batch$b.hap.gz tmp.ch$ch.batch$b.samples /mnt/archgen/users/yilei/Data/Hapmap/GERMLINE/germline_map_chr$ch.txt batch$b.ch$ch
    done
    rm tmp.*
    cd ../../../
done