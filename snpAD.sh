#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N snpAD #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -pe smp 20 #needs 8 CPU cores
#$ -l h_vmem=100G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:22:1
# -tc 11

id="I3949"
echo doing$id
mkdir ./genotypes/$id
cd ./genotypes/$id
ch=$SGE_TASK_ID
mkdir chr$ch
cd chr$ch

ref="/mnt/archgen/users/yilei/Data/1000G_release/raw_data/hs37d5.fa"
#path2bam="../../../BAM/$id.bam"
path2bam="/mnt/ancient/AGDP/release_v0/bams/all/I3949.bam"
Bam2snpAD -r $ch -f $ref -Q 30 -q 30 $path2bam > $id.chr$ch.snpAD
/mnt/archgen/users/yilei/bin/intersectbed.pl $id.chr$ch.snpAD ../../../map_track/chr$ch.bed > $id.chr$ch-2.snpAD
mv $id.chr$ch-2.snpAD $id.chr$ch.snpAD
snpAD -c 20 -o priors.txt -O errors.txt $id.chr$ch.snpAD
snpADCall -N $id -e errors.txt -p "`cat priors.txt`" $id.chr$ch.snpAD > $id.chr$ch.vcf

