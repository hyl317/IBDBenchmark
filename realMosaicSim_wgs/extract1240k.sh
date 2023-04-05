#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N postprocess #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=50G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID

## Extract glimpse imputed genotypes on chr3 for the 1240k sites to establish "groundtruth genotypes" for these samples

path2sites="/mnt/archgen/users/hringbauer/data/GRG/vcf/pos/chr3.sites"

for i in {1..52};
do
    iid="$(sed "${i}q;d" iid.list)"
    echo postprocessing $iid
    path2vcf="/mnt/archgen/users/yilei/IBDsim/genotypes/$iid/impute_trim5/ch3/$iid.chr3.merged.vcf.gz"
    vcftools --gzvcf $path2vcf --positions $path2sites  --recode --out ./genotypes/$iid.chr3.1240k; 
    mv ./genotypes/$iid.chr3.1240k.recode.vcf ./genotypes/$iid.chr3.1240k.vcf
    bgzip ./genotypes/$iid.chr3.1240k.vcf
    bcftools index ./genotypes/$iid.chr3.1240k.vcf.gz
done