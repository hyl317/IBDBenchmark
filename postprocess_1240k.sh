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
#$ -t 1:12:1

# post process vcf produced by snpAD
i=$SGE_TASK_ID
iid="$(sed "${i}q;d" iid.list)"
echo postprocessing $iid

DP=20

cd genotypes/$iid
for ch in {1..22}; 
do 
    cd ./chr$ch; 
    vcftools --vcf $iid.chr$ch.vcf --minDP $DP --positions /mnt/archgen/users/hringbauer/data/GRG/vcf/pos/chr$ch.sites --recode --out $iid.chr$ch.minDP$DP.1240k; 
    mv $iid.chr$ch.minDP$DP.1240k.recode.vcf ../;
    cd ../
done

for ch in {1..22};
do
    bgzip $iid.chr$ch.minDP$DP.1240k.recode.vcf
    bcftools index $iid.chr$ch.minDP$DP.1240k.recode.vcf.gz
done

bcftools concat -O v -o $iid.minDP$DP.1240k.all.vcf $iid.chr{1..22}.minDP$DP.1240k.recode.vcf.gz
bgzip $iid.minDP$DP.1240k.all.vcf
bcftools index $iid.minDP$DP.1240k.all.vcf.gz
rm $iid.chr*.vcf.gz
rm $iid.chr*.vcf.gz.csi
