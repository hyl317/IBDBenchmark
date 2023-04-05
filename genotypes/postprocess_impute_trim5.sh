#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N glimipse_post #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=45G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID

for iid in I2105 I3950 I5273 I5279;
do
    cd ./$iid/impute_trim5
    for ch in {1..22};
    do
        rm ./ch$ch/$iid.chr$ch.vcf.gz
        rm ./ch$ch/$iid.chr$ch.vcf.gz.csi
        rm ./ch$ch/chunks.chr$ch.txt
        bcftools annotate -a ./ch$ch/GLIMPSE_ligated/$iid.chr$ch.merged.bcf -c FORMAT/GP ./ch$ch/GLIMPSE_phased/$iid.chr$ch.phased.bcf -Oz -o $iid.chr$ch.merged.vcf.gz
        bcftools index $iid.chr$ch.merged.vcf.gz
    done
    bcftools concat -O b $iid.chr{1..22}.merged.vcf.gz | bcftools view -O z -o $iid.glimpse.GP99.MAF2.vcf.gz -i 'MAX(FORMAT/GP)>=0.99&&INFO/INFO>=0.8&&INFO/RAF>=0.02&&INFO/RAF<=0.98'
    bcftools index $iid.glimpse.GP99.MAF2.vcf.gz
    # remove temporary file
    rm $iid.chr*.merged.vcf.gz
    rm $iid.chr*.merged.vcf.gz.csi
    cd ../../
done