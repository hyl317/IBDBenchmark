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
#$ -t 6:6:1

i=$SGE_TASK_ID
iid="$(sed "${i}q;d" iid.list)"
cd ./$iid/impute_atlas
for ch in {1..22};
do
    #rm ./ch$ch/$iid.chr$ch.vcf.gz
    #rm ./ch$ch/$iid.chr$ch.vcf.gz.csi
    rm ./ch$ch/chunks.chr$ch.txt
    # create (temporary) symbolic link for ease of postprocessing
    bcftools reheader -s ./reheader -o $iid.chr$ch.merged.bcf ./ch$ch/GLIMPSE_ligated/$iid.chr$ch.merged.bcf
    bcftools index $iid.chr$ch.merged.bcf
done

bcftools concat -O b $iid.chr{1..22}.merged.bcf | bcftools view -O z -o $iid.glimpse.GP99.MAF1.vcf.gz -i 'MAX(FORMAT/GP)>=0.99&&INFO/INFO>=0.8&&INFO/RAF>=0.01&&INFO/RAF<=0.99'
bcftools index $iid.glimpse.GP99.MAF1.vcf.gz
# remove symbolic link
rm $iid.chr*.merged.bcf
rm $iid.chr*.merged.bcf.csi