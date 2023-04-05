#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N ibis #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=50G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID

##################################### prepare input for IBIS ##################################

# 1. Merge vcf from all samples into one
bcftools merge --file-list vcflist | bcftools view -O z -o merged.wgs.vcf.gz -g ^miss -i 'MAC>=1'


# # 3. convert to plink binary format
# plink --vcf bronze_age.all.nomiss.vcf --make-bed --out bronze_age
# # add map to .bim file
# /mnt/archgen/users/yilei/bin/ibis/add-map-plink.pl bronze_age.bim /mnt/archgen/users/yilei/Data/Hapmap/genetic_map_GRCh37_chr{1..22}.txt > new.bim
# # 4. run IBIS
# /mnt/archgen/users/yilei/bin/ibis/ibis -bfile bronze_age -noFamID -o ibis.BA