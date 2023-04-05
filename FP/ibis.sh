#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N ibis #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=45G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID

# bcftools merge -Ou -l vcf.file.list | bcftools reheader -s reheader.txt | bcftools view -Oz -o AGDP.ch3.GP99.MAF5.vcf.gz -g ^miss  
# bcftools index AGDP.ch3.GP99.MAF5.vcf.gz
# plink --vcf AGDP.ch3.GP99.MAF5.vcf.gz --make-bed --out AGDP.ch3
# /mnt/archgen/users/yilei/bin/ibis/add-map-plink.pl AGDP.ch3.bim /mnt/archgen/users/yilei/Data/Hapmap/genetic_map_GRCh37_chr3.txt > new.bim
# rm AGDP.ch3.bim
# mv new.bim AGDP.ch3.bim
# run IBIS
/mnt/archgen/users/yilei/bin/ibis/ibis -bfile AGDP.ch3 -noFamID -o AGDP.ch3 -min_l 4