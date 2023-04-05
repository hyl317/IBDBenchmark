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
# bcftools merge --file-list vcflist | bcftools view -O v -g ^miss | awk -F '\t' '($1~/^#/) || !(($4=="A"&&$5=="G")||($4=="G"&&$5=="A")||($4=="C"&&$5=="T")||($4=="T"&&$5=="C"))' > merged.glimpse.GP99.transversion_only.vcf
# bgzip merged.glimpse.GP99.transversion_only.vcf
# bcftools index merged.glimpse.GP99.transversion_only.vcf.gz


# # 3. convert to plink binary format
plink --vcf merged.glimpse.GP99.transversion_only.vcf.gz --make-bed --out merged.glimpse.GP99.transversion_only
# # add map to .bim file
/mnt/archgen/users/yilei/bin/ibis/add-map-plink.pl merged.glimpse.GP99.transversion_only.bim /mnt/archgen/users/yilei/Data/Hapmap/genetic_map_GRCh37_chr{1..22}.txt > new.bim
mv new.bim merged.glimpse.GP99.transversion_only.bim
# # 4. run IBIS
/mnt/archgen/users/yilei/bin/ibis/ibis -bfile merged.glimpse.GP99.transversion_only -noFamID -o merged.glimpse.GP99.transversion_only -mL 4