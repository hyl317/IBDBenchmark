#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N germline #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=120G #request 4Gb of memory
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

    # rm -r ./GERMLINE
    if [ ! -d GERMLINE ]; then
        mkdir GERMLINE
    fi
    cd GERMLINE

    # bcftools view -m2 -M2 -v snps -Oz -o batch$b.subset.vcf.gz -s I2105,I3950,I5273,I5279 ../IBDseq/batch$b.merged.ibdseq.vcf.gz

    # vcftools --gzvcf batch$b.subset.vcf.gz --plink --out batch$b

    # plink --vcf batch$b.subset.vcf.gz --make-just-bim --out batch$b
    # /mnt/archgen/users/yilei/bin/ibis/add-map-plink.pl batch$b.bim /mnt/archgen/users/yilei/Data/Hapmap/genetic_map_GRCh37_chr{1..22}.txt > new.bim
    # rm batch$b.bim
    # mv new.bim batch$b.bim

    # ## update the map file
    # awk 'BEGIN{OFS="\t"}{print $1,$1":"$4,$3,$4}' batch$b.bim > batch$b.map

    /mnt/archgen/users/yilei/bin/germline-1-5-3/germline -input batch$b.ped batch$b.map -output batch$b.germline -g_extend -err_het 10 -err_hom 10

    cd ../../../
done