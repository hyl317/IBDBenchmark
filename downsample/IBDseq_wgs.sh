#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N IBDseq #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
# -pe smp 12 #needs 8 CPU cores
#$ -l h_vmem=25G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:50:1
#$ -l 'h=!bionode0[1-9]'


b=$SGE_TASK_ID
echo "running IBD seq for batch$b"


# for id in I2105 I3950 I5273 I5279;
# do
#     echo doing $id
#     for cov in cov5 cov2 cov1 cov3over4 cov1over2 cov1over4 cov1over10;
#     do
#         cd ./wgs/$id/$cov/batch$b
#         # for ch in {1..22};
#         # do
#         #     bcftools annotate -a ./GLIMPSE_ligated/$id.batch$b.chr$ch.merged.bcf -c FORMAT/GP ./GLIMPSE_phased/$id.batch$b.chr$ch.phased.bcf -Oz -o $id.chr$ch.merged.vcf.gz
#         #     bcftools index $id.chr$ch.merged.vcf.gz
#         # done
#         # echo $id > reheader
#         # bcftools concat -Ob $id.chr{1..22}.merged.vcf.gz | bcftools reheader -s reheader | bcftools view -Oz -o $id.batch$b.phased.vcf.gz
#         # bcftools index -f $id.batch$b.phased.vcf.gz
#         # # remove temporary file
#         # rm reheader
#         # rm $id.chr*.merged.vcf.gz
#         # rm $id.chr*.merged.vcf.gz.csi
#         # bcftools view -O z -o $id.batch$b.phased.filter1.vcf.gz -i 'MAX(FORMAT/GP)>=0.99&&INFO/INFO>=0.8&&INFO/RAF>=0.05&&INFO/RAF<=0.95' $id.batch$b.phased.vcf.gz
#         # bcftools index -f $id.batch$b.phased.filter1.vcf.gz
#         bcftools view -O z -o $id.batch$b.phased.ibdseq.vcf.gz -i 'INFO/RAF>=0.01&&INFO/RAF<=0.99&&INFO/INFO>=0.8' $id.batch$b.phased.vcf.gz
#         bcftools index -f $id.batch$b.phased.ibdseq.vcf.gz
#         cd ../../../../
#     done
# done




# base=/mnt/archgen/users/yilei/IBDsim/downsample/wgs
# cd callIBD/wgs

# for cov in cov5 cov2 cov1 cov3over4 cov1over2 cov1over4 cov1over10;
# do
#     cd ./$cov
#     if [ ! -d batch$b ]; then
#         mkdir batch$b
#     fi
#     cd batch$b

#     if [ ! -d IBDseq ]; then
#         mkdir IBDseq
#     fi
#     cd IBDseq

#     bcftools merge $base/I2105/$cov/batch$b/I2105.batch$b.phased.ibdseq.vcf.gz \
#     $base/I3950/$cov/batch$b/I3950.batch$b.phased.ibdseq.vcf.gz \
#     $base/I5273/$cov/batch$b/I5273.batch$b.phased.ibdseq.vcf.gz \
#     $base/I5279/$cov/batch$b/I5279.batch$b.phased.ibdseq.vcf.gz \
#     /mnt/archgen/users/yilei/Data/1000G_release/raw_data/1kg.all.autosome.vcf.gz \
#     | bcftools view -O u -S /mnt/archgen/users/yilei/IBDsim/downsample/callIBD/IBDseq.iids -e 'FORMAT/GT[*]="./."' \
#     | bcftools view -O z -o batch$b.merged.ibdseq.vcf.gz -i 'MAC>1'

#     bcftools index batch$b.merged.ibdseq.vcf.gz

#     # for ch in {1..22};
#     # do
#     #     java -Xmx75g -jar /mnt/archgen/users/yilei/bin/ibdseq.r1206.jar gt=batch$b.merged.ibdseq.vcf.gz out=batch$b.IBDseq.ch$ch chrom=$ch nthreads=5
#     # done
#     # rm batch$b.merged.ibdseq.vcf.gz # to save space
#     # rm batch$b.merged.ibdseq.vcf.gz.csi

#     cd ../../../
# done

# cd callIBD/wgs
# for cov in cov5 cov2 cov1 cov3over4 cov1over2 cov1over4 cov1over10;
# do
#     cd ./$cov/batch$b/IBDseq
#     for ch in {1..22};
#     do
#         if [ -s batch$b.IBDseq.ch$ch.log ]
#         then
#             echo log file of ch$ch already exists, skipping...
#         else
#             java -Xmx110g -jar /mnt/archgen/users/yilei/bin/ibdseq.r1206.jar gt=batch$b.merged.ibdseq.vcf.gz out=batch$b.IBDseq.ch$ch chrom=$ch nthreads=12
#         fi
#     done
#     cd ../../../
# done

####### clean up
cd callIBD/wgs
for cov in cov5 cov2 cov1 cov3over4 cov1over2 cov1over4 cov1over10;
do
    cd ./$cov/batch$b/IBDseq
    rm *.hbd
    rm *.log
    rm *.r2.filtered
    for ch in {1..22};
    do
        grep -v HG batch$b.IBDseq.ch$ch.ibd | grep -v NA > ch$ch.tmp.ibd
    done
    cat ch*.tmp.ibd > batch$b.IBDseq.ibd
    rm *.tmp.ibd
    rm batch$b.IBDseq.ch*.ibd
    cd ../../../
don