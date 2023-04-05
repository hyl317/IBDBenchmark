#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N hapIBD #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -pe smp 5 #needs 8 CPU cores
#$ -l h_vmem=80G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:50:1


b=$SGE_TASK_ID
echo "running hapIBD for batch$b"

# for id in I2105 I3950 I5273 I5279;
# do
#     for cov in cov5 cov2 cov1 cov3over4 cov1over2 cov1over4 cov1over10;
#     do
#         cd ./wgs/$id/$cov/batch$b
#         bcftools view -O z -o $id.batch$b.phased.filter2.vcf.gz -i 'MAX(FORMAT/GP)>=0.99&&INFO/INFO>=0.8&&INFO/RAF>=0.02&&INFO/RAF<=0.98' $id.batch$b.phased.vcf.gz
#         bcftools index $id.batch$b.phased.filter2.vcf.gz
#         cd ../../../../
#     done
# done


base=/mnt/archgen/users/yilei/IBDsim/downsample/wgs
cd callIBD/wgs

for cov in cov5 cov2 cov1 cov3over4 cov1over2 cov1over4 cov1over10;
do
    cd ./$cov
    if [ ! -d batch$b ]; then
        mkdir batch$b
    fi
    cd batch$b

    # rm -r ./hapIBD
    if [ ! -d hapIBD ]; then
        mkdir hapIBD
    fi
    cd hapIBD

    # bcftools merge $base/I2105/$cov/batch$b/I2105.batch$b.phased.filter1.vcf.gz \
    # $base/I3950/$cov/batch$b/I3950.batch$b.phased.filter1.vcf.gz \
    # $base/I5273/$cov/batch$b/I5273.batch$b.phased.filter1.vcf.gz \
    # $base/I5279/$cov/batch$b/I5279.batch$b.phased.filter1.vcf.gz \
    # | bcftools view -Oz -o batch$b.merged.hapIBD.vcf.gz -e 'FORMAT/GT[*]="./."'

    # bcftools index batch$b.merged.hapIBD.vcf.gz

    ## Running hapIBD
    java -Xmx40g -jar /mnt/archgen/users/yilei/bin/hap-ibd.jar min-output=2 gt=batch$b.merged.hapIBD.vcf.gz map=/mnt/archgen/users/yilei/Data/Hapmap/plink.all.GRCh37.map out=batch$b min-seed=0.1 min-extend=0.05 min-mac=1 max-gap=500000

    cd ../../../
done