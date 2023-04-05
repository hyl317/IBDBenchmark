#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N ancIBD_MAF5 #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
# -pe smp 4 #needs 8 CPU cores
#$ -l h_vmem=100G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:50:1

b=$SGE_TASK_ID
echo "running for batch$b"

############################### now merge vcf of the 5 samples into a single vcf ############################

base=/mnt/archgen/users/yilei/IBDsim/downsample/1240k
cd callIBD/1240k
for cov in cov2 cov1 cov3over4 cov1over2;
do
    cd ./$cov
    if [ ! -d batch$b ]; then
        mkdir batch$b
    fi
    cd batch$b
    
    rm -r ./raw
    rm -r ./processed_1KG_MAF5 # delete results from previous trial

    bcftools merge $base/I2105/$cov/batch$b/I2105.batch$b.phased.filter1.vcf.gz \
        $base/I3950/$cov/batch$b/I3950.batch$b.phased.filter1.vcf.gz \
        $base/I5273/$cov/batch$b/I5273.batch$b.phased.filter1.vcf.gz \
        $base/I5279/$cov/batch$b/I5279.batch$b.phased.filter1.vcf.gz | bcftools view -O z -o batch$b.merged.filter1.vcf.gz -e 'FORMAT/GT[*]="./."'
    bcftools index batch$b.merged.filter1.vcf.gz

    ########## now that vcfs are merged, transform it to the appropriate hdf5 format ##########################
    mkdir raw
    mkdir processed_1KG_MAF5
    for ch in {1..22};
    do
        bcftools view -r $ch -O v -o chr$ch.temp.vcf batch$b.merged.filter1.vcf.gz
        python3 /mnt/archgen/users/yilei/LAI/phasing_experiment/AfanasievoFamily/downsample/vcf2hdf5.py --vcf chr$ch.temp.vcf --ch $ch --out processed_1KG_MAF5
        rm chr$ch.temp.vcf
    done
    rm -r ./raw
    cd ../../
done

########################################### now run hapBLOCK using 1000g snps ################################

cd /mnt/archgen/users/yilei/IBDsim/downsample/callIBD/1240k
for cov in cov2 cov1 cov3over4 cov1over2;
do
    cd ./$cov/batch$b
    python3 /mnt/archgen/users/yilei/IBDsim/downsample/hapBLOCK.py --min_cm1 5 --infolder processed_1KG_MAF5 --out ancIBD_1KG_MAF5
    cd ../../
done