#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N ibis #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
# -pe smp 2 #needs 8 CPU cores
#$ -l h_vmem=50G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:50:1

b=$SGE_TASK_ID
echo "running IBIS for batch$b"


# do some variants filtering

# for id in I2105 I3950 I5273 I5279;
# do
#     for cov in cov5 cov2 cov1 cov3over4 cov1over2;
#     do
#         cd ./$id/$cov/batch$b
#         rm $id.batch$b.phased.filter*.vcf.gz
#         rm $id.batch$b.phased.filter*.vcf.gz.csi
#         bcftools view -O z -o $id.batch$b.phased.filter1.vcf.gz -i 'MAX(FORMAT/GP)>=0.99&&INFO/INFO>=0.8&&INFO/RAF>=0.05&&INFO/RAF<=0.95' $id.batch$b.phased.vcf.gz
#         bcftools index $id.batch$b.phased.filter1.vcf.gz
#         cd ../../../
#     done
# done


base=/mnt/archgen/users/yilei/IBDsim/downsample/1240k
cd callIBD/1240k

for cov in cov2 cov1 cov3over4 cov1over2;
do
    cd ./$cov
    if [ ! -d batch$b ]; then
        mkdir batch$b
    fi
    cd batch$b

    if [ ! -d IBIS_MAF5 ]; then
        mkdir IBIS_MAF5
    fi
    cd IBIS_MAF5

    echo "I2105.batch$b.bam    I2105" > reheader
    echo "I3950.batch$b.bam    I3950" >> reheader
    echo "I5273.batch$b.bam    I5273" >> reheader
    echo "I5279.batch$b.bam    I5279" >> reheader

    outfile="batch$b.merged.GP99.MAF5.vcf.gz"
    bcftools merge $base/I2105/$cov/batch$b/I2105.batch$b.phased.filter1.vcf.gz \
    $base/I3950/$cov/batch$b/I3950.batch$b.phased.filter1.vcf.gz \
    $base/I5273/$cov/batch$b/I5273.batch$b.phased.filter1.vcf.gz \
    $base/I5279/$cov/batch$b/I5279.batch$b.phased.filter1.vcf.gz | bcftools reheader -s reheader \
    | bcftools view -O z -s I2105,I5273,I5279,I3950 -e 'FORMAT/GT[*]="./."' -o $outfile
    bcftools index $outfile
    rm reheader

    # bcftools merge $base/I2105/$cov/batch$b/I2105.batch$b.phased.filter1.vcf.gz \
    # $base/I3950/$cov/batch$b/I3950.batch$b.phased.filter1.vcf.gz \
    # $base/I5273/$cov/batch$b/I5273.batch$b.phased.filter1.vcf.gz \
    # $base/I5279/$cov/batch$b/I5279.batch$b.phased.filter1.vcf.gz | bcftools reheader -s reheader \
    # | bcftools view -O v -s I2105,I5273,I5279,I3950 -e 'FORMAT/GT[*]="./."' \
    # | awk -F '\t' '($1~/^#/) || !(($4=="A"&&$5=="G")||($4=="G"&&$5=="A")||($4=="C"&&$5=="T")||($4=="T"&&$5=="C"))' > batch$b.merged.GP99.transversion_only.vcf
    
    # bgzip batch$b.merged.GP99.transversion_only.vcf
    # bcftools index batch$b.merged.GP99.transversion_only.vcf.gz
    # rm reheader

    plink --vcf $outfile --make-bed --out batch$b
    /mnt/archgen/users/yilei/bin/ibis/add-map-plink.pl batch$b.bim /mnt/archgen/users/yilei/Data/Hapmap/genetic_map_GRCh37_chr{1..22}.txt > new.bim
    rm batch$b.bim
    mv new.bim batch$b.bim
    # run IBIS
    /mnt/archgen/users/yilei/bin/ibis/ibis -bfile batch$b -noFamID -o batch$b -min_l 5

    # clean up, to save some storage space
    rm $outfile
    rm $outfile.csi

    cd ../../../
done


