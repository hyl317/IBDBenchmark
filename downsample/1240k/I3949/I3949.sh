#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N I3949 #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
# -pe smp 2 #needs 8 CPU cores
#$ -l h_vmem=45G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 201:300:1
#$ -tc 100
#$ -l 'h=!bionode0[1-9]'


b=$SGE_TASK_ID

j=$(($SGE_TASK_ID-1))
i=$(($j/50))
b=$(($j-50*$i))
b=$((1+$b))
id="I3949"
echo doing$id
echo batch$b

# coverage of I2105 is ~20.4x
covs=(cov2 cov1 cov3over4 cov1over2 cov1over4 cov1over10)
targets=(2.0 1.0 0.75 0.5 0.25 0.1)
path2bam=/mnt/archgen/users/yilei/IBDsim/BAM/1240k/$id.bam
path2bed=/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/target/1240k.all.target.bed

cov=${covs[$i]}
targetCov=${targets[$i]}
frac=`samtools depth -a -b $path2bed -q 30 -Q 30 $path2bam | awk '{print $3}' | ~/stat.pl | awk -v cov=$targetCov '($1=="avg:"){print cov/$2}'`
echo doing coverage $cov with sampling fraction $frac



if [ ! -d $cov ]; then
    mkdir $cov
fi
cd ./$cov

if [ ! -d batch$b ]; then
    mkdir batch$b
fi
cd batch$b
pwd

# actual work start here
samtools view --subsample-seed $b --subsample $frac -b -o $id.batch$b.bam $path2bam
samtools index $id.batch$b.bam
for ch in {1..22}; 
do
    echo processing chr$ch
    #prepare GL file (with bcftools) for GLIMPSE
    BAM=$id.batch$b.bam
    VCF=/mnt/archgen/users/yilei/Data/1000G_release/reference_panel/1000GP.chr$ch.sites.vcf.gz
    TSV=/mnt/archgen/users/yilei/Data/1000G_release/reference_panel/1000GP.chr$ch.sites.tsv.gz
    REFGEN=/mnt/archgen/users/yilei/Data/1000G_release/raw_data/hs37d5.fa
    OUT=$id.batch$b.chr$ch.vcf.gz

    bcftools mpileup -f ${REFGEN} --ignore-RG -I -E -a 'FORMAT/DP' -T ${VCF} -r $ch -q 30 -Q 30 ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
    bcftools index -f ${OUT}

    # running GLIMPSE
    path2glimpse="/mnt/archgen/users/yilei/bin/GLIMPSE"
    VCF=$id.batch$b.chr$ch.vcf.gz
    $path2glimpse/chunk/bin/GLIMPSE_chunk --input $VCF --region $ch --window-size 2000000 --buffer-size 200000 --output chunks.chr$ch.txt
    REF=/mnt/archgen/users/yilei/Data/1000G_release/reference_panel/1000GP.chr$ch.bcf
    MAP=$path2glimpse/maps/genetic_maps.b37/chr$ch.b37.gmap.gz

    if [ ! -d GLIMPSE_imputed ]; then
        mkdir GLIMPSE_imputed
    fi

    while IFS="" read -r LINE || [ -n "$LINE" ];
    do
        printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)
        OUT=GLIMPSE_imputed/$id.batch$b.chr$ch.${ID}.bcf
        $path2glimpse/phase/bin/GLIMPSE_phase --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
        bcftools index -f ${OUT}
    done < chunks.chr$ch.txt

    if [ ! -d GLIMPSE_ligated ]; then
        mkdir GLIMPSE_ligated
    fi
    LST=GLIMPSE_ligated/list.chr$ch.txt
    ls GLIMPSE_imputed/$id.batch$b.chr$ch.*.bcf > ${LST}
    OUT=GLIMPSE_ligated/$id.batch$b.chr$ch.merged.bcf
    $path2glimpse/ligate/bin/GLIMPSE_ligate --input ${LST} --output $OUT
    bcftools index -f ${OUT}
    rm -r ./GLIMPSE_imputed

    if [ ! -d GLIMPSE_phased ]; then
        mkdir GLIMPSE_phased
    fi
    VCF=GLIMPSE_ligated/$id.batch$b.chr$ch.merged.bcf
    OUT=GLIMPSE_phased/$id.batch$b.chr$ch.phased.bcf
    $path2glimpse/sample/bin/GLIMPSE_sample --input ${VCF} --solve --output ${OUT}
    bcftools index -f ${OUT}

    # merge the phasing and imputation file together
    bcftools annotate -a ./GLIMPSE_ligated/$id.batch$b.chr$ch.merged.bcf -c FORMAT/GP ./GLIMPSE_phased/$id.batch$b.chr$ch.phased.bcf -Oz -o $id.chr$ch.merged.vcf.gz
    bcftools index $id.chr$ch.merged.vcf.gz
done
rm chunks.chr*.txt
rm $id.batch$b.chr*.vcf.gz
rm $id.batch$b.chr*.vcf.gz.csi
bcftools concat -Oz -o $id.batch$b.phased.vcf.gz $id.chr{1..22}.merged.vcf.gz
bcftools index $id.batch$b.phased.vcf.gz
rm $id.chr*.merged.vcf.gz
rm $id.chr*.merged.vcf.gz.csi
vcftools --gzvcf $id.batch$b.phased.vcf.gz --positions /mnt/archgen/users/hringbauer/data/GRG/vcf/pos/1240k.all.sites --out $id.batch$b.phased.1240k --recode
mv $id.batch$b.phased.1240k.recode.vcf $id.batch$b.phased.1240k.vcf
bgzip $id.batch$b.phased.1240k.vcf
bcftools index $id.batch$b.phased.1240k.vcf.gz
# delete BAM file because they are too large, takes too much space
# I used seed when downsampling, so it's easy to reproduce
rm $id.batch$b.bam
rm $id.batch$b.bam.bai

######################## compute phasing error rate on 1240k snps #############################
bcftools view -Oz -o tmp.vcf.gz -r 3 $id.batch$b.phased.1240k.vcf.gz
bcftools index tmp.vcf.gz
echo "I3949" > reheader
bcftools reheader -s reheader tmp.vcf.gz | bcftools view -Oz -o ch3.1240k.vcf.gz
vcftools --gzdiff ch3.1240k.vcf.gz --gzvcf /mnt/archgen/users/yilei/IBDsim/downsample/phasingSwitch/AF.ch3.trioPhased.1240k.vcf.gz --diff-switch-error 
rm ch3.1240k.vcf.gz
rm tmp.vcf.gz
rm tmp.vcf.gz.csi
mv out.diff.indv.switch out.diff.indv.switch.1240k
mv out.diff.switch out.diff.switch.1240k


######################## compute phasing error rate on 1kg snps #############################
bcftools index $id.batch$b.phased.vcf.gz
bcftools view -Oz -o tmp.vcf.gz -r 3 $id.batch$b.phased.vcf.gz
bcftools index tmp.vcf.gz
bcftools reheader -s reheader tmp.vcf.gz | bcftools view -Oz -o ch3.vcf.gz
vcftools --gzdiff ch3.vcf.gz --gzvcf /mnt/archgen/users/yilei/IBDsim/downsample/phasingSwitch/AF.ch3.trioPhased.1kg.vcf.gz --diff-switch-error 
rm ch3.vcf.gz
rm tmp.vcf.gz
rm tmp.vcf.gz.csi
mv out.diff.indv.switch out.diff.indv.switch.1kg
mv out.diff.switch out.diff.switch.1kg

rm reheader

##############################################################################################################

######################## compute phasing error rate on 1240k tv snps #############################
bcftools view -Oz -o tmp.vcf.gz -r 3 $id.batch$b.phased.1240k.vcf.gz
bcftools index tmp.vcf.gz
echo "I3949" > reheader
bcftools reheader -s reheader tmp.vcf.gz | bcftools view -Oz -o ch3.1240k.vcf.gz
vcftools --gzdiff ch3.1240k.vcf.gz --gzvcf /mnt/archgen/users/yilei/IBDsim/downsample/phasingSwitch/AF.ch3.trioPhased.1240k.tv.vcf.gz --diff-switch-error 
rm ch3.1240k.vcf.gz
rm tmp.vcf.gz
rm tmp.vcf.gz.csi
mv out.diff.indv.switch out.diff.indv.switch.1240k.tv
mv out.diff.switch out.diff.switch.1240k.tv


######################## compute phasing error rate on 1kg tv snps #############################
bcftools index $id.batch$b.phased.vcf.gz
bcftools view -Oz -o tmp.vcf.gz -r 3 $id.batch$b.phased.vcf.gz
bcftools index tmp.vcf.gz
bcftools reheader -s reheader tmp.vcf.gz | bcftools view -Oz -o ch3.vcf.gz
vcftools --gzdiff ch3.vcf.gz --gzvcf /mnt/archgen/users/yilei/IBDsim/downsample/phasingSwitch/AF.ch3.trioPhased.1kg.tv.vcf.gz --diff-switch-error 
rm ch3.vcf.gz
rm tmp.vcf.gz
rm tmp.vcf.gz.csi
mv out.diff.indv.switch out.diff.indv.switch.1kg.tv
mv out.diff.switch out.diff.switch.1kg.tv

rm reheader