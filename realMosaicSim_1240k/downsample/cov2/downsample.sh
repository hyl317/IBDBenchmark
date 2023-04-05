#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N glimipse_ch3 #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=45G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:49:1

ch=3
i=$SGE_TASK_ID
id="$(sed "${i}q;d" /mnt/archgen/users/yilei/IBDsim/BAM/1240k/iid.list)"

mkdir -p $id/ch3/
cd $id/ch3

path2BAM="/mnt/archgen/users/yilei/IBDsim/BAM/1240k/$id.bam"
samtools index $path2BAM

#### downsample original BAM file to desired coverage

### determine downsampling frac

path2BED="/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/target/1240kChr3.target.bed"
frac=`samtools depth -a -b $path2BED -q 30 -Q 30 -r $ch $path2BAM | awk '{print $3}' | ~/stat.pl | awk '($1=="avg:"){print 2/$2}'`
echo fraction of downsampling is $frac


samtools view --subsample-seed 1 --subsample $frac -b -o $id.downsampled.ch3.bam $path2BAM $ch
samtools index $id.downsampled.ch3.bam

#prepare GL file (with bcftools) for GLIMPSE
VCF=/mnt/archgen/users/yilei/Data/1000G_release/reference_panel/1000GP.chr$ch.sites.vcf.gz
TSV=/mnt/archgen/users/yilei/Data/1000G_release/reference_panel/1000GP.chr$ch.sites.tsv.gz
REFGEN=/mnt/archgen/users/yilei/Data/1000G_release/raw_data/hs37d5.fa
OUT=$id.chr$ch.vcf.gz

bcftools mpileup -f ${REFGEN} -I -B -a 'FORMAT/DP' -r $ch -T ${VCF} -q 30 -Q 30 $id.downsampled.ch3.bam -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
bcftools index -f ${OUT}

# running GLIMPSE
path2glimpse="/mnt/archgen/users/yilei/bin/GLIMPSE"
VCF=$id.chr$ch.vcf.gz
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
        OUT=GLIMPSE_imputed/$id.chr$ch.${ID}.bcf
        $path2glimpse/phase/bin/GLIMPSE_phase --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
        bcftools index -f ${OUT}
    done < chunks.chr$ch.txt

if [ ! -d GLIMPSE_ligated ]; then
    mkdir GLIMPSE_ligated
fi
LST=GLIMPSE_ligated/list.chr$ch.txt
ls GLIMPSE_imputed/$id.chr$ch.*.bcf > ${LST}
OUT=GLIMPSE_ligated/$id.chr$ch.merged.bcf
$path2glimpse/ligate/bin/GLIMPSE_ligate --input ${LST} --output $OUT
bcftools index -f ${OUT}
rm -r ./GLIMPSE_imputed

if [ ! -d GLIMPSE_phased ]; then
    mkdir GLIMPSE_phased
fi
VCF=GLIMPSE_ligated/$id.chr$ch.merged.bcf
OUT=GLIMPSE_phased/$id.chr$ch.phased.bcf
$path2glimpse/sample/bin/GLIMPSE_sample --input ${VCF} --solve --output ${OUT}
bcftools index -f ${OUT}

path2sites="/mnt/archgen/users/hringbauer/data/GRG/vcf/pos/chr$ch.sites"
vcftools --bcf $VCF --positions $path2sites --recode --out ./GLIMPSE_ligated/$id.chr3.1240k; 
mv ./GLIMPSE_ligated/$id.chr3.1240k.recode.vcf ./GLIMPSE_ligated/$id.chr3.1240k.vcf