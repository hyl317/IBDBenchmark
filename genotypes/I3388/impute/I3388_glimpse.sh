#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N glimipse_I3388 #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=45G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:22:1

ch=$SGE_TASK_ID
echo imputing ch$ch

if [ ! -d ch$ch ]; then
    mkdir ch$ch
fi

cd ch$ch

id=I3388
path2BAM="/mnt/archgen/users/yilei/IBDsim/BAM/$id.bam"


#prepare GL file (with bcftools) for GLIMPSE
VCF=/mnt/archgen/users/yilei/Data/1000G_release/reference_panel/1000GP.chr$ch.sites.vcf.gz
TSV=/mnt/archgen/users/yilei/Data/1000G_release/reference_panel/1000GP.chr$ch.sites.tsv.gz
REFGEN=/mnt/archgen/users/yilei/Data/1000G_release/raw_data/hs37d5.fa
OUT=$id.chr$ch.vcf.gz

bcftools mpileup -f ${REFGEN} -I -B -a 'FORMAT/DP' -r $ch -T ${VCF} -q 30 -Q 30 ${path2BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
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

