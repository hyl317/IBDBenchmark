#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N glimpse_ch3 #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=45G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:100:1
#$ -l 'h=!bionode0[1-9]'
#$ -tc 25


ch=3
b=$SGE_TASK_ID

mkdir batch$b
cd batch$b

iids=(I4893 I4596 I1583 I2978 I5838 I1507 I2861 I2520 I3758 I5077 I0708 I5233 I3123)
covs=(5.140 4.271 9.208 4.476 5.288 4.704 4.771 5.758 4.699 3.209 7.473 2.728 5.120)
targetCoverage=1.5

for i in {1..13}
do
    i=$(($i-1))
    iid=${iids[$i]}
    cov=${covs[$i]}
    fraction=$(echo "scale=8;$targetCoverage/$cov" | bc)
    echo downsampling $iid to target coverage $targetCoverage with sampling fraction $fraction
    path2BAM="/mnt/archgen/users/yilei/IBDsim/BAM/1240k/$iid.bam"
    samtools view --subsample-seed $b --subsample $fraction -b -o $iid.chr$ch.bam -h $path2BAM $ch
    samtools index $iid.chr$ch.bam

    #prepare GL file (with bcftools) for GLIMPSE
    VCF=/mnt/archgen/users/yilei/Data/1000G_release/reference_panel/1000GP.chr$ch.sites.vcf.gz
    TSV=/mnt/archgen/users/yilei/Data/1000G_release/reference_panel/1000GP.chr$ch.sites.tsv.gz
    REFGEN=/mnt/archgen/users/yilei/Data/1000G_release/raw_data/hs37d5.fa
    OUT=$iid.chr$ch.vcf.gz

    #bcftools mpileup -f ${REFGEN} -I -B -a 'FORMAT/DP' -r $ch -T ${VCF} -q 30 -Q 30 $iid.chr$ch.bam -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
    echo "SM    $iid" > reheader.txt
    bcftools mpileup -f ${REFGEN} -I -B -a 'FORMAT/DP' -r $ch -T ${VCF} -q 30 -Q 30 $iid.chr$ch.bam -Ou | bcftools call -Aim -C alleles -T ${TSV} -Ou | bcftools reheader -s reheader.txt | bcftools view -Oz -o ${OUT}
    bcftools index -f ${OUT}
    rm -r reheader.txt

    # ####### delete bam file to save storage space
    rm $iid.chr$ch.bam
    rm $iid.chr$ch.bam.bai
    # ###################################################

    # running GLIMPSE
    path2glimpse="/mnt/archgen/users/yilei/bin//GLIMPSE-glimpse1"
    VCF=$iid.chr$ch.vcf.gz
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
            OUT=GLIMPSE_imputed/$iid.chr$ch.${ID}.bcf
            $path2glimpse/phase/bin/GLIMPSE_phase --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
            bcftools index -f ${OUT}
        done < chunks.chr$ch.txt

    if [ ! -d GLIMPSE_ligated ]; then
        mkdir GLIMPSE_ligated
    fi
    LST=GLIMPSE_ligated/list.chr$ch.txt
    ls GLIMPSE_imputed/$iid.chr$ch.*.bcf > ${LST}
    OUT=GLIMPSE_ligated/$iid.chr$ch.merged.bcf
    $path2glimpse/ligate/bin/GLIMPSE_ligate --input ${LST} --output $OUT
    bcftools index -f ${OUT}
    rm -r ./GLIMPSE_imputed

    if [ ! -d GLIMPSE_phased ]; then
        mkdir GLIMPSE_phased
    fi
    VCF=GLIMPSE_ligated/$iid.chr$ch.merged.bcf
    OUT=GLIMPSE_phased/$iid.chr$ch.phased.bcf
    $path2glimpse/sample/bin/GLIMPSE_sample --input ${VCF} --solve --output ${OUT}
    bcftools index -f ${OUT}

    bcftools annotate -a ./GLIMPSE_ligated/$iid.chr$ch.merged.bcf -c FORMAT/GP ./GLIMPSE_phased/$iid.chr$ch.phased.bcf -Oz -o $iid.chr$ch.merged.vcf.gz
    bcftools index $iid.chr$ch.merged.vcf.gz

    #### final clean-up
    rm $iid.chr$ch.vcf.gz
    rm $iid.chr$ch.vcf.gz.csi
    rm -r ./GLIMPSE_ligated
    rm -r ./GLIMPSE_phased
    rm chunks.chr$ch.txt
done

bcftools merge -Oz -o batch$b.all.chr$ch.vcf.gz I4893.chr$ch.merged.vcf.gz I4596.chr$ch.merged.vcf.gz \
        I1583.chr$ch.merged.vcf.gz I2978.chr$ch.merged.vcf.gz I5838.chr$ch.merged.vcf.gz \
        I1507.chr$ch.merged.vcf.gz I2861.chr$ch.merged.vcf.gz I2520.chr$ch.merged.vcf.gz \
        I3758.chr$ch.merged.vcf.gz I5077.chr$ch.merged.vcf.gz I0708.chr$ch.merged.vcf.gz \
        I5233.chr$ch.merged.vcf.gz I3123.chr$ch.merged.vcf.gz

python3 /mnt/archgen/users/yilei/IBDsim/FP/downsample/run_ancIBD.py --vcf batch$b.all.chr$ch.vcf.gz --prefix batch$b
rm *.merged.vcf.gz
rm *.merged.vcf.gz.csi