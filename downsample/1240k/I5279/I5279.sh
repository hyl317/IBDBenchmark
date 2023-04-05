#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N I5279 #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
# -pe smp 2 #needs 8 CPU cores
#$ -l h_vmem=45G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:50:1
#$ -tc 25

b=$SGE_TASK_ID
id="I5279"
echo doing$id
echo batch$b

covs=(cov2 cov1 cov3over4 cov1over2)

targets=(2.0 1.0 0.75 0.5)
path2bam=/mnt/archgen/users/yilei/IBDsim/BAM/1240k/$id.bam
path2bed=/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/target/1240k.all.target.bed

for i in {3..0};
do
    cov=${covs[$i]}
    targetCov=${targets[$i]}
    echo sampling with coverage $targetCov
    frac=`samtools depth -a -b $path2bed -q 30 -Q 30 $path2bam | awk '{print $3}' | ~/stat.pl | awk -v cov=$targetCov '($1=="avg:"){print cov/$2}'`
    echo downsampling frac $frac

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
        bcftools merge --force-samples -O v -o $id.batch$b.chr$ch.tmp.vcf ./GLIMPSE_ligated/$id.batch$b.chr$ch.merged.bcf ./GLIMPSE_phased/$id.batch$b.chr$ch.phased.bcf
    done
    rm chunks.chr*.txt
    rm $id.batch$b.chr*.vcf.gz
    rm $id.batch$b.chr*.vcf.gz.csi
    bcftools concat -O v -o $id.batch$b.phased.vcf $id.batch$b.chr{1..22}.tmp.vcf
    rm $id.batch$b.chr*.tmp.vcf
    vcftools --vcf $id.batch$b.phased.vcf --positions /mnt/archgen/users/hringbauer/data/GRG/vcf/pos/1240k.all.sites --out $id.batch$b.phased.1240k --recode
    mv $id.batch$b.phased.1240k.recode.vcf $id.batch$b.phased.1240k.vcf
    # delete BAM file because they are too large, takes too much space
    # I used seed when downsampling, so it's easy to reproduce
    rm $id.batch$b.bam
    rm $id.batch$b.bam.bai

    cd ../../
done




