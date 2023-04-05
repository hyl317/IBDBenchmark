#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N IBDseq #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -pe smp 12 #needs 8 CPU cores
#$ -l h_vmem=125G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:50:1


b=$SGE_TASK_ID
echo "running IBD seq for batch$b"




base=/mnt/archgen/users/yilei/IBDsim/downsample/wgs
cd callIBD/wgs

for cov in cov5 cov2 cov1 cov3over4 cov1over2 cov1over4 cov1over10;
do
    cd ./$cov
    if [ ! -d batch$b ]; then
        mkdir batch$b
    fi
    cd batch$b

    if [ ! -d IBDseq_tvOnly ]; then
        mkdir IBDseq_tvOnly
    fi
    cd IBDseq_tvOnly

    zcat ../IBDseq/batch$b.merged.ibdseq.vcf.gz | awk -F '\t' '($1~/^#/) || !(($4=="A"&&$5=="G")||($4=="G"&&$5=="A")||($4=="C"&&$5=="T")||($4=="T"&&$5=="C"))' > batch$b.merged.ibdseq.tvOnly.vcf
    bgzip batch$b.merged.ibdseq.tvOnly.vcf
    bcftools index batch$b.merged.ibdseq.tvOnly.vcf.gz

    for ch in {1..22};
    do
        java -Xmx110g -jar /mnt/archgen/users/yilei/bin/ibdseq.r1206.jar gt=batch$b.merged.ibdseq.tvOnly.vcf.gz out=batch$b.IBDseq.ch$ch chrom=$ch nthreads=12
    done

    cd ../../../
done


####### clean up
cd callIBD/wgs
for cov in cov5 cov2 cov1 cov3over4 cov1over2 cov1over4 cov1over10;
do
    cd ./$cov/batch$b/IBDseq_tvOnly
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
done