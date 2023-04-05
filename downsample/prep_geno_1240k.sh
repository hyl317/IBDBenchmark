#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N prep_geno #Name of the command that will be listed in the queue
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
echo "merging results for batch$b"


# # bgzip and index each file
# for id in I2105 I3950 I5273 I5279;
# do
#     for cov in cov2 cov1 cov3over4 cov1over2;
#     do
#         cd ./1240k/$id/$cov/batch$b
#         bgzip $id.batch$b.phased.1240k.vcf
#         bcftools index $id.batch$b.phased.1240k.vcf.gz
#         bgzip $id.batch$b.phased.vcf
#         bcftools index $id.batch$b.phased.vcf.gz
#         cd ../../../../
#     done
# done

############################### now merge vcf of the 5 samples into a single vcf ############################

# base=/mnt/archgen/users/yilei/IBDsim/downsample/1240k
# cd ../callIBD/1240k
# for cov in cov2 cov1 cov3over4 cov1over2;
# do
#     cd ./$cov
#     if [ ! -d batch$b ]; then
#         mkdir batch$b
#     fi
#     cd batch$b

#     rm -r raw
#     rm -r processed_1240k
    
#     bcftools merge $base/I2105/$cov/batch$b/I2105.batch$b.phased.1240k.vcf.gz \
#         $base/I3950/$cov/batch$b/I3950.batch$b.phased.1240k.vcf.gz \
#         $base/I5273/$cov/batch$b/I5273.batch$b.phased.1240k.vcf.gz \
#         $base/I5279/$cov/batch$b/I5279.batch$b.phased.1240k.vcf.gz | bcftools view -O z -o batch$b.merged.1240k.vcf.gz -e 'FORMAT/GT[*]="./."'
#     bcftools index batch$b.merged.1240k.vcf.gz

#     ########## now that vcfs are merged, transform it to the appropriate hdf5 format ##########################
#     mkdir raw
#     mkdir processed_1240k
#     for ch in {1..22};
#     do
#         bcftools view -r $ch -O v -o chr$ch.temp.vcf batch$b.merged.1240k.vcf.gz
#         python3 /mnt/archgen/users/yilei/LAI/phasing_experiment/AfanasievoFamily/downsample/prep_geno.py --vcf chr$ch.temp.vcf --ch $ch --out processed_1240k
#         rm chr$ch.temp.vcf
#     done
#     rm -r ./raw
#     cd ../../
# done

#################################### finally let's run hapBLOCK #######################################

# cd /mnt/archgen/users/yilei/IBDsim/downsample/callIBD/1240k
# for cov in cov2 cov1 cov3over4 cov1over2;
# do
#     cd ./$cov/batch$b
#     python3 /mnt/archgen/users/yilei/IBDsim/downsample/hapBLOCK.py --min_cm1 5
#     cd ../../
# done

########################################################################################################

######################################  run IBIS on imputed genotypes ###################################

# cd /mnt/archgen/users/yilei/IBDsim/downsample/callIBD/1240k
# for cov in cov2 cov1 cov3over4 cov1over2;
# do
#     cd ./$cov/batch$b
#     mkdir IBIS_1240k
#     cd ./IBIS_1240k
#     # Prepare input for IBIS
#     echo "I2105.batch$b.bam    I2105" >> reheader
#     echo "I3950.batch$b.bam    I3950" >> reheader
#     echo "I5273.batch$b.bam    I5273" >> reheader
#     echo "I5279.batch$b.bam    I5279" >> reheader

#     bcftools reheader -s reheader ../batch$b.merged.1240k.vcf.gz | bcftools view -s I2105,I5273,I5279,I3950 -O z -o tmp.vcf.gz
#     plink --vcf tmp.vcf.gz --make-bed --out batch$b
#     /mnt/archgen/users/yilei/bin/ibis/add-map-plink.pl batch$b.bim /mnt/archgen/users/yilei/Data/Hapmap/genetic_map_GRCh37_chr{1..22}.txt > new.bim
#     rm batch$b.bim
#     mv new.bim batch$b.bim
#     rm tmp.vcf.gz
#     rm reheader

#     # run IBIS
#     /mnt/archgen/users/yilei/bin/ibis/ibis -bfile batch$b -noFamID -o batch$b -min_l 5

#     cd ../../../
# done


################################### run hapIBD on Phased Genotypes #########################################
# for cov in cov5 cov2 cov1 cov3over4 cov1over2;
# do
#     cd ./$cov/batch$b
#     mkdir hapIBD
#     cd ./hapIBD
#     # Prepare input for IBIS
#     echo "2:I2105.batch$b.bam    I2105" >> reheader
#     echo "2:I3388.batch$b.bam    I3388" >> reheader
#     echo "2:I3950.batch$b.bam    I3950" >> reheader
#     echo "2:I5273.batch$b.bam    I5273" >> reheader
#     echo "2:I5279.batch$b.bam    I5279" >> reheader

#     bcftools reheader -s reheader ../batch$b.merged.1240k.vcf.gz | bcftools view -s I2105,I3388,I5273,I5279 -O z -o tmp.vcf.gz
#     rm reheader

#     # run hapIBD
#     java -Xmx40g -jar /mnt/archgen/users/yilei/bin/hap-ibd.jar min-output=2 gt=tmp.vcf.gz map=/mnt/archgen/users/yilei/Data/Hapmap/plink.all.GRCh37.map out=batch$b min-seed=1.0 min-extend=0.5 min-mac=1 max-gap=7500
#     # rm tmp.vcf.gz

#     cd ../../../
# done


############################# run IBIS using a more relaxed threshold #################

cd /mnt/archgen/users/yilei/IBDsim/downsample/callIBD/1240k
for cov in cov2 cov1 cov3over4 cov1over2;
do
    cd ./$cov/batch$b
    mkdir IBIS_1240k_relaxed
    cd ./IBIS_1240k_relaxed
    # Prepare input for IBIS
    # echo "I2105.batch$b.bam    I2105" >> reheader
    # echo "I3950.batch$b.bam    I3950" >> reheader
    # echo "I5273.batch$b.bam    I5273" >> reheader
    # echo "I5279.batch$b.bam    I5279" >> reheader

    # bcftools reheader -s reheader ../batch$b.merged.1240k.vcf.gz | bcftools view -s I2105,I5273,I5279,I3950 -O z -o tmp.vcf.gz
    # plink --vcf tmp.vcf.gz --make-bed --out batch$b
    # /mnt/archgen/users/yilei/bin/ibis/add-map-plink.pl batch$b.bim /mnt/archgen/users/yilei/Data/Hapmap/genetic_map_GRCh37_chr{1..22}.txt > new.bim
    # rm batch$b.bim
    # mv new.bim batch$b.bim
    # rm tmp.vcf.gz
    # rm reheader

    # run IBIS
    /mnt/archgen/users/yilei/bin/ibis/ibis -bfile batch$b -noFamID -o batch$b -min_l 5 -er 0.5

    cd ../../../
done
