#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N ibis #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
#$ -l h_vmem=50G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID

DP=15
base=/mnt/archgen/users/yilei/IBDsim/genotypes

for id1 in I2105 I3950 I5273 I5279;
do
    for id2 in I2105 I3950 I5273 I5279;
    do
        if [[ "$id1" == "$id2" || "$id1" > "$id2" ]]; then
            continue
        else
            echo "finding IBD between $id1 and $id2"
            if [ ! -d "${id1}_${id2}" ]; then
                mkdir "${id1}_${id2}"
            fi
            cd "${id1}_${id2}"
            
            bcftools merge -O u $base/$id1/$id1.minDP$DP.all.vcf.gz \
                $base/$id2/$id2.minDP$DP.all.vcf.gz | bcftools view -O v -g ^miss -i 'MAC>=1'\
                | awk -F '\t' '($1~/^#/) || !(($4=="A"&&$5=="G")||($4=="G"&&$5=="A")||($4=="C"&&$5=="T")||($4=="T"&&$5=="C"))' > "${id1}_${id2}.vcf"
            bgzip "${id1}_${id2}.vcf"
            plink --vcf "${id1}_${id2}.vcf.gz" --make-bed --out "${id1}_${id2}"

            # prepare input for ibis and run ibis
            /mnt/archgen/users/yilei/bin/ibis/add-map-plink.pl "${id1}_${id2}.bim" /mnt/archgen/users/yilei/Data/Hapmap/genetic_map_GRCh37_chr{1..22}.txt > new.bim
            rm "${id1}_${id2}.bim"
            mv new.bim "${id1}_${id2}.bim"
            /mnt/archgen/users/yilei/bin/ibis/ibis -bfile "${id1}_${id2}" -noFamID -o "${id1}_${id2}" -mL 4


            cd ../
        fi
    done
done