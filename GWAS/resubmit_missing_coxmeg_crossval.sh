#!/bin/bash

fold=$1

for chr in {1..22}
do
    for part in {00..99}
    do
        if [ ! -f "./result_coxmeg_crossval/${fold}/${chr}_part${part}" ];
        then
            echo ${chr}_part${part}
            part_num=$(echo $part | sed 's/^0*//')
            num_part=$((part_num + 1))
            bsub -q rerunnable -e "output/fold${fold}.${chr}.${part}.err" -o "output/fold${fold}.${chr}.${part}.out" Rscript OS_GWAS_coxmeg_crossvalid.R $chr $num_part $fold
        fi
    done
done
