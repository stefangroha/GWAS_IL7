#!/bin/bash

fold=$1

for chr in {1..22}
do
    for part in {0..99}
    do
        if [ ! -f "./cross_validation/coxme/fold${fold}/${chr}.part${part}" ];
        then
            echo ${chr}.part${part}
            num_part=$((part + 1))
            bsub -q rerunnable -e "output/fold${fold}.${chr}.${part}.err" -o "output/fold${fold}.${chr}.${part}.out" Rscript OS_GWAS_coxme_crossvalid.R $chr $num_part $fold
        fi
    done
done
