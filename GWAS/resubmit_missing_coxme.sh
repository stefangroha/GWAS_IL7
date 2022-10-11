#!/bin/bash

for chr in {1..22}
do
    for part in {00..99}
    do
        if [ ! -f "result_coxme/${chr}.part${part}" ];
        then
            echo ${chr}.part${part}
            bsub -q long -o output/${chr}.${part}.out -e output/${chr}.${part}.err Rscript OS_GWAS_coxme.R ${chr} ${part}
        fi
    done
done
