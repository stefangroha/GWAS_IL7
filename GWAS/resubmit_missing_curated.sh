#!/bin/bash

declare -a missing

for chr in {1..22}
do
    for part in {00..99}
    do
        if [ ! -f "result_curated_coxmeg/${chr}_part${part}" ];
        then
            echo ${chr}.part${part}
            bsub -q rerunnable -o output/${chr}.${part}.out -e output/${chr}.${part}.err Rscript OS_GWAS_curated_coxmeg.R ${chr} ${part}
        fi
    done
done
