#!/bin/bash

declare -a missing

for chr in {1..22}
do
    for part in {00..99}
    do
        if [ ! -f "result_coxmeg_truenegcorrected/${chr}_part${part}" ];
        then
            echo ${chr}.part${part}
            bsub -q rerunnable -o output/${chr}.${part}.out -e output/${chr}.${part}.err ~/anaconda3/envs/R/bin/Rscript OS_GWAS_coxmeg.R ${chr} ${part}
        fi
    done
done
