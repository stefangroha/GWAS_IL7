#!/bin/bash

for chr in {1..22}
do
    for part in {00..99}
    do
        bsub -q long -n 1 -R 'rusage[mem=8000]' -e output/${chr}.${part}.err -o output/${chr}.${part}.out -n 1 Rscript OS_GWAS_coxme.R $chr $part
    done
done
