#!/bin/bash

for chr in {1..22}
do
    for part in {00..99}
    do
        bsub -q rerunnable -e output/${chr}.${part}.err -o output/${chr}.${part}.out -n 1 ~/anaconda3/envs/R/bin/Rscript OS_GWAS_coxmeg.R $chr $part
    done
done
