#!/bin/bash

for fold in {1..5}
do
    for chr in {1..22}
    do
        bsub -q rerunnable -e "output/fold${fold}.${chr}.%I.err" -o "output/fold${fold}.${chr}.%I.out" -J "fold${fold}.${chr}.part[1-100]" "Rscript OS_GWAS_coxme_crossvalid.R $chr \$LSB_JOBINDEX $fold"
    done
done
