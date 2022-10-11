#!/bin/bash

for cross in {1..5}
do
    echo crossval: $cross
    Rscript make_prs.R $cross
done
