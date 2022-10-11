suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(caret))

set.seed(42)

sample_names <- c(t(suppressWarnings(read.table("/data/gusev/PROFILE/GWAS_Surv/IO_COHORT/KINSHIP/kinship.matrix.rel.id"))))

folds <- createFolds(sample_names, k = 5)
sample_folds <- lapply(folds, function(ind, dat) dat[ind], dat = sample_names)

for (i in 1:5) {
    write(c(unlist(sample_folds[setdiff(c(1:5),c(i))])),paste("crossval_folds/",i,"_train",sep=""))
    write(sample_folds[[i]],paste("crossval_folds/",i,"_val",sep=""))
}
