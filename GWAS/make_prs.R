suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)
fold = args[[1]]

result <- data.frame()

files <- list.files(paste("result_coxmeg_crossval/",fold,"/",sep=""))
pb <- txtProgressBar(min = 1, max = length(files), style = 3)
count <- 1
for (file in files) {
    res <- fread(paste("result_coxmeg_crossval/",file,sep=""))
    result <- rbind(result,res)
    setTxtProgressBar(pb, count)  
    count <- count + 1
}

snptrans <- fread("cat /data/gusev/PROFILE/STITCH/1000GP_Phase3_chr*.snpid")
snptrans <- snptrans %>% select(V1,V3) %>% mutate(V1=gsub(":","_",V1)) %>% rename(allele=V1,SNP=V3)

snpids <- data.frame()
count <- 1
for (chr in 1:22) {
        snpids <- rbind(snpids,fread(paste("/data/gusev/PROFILE/STITCH/1000GP_Phase3_chr",chr,".snpid",sep="")))
}
names(snpids) <- c("chr_pos","CHR","SNP","POS","A1","A2")

prs <- result %>% mutate(CHR=as.integer(gsub("([0-9]{1,2})_.*","\\1",allele)),POS=as.integer(gsub("[0-9]{1,2}_(.*)","\\1",allele)))
prs <- left_join(prs,snpids)

ids <- fread(paste("./crossval_folds/",fold,"_train",sep=""))
print(dim(ids)[[1]])
# prs <- prs %>% select(CHR,POS,SNP,A1,A2,beta,p) %>% rename(PVAL=p,EFF=beta,RS=SNP) %>% mutate(NCOL=c(1666))
