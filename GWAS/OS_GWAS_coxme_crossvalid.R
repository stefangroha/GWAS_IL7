suppressPackageStartupMessages(library(mltools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(coxme))

# input 
args = commandArgs(trailingOnly=TRUE)
chr_number = args[[1]]
part = as.integer(args[[2]])-1
fold = args[[3]]

print(paste("chromosome: ",chr_number,sep=""))
print(paste("part: ",part,sep=""))
print(paste("fold: ",fold,sep=""))

pid2bl <- fread("../../../../OncDRS_data/ALL/REQ_AG618_103317_1_GENOMIC_SPECIMEN.csv",fill=TRUE) %>% select(PATIENT_ID,SAMPLE_ACCESSION_NBR) %>% rename(BLNUM=SAMPLE_ACCESSION_NBR) %>% unique
bl2cbio <- fread("../../../../OncDRS_data/ALL/sample_id_to_cbio.txt") %>% select(SampleID,cBio_PatientID,cBio_SampleID) %>% rename(BLNUM=SampleID,CBIO_PATIENT=cBio_PatientID,CBIO_SAMPLE=cBio_SampleID) %>% unique
bl2cbio <- full_join(fread("../../../../OncDRS_data/LUNG/raw/PTID_CBIO.csv") %>% select(BLNUM,CBIO_PATIENT,CBIO_SAMPLE),bl2cbio,by=c("BLNUM","CBIO_PATIENT","CBIO_SAMPLE")) %>% distinct
pid2mrn <- fread("../../../../OncDRS_data/ALL/REQ_AG618_103317_1_PT_INFO_STATUS_REGISTRATION.csv",fill=TRUE) %>% select(PATIENT_ID,DFCI_MRN) %>% rename(MRN=DFCI_MRN) %>% unique
temp <- right_join(bl2cbio,pid2bl,by="BLNUM") %>% group_by(PATIENT_ID) %>% summarise(CBIO_PATIENT=as.character(ifelse(all(is.na(CBIO_PATIENT)),NA,first(CBIO_PATIENT[!is.na(CBIO_PATIENT)]))))
translate_pid <- merge(temp,pid2mrn,by="PATIENT_ID") %>% select(PATIENT_ID,CBIO_PATIENT,MRN) %>% unique
translate_sid <- right_join(bl2cbio,pid2bl,by="BLNUM") %>% select(PATIENT_ID,BLNUM,CBIO_PATIENT,CBIO_SAMPLE) %>% unique

OS <- fread("../../OS/OS.csv",na.string="")
OS <- OS %>% mutate(LAST_DATE=as.Date(LAST_DATE),IO_START=as.Date(IO_START),IO_STOP=as.Date(IO_STOP))

OS$PANEL_VERSION <- factor(OS$PANEL_VERSION)
OS$AGE_AT_TREATMENTSTART_BINNED <- factor(OS$AGE_AT_TREATMENTSTART_BINNED)
OS$ANCESTRY_SELFREPORTED <- factor(OS$ANCESTRY_SELFREPORTED)
OS$ANCESTRY_SELFREPORTED <- relevel(OS$ANCESTRY_SELFREPORTED,"White")
OS$ANCESTRY_GENETICS <- factor(OS$ANCESTRY_GENETICS)
OS$ANCESTRY_GENETICS <- relevel(OS$ANCESTRY_GENETICS,"White")
OS$TMB_binned <- factor(OS$TMB_binned)
OS$TMB_binned <- relevel(OS$TMB_binned,"0")
OS$PROFILE_AFTER_TREATMENT <- factor(OS$PROFILE_AFTER_TREATMENT)
OS$TREATMENT_CUM <- factor(OS$TREATMENT_CUM)
levels(OS$TREATMENT_CUM) <- gsub("/","_",levels(OS$TREATMENT_CUM))
OS$TREATMENT_CUM <- factor(OS$TREATMENT_CUM)
OS$TREATMENT_CUM <- relevel(OS$TREATMENT_CUM,"PD1_PDL1")
OS$MEDICATION_CUM <- factor(OS$MEDICATION_CUM)
OS$MEDICATION_CUM <- relevel(OS$MEDICATION_CUM,"pembrolizumab")
OS$GENDER <- factor(OS$GENDER)
OS$PANEL_VERSION <- factor(OS$PANEL_VERSION)
OS$PANEL_VERSION <- relevel(OS$PANEL_VERSION,"3")
OS$START_YEAR <- factor(OS$START_YEAR)
levels(OS$START_YEAR) <- gsub(">","after",gsub("<","before",gsub("/","_",levels(OS$START_YEAR))))
OS$START_YEAR <- relevel(OS$START_YEAR,"2016_2017")
OS$PRIOR_TREATMENT <- factor(OS$PRIOR_TREATMENT)
levels(OS$PRIOR_TREATMENT) <- gsub("\\(","",gsub("\\)","",gsub(" ","_",levels(OS$PRIOR_TREATMENT))))
OS$TPLAN_GOAL <- factor(OS$TPLAN_GOAL)
OS$TPLAN_GOAL <- relevel(OS$TPLAN_GOAL,"PALLIATIVE")
OS$STAGE <- factor(OS$STAGE)
OS$STAGE <- relevel(OS$STAGE,"STAGE1_2")

allele_in <- fread(paste("../../../../GWAS/IO_COHORT/DOSAGES/",chr_number,".part",part,sep=""),data.table=F) %>% rename(allele=V1)

patients <- gsub("_S[0-9]*$","",(fread("train_fivefold.cross",na.strings = "") %>% select(matches(paste("fold_",fold,sep=""))) %>% drop_na)[[paste("fold_",fold,sep="")]])

OS_C <- OS %>% filter(CBIO_PATIENT %in% patients)
vars=c("AGE_AT_TREATMENTSTART","GENDER","START_YEAR","PROFILE_AFTER_TREATMENT","TREATMENT_CUM","PANEL_VERSION","LINE_OF_TREATMENT","ANCESTRY_PC_NS","ANCESTRY_PC_ASH")
OS_C <- OS_C %>% select(MRN,CBIO_PATIENT,tstart,tstop,event,CANCER_TYPE,one_of(vars)) %>% drop_na
mm <- model.matrix(as.formula(paste("~ ",paste(vars,collapse=" + "),sep="")), OS_C)  
OS_C <- OS_C %>% dplyr::select(-one_of(vars))
OS_C <- cbind(OS_C,mm)
vars <- c(unlist(names(OS_C)[grepl(paste(vars,collapse="|"),names(OS_C))]))
# merge with first variant to remove samples which are NaN
df <- data.frame(CBIO_SAMPLE=names(allele_in),dosage=as.vector(as.matrix(allele_in)[1,]))
df <- df %>% mutate(CBIO_PATIENT=gsub("_S1","",CBIO_SAMPLE)) %>% dplyr::select(-CBIO_SAMPLE)
OS_C <- merge(OS_C,df) %>% dplyr::select(-dosage)

prophazard_pval_cutoff <- 0.01
formula <- paste("Surv(tstart,tstop,event) ~ strata(CANCER_TYPE) + ",paste(vars,collapse=" + "),sep="")
fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
nonfinite_vars <- names(fit$coefficients)[unlist(is.na(fit$coefficients))]
vars <- vars[!vars %in% nonfinite_vars]
formula <- paste("Surv(tstart,tstop,event) ~ strata(CANCER_TYPE) + ",paste(vars,collapse=" + "),sep="")
fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
strat_vars <- c()
while(length(vars[unlist(lapply(vars,function(x) grepl(x,paste(unlist((setDT(data.frame(suppressWarnings(cox.zph(fit))$table), keep.rownames = TRUE)[] %>% filter(p<prophazard_pval_cutoff))$rn),collapse=""))))] > 0)) {
    df <- data.frame(suppressWarnings(cox.zph(fit))$table)
    strat_vars <- c(strat_vars,vars[unlist(lapply(vars,function(x) grepl(x,paste(unlist((setDT(df, keep.rownames = TRUE)[] %>% filter(p<prophazard_pval_cutoff))$rn),collapse=""))))])
    if (length(strat_vars)==0) {
        break
    }
    vars <- vars[!vars %in% strat_vars]
    cont_strat_vars <- strat_vars[(sapply(OS_C[strat_vars], is.double))] 
    if (length(cont_strat_vars)>0) {
        strat_vars <- strat_vars[!strat_vars %in% cont_strat_vars]
        new_strat <- data.table(OS_C[cont_strat_vars])
        for (csv in cont_strat_vars) {
            new_strat[,(paste(csv,"_binned",sep="")):=findInterval(get(csv),quantile(get(csv),prob=c(1/3.,2/3.)),rightmost.closed=TRUE)]
        }
        new_strat <- new_strat %>% mutate_all(as.factor)                                             
        x <- apply(model.matrix(as.formula(paste("~ ",paste(paste(cont_strat_vars,"_binned",sep=""),collapse=" + "),sep="")), new_strat),2,as.integer)         
        x <- x[,colnames(x)!="(Intercept)"]
        vars <- c(vars,colnames(x))
        OS_C <- cbind(OS_C,x)
    }
    if (length(strat_vars)==0) {
        formula <- paste("Surv(tstart,tstop,event) ~ strata(CANCER_TYPE) + ",paste(vars,collapse=" + "),sep="")
    } else {
        formula <- paste("Surv(tstart,tstop,event) ~ strata(CANCER_TYPE) + ",paste(vars,collapse=" + ")," + strata(",paste(strat_vars,collapse=","),")",sep="")
    } 
    fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
}                                          
if (length(strat_vars)>0) {
    form <- paste("Surv(tstart,tstop,event) ~ strata(CANCER_TYPE) + (1|MRN) + dosage + ",paste(vars,collapse=" + ")," + strata(",paste(strat_vars,collapse=","),")",sep="")
} else {
    form <- paste("Surv(tstart,tstop,event) ~ strata(CANCER_TYPE) + (1|MRN) + dosage + ",paste(vars,collapse=" + "),sep="")
}                                                  
                                                  
result_cancer_chrom <- data.frame()
len_of_loop <- length(rownames(allele_in))
# len_of_loop <- 100
pb <- txtProgressBar(min = 1, max = len_of_loop, style = 3)

for (i in 1:len_of_loop) {
    df <- setDT(data.frame(t(allele_in[i,])),keep.rownames = T)[] %>% filter(rn!="allele")
    df <- df %>% rename(CBIO_SAMPLE=rn,dosage=contains("X")) %>% mutate(dosage=as.numeric(as.character(dosage)))
    df <- df %>% mutate(CBIO_PATIENT=gsub("_S1","",CBIO_SAMPLE)) %>% dplyr::select(-CBIO_SAMPLE)
    OS_df <- suppressWarnings(inner_join(OS_C,df,by="CBIO_PATIENT"))
    fit <- coxme(as.formula(form),OS_df)
    haz <- fit$coefficients[["dosage"]]
    var <- diag(vcov(fit))[[1]]
    p_wald <- 2*pnorm(-abs(haz/sqrt(var)))
    
    temp <- data.frame(allele=rownames(allele_in)[[i]],loghazard=haz,var=var,p_wald=p_wald)
    result_cancer_chrom <- rbind(result_cancer_chrom,temp)
    setTxtProgressBar(pb, i)    
}
fwrite(result_cancer_chrom,paste("./cross_validation/coxme/fold",fold,"/",chr_number,".part",part,sep=""))
