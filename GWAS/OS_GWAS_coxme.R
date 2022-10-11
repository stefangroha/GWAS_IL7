suppressPackageStartupMessages(library(mltools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(coxme))
suppressPackageStartupMessages(library(mstate))

# input 
args = commandArgs(trailingOnly=TRUE)
chr_number = args[[1]]
part = args[[2]]

make_multistate <- function(OS_temp) {
    temp <- merge(OS_temp,fread("../adverse_events_predicted.csv",na.strings="") %>% select(MRN,ADV_EVENT,ADV_EVENT_DT),by=c("MRN")) %>% mutate(irAE_time=as.integer(ifelse(ADV_EVENT,as.integer(difftime(ADV_EVENT_DT,IO_START,unit="days")),tstop)),irAE=as.numeric(ifelse(ADV_EVENT,1,0))) %>% select(-ADV_EVENT_DT,-ADV_EVENT)
    temp <- temp %>% mutate(irAE_time=as.integer(ifelse(irAE_time==0,1,irAE_time)),tstop=as.integer(ifelse(tstop<=irAE_time,irAE_time,tstop)))
    tmat <- transMat(x = list(c(2,3), c(3), c()), names = c("Treatment", "AE", "Death"))
    OS_mult <- msprep(data = temp, trans = tmat, time = c(NA, "irAE_time", "tstop"), status = c(NA, "irAE", "event"), id="MRN")
    OS_mult <- merge(OS_mult,temp,by="MRN") %>% rowwise %>% mutate(Tstart=max(tstart,Tstart),time=Tstop-Tstart) %>% select(-irAE_time,-tstart,-tstop,-event,-irAE)
    OS_mult <- OS_mult %>% filter(Tstop>Tstart)

    attributes(OS_mult)$class <- c("msdata","data.frame")
    attributes(OS_mult)$trans <- tmat

    print(events(OS_mult))
    covs <- names(OS_mult %>% select(-PBP,-MRN,-from,-to,-trans,-Tstart,-Tstop,-time,-status,-IO_START))
    OS_mult <- expand.covs(OS_mult, covs, longnames = T)
    return(OS_mult)
}

pid2bl <- fread("/PHShome/szg99/immuno/Profile/OncDRS_data/ALL/GENOMIC_SPECIMEN.csv",fill=TRUE) %>% select(PATIENT_ID,SAMPLE_ACCESSION_NBR) %>% rename(BLNUM=SAMPLE_ACCESSION_NBR) %>% unique
bl2cbio <- fread("/PHShome/szg99/immuno/Profile/OncDRS_data/ALL/sample_id_to_cbio.txt") %>% select(SampleID,cBio_PatientID,cBio_SampleID) %>% rename(BLNUM=SampleID,CBIO_PATIENT=cBio_PatientID,CBIO_SAMPLE=cBio_SampleID) %>% unique
bl2cbio <- full_join(fread("/PHShome/szg99/immuno/Profile/OncDRS_data/LUNG/raw/PTID_CBIO.csv") %>% select(BLNUM,CBIO_PATIENT,CBIO_SAMPLE),bl2cbio,by=c("BLNUM","CBIO_PATIENT","CBIO_SAMPLE")) %>% distinct
pid2mrn <- fread("/PHShome/szg99/immuno/Profile/OncDRS_data/ALL/PT_INFO_STATUS_REGISTRATION.csv",fill=TRUE) %>% select(PATIENT_ID,DFCI_MRN) %>% rename(MRN=DFCI_MRN) %>% unique
pbp2bl <- fread("/data/gusev/PROFILE/2020_04/gusev-pseudonymize-output.txt") %>% rename(BLNUM=bl_number,PBP=pbp_number) %>% select(BLNUM,PBP)
temp <- right_join(bl2cbio,pid2bl,by="BLNUM") %>% group_by(PATIENT_ID) %>% summarise(CBIO_PATIENT=as.character(ifelse(all(is.na(CBIO_PATIENT)),NA,first(CBIO_PATIENT[!is.na(CBIO_PATIENT)]))))
translate_pid <- merge(temp,pid2mrn,by="PATIENT_ID") %>% select(PATIENT_ID,CBIO_PATIENT,MRN) %>% unique
translate_sid <- full_join(pbp2bl,right_join(bl2cbio,pid2bl,by="BLNUM"),by="BLNUM") %>% select(PATIENT_ID,BLNUM,PBP,CBIO_PATIENT,CBIO_SAMPLE) %>% unique

OS <- fread("../../OS/OS.csv",na.string="")
OS <- OS %>% mutate(LAST_DATE=as.Date(LAST_DATE),IO_START=as.Date(IO_START))

OS$PANEL_VERSION <- factor(OS$PANEL_VERSION)
OS$AGE_AT_TREATMENTSTART_BINNED <- factor(OS$AGE_AT_TREATMENTSTART_BINNED)
OS$ANCESTRY_SELFREPORTED <- factor(OS$ANCESTRY_SELFREPORTED)
OS$ANCESTRY_SELFREPORTED <- relevel(OS$ANCESTRY_SELFREPORTED,"White")
OS$ANCESTRY_GENETICS <- factor(OS$ANCESTRY_GENETICS)
OS$ANCESTRY_GENETICS <- relevel(OS$ANCESTRY_GENETICS,"EUR")
OS$TMB_binned <- factor(OS$TMB_binned)
OS$TMB_binned <- relevel(OS$TMB_binned,"0")
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
OS$START_YEAR <- relevel(OS$START_YEAR,"2016")
OS$PRIOR_TREATMENT <- factor(OS$PRIOR_TREATMENT)
levels(OS$PRIOR_TREATMENT) <- gsub("\\(","",gsub("\\)","",gsub(" ","_",levels(OS$PRIOR_TREATMENT))))
OS$TPLAN_GOAL <- factor(OS$TPLAN_GOAL)
OS$TPLAN_GOAL <- relevel(OS$TPLAN_GOAL,"PALLIATIVE")
OS$STAGE <- factor(OS$STAGE)
OS$STAGE <- relevel(OS$STAGE,"STAGE1_2")
OS$CANCER_TYPE <- gsub(" |-","_",OS$CANCER_TYPE)

allele_in <- fread(paste("/data/gusev/PROFILE/GWAS_Surv/IO_COHORT/DOSAGES_split/PROF_2020_04.",chr_number,".QC.part",part,".traw",sep=""),data.table=F)
allele_in <- allele_in %>% dplyr::select(-CHR,-`(C)M`,-POS,-COUNTED,-ALT)

OS_C <- OS# %>% filter(CANCER_TYPE==gsub("_"," ",cancer))
vars=c("CANCER_TYPE","AGE_AT_TREATMENTSTART","GENDER","START_YEAR","PROFILE_AFTER_TREATMENT","TREATMENT_CUM","PANEL_VERSION","LINE_OF_TREATMENT","BIOPSY_SITE_TYPE","ANCESTRY_PC_NS","ANCESTRY_PC_ASH")
OS_C <- OS_C %>% select(MRN,tstart,tstop,event,IO_START,CANCER_TYPE,one_of(vars)) %>% drop_na
mm <- model.matrix(as.formula(paste("~ ",paste(vars,collapse=" + "),"-1",sep="")), OS_C)  
OS_C <- OS_C %>% dplyr::select(-one_of(vars)) #%>% cbind(x)
OS_C <- cbind(OS_C,mm)
vars <- c(unlist(names(OS_C)[grepl(paste(vars,collapse="|"),names(OS_C))]))
# merge with first variant to remove samples which are NaN
df <- data.frame(PBP=names(allele_in)[names(allele_in)!="SNP"],dosage=as.vector(as.matrix(allele_in)[1,])[-1])
df <- df %>% mutate(PBP=gsub("^0_","",PBP))
OS_C <- inner_join(inner_join(inner_join(OS_C,translate_pid %>% dplyr::select(MRN,PATIENT_ID),by="MRN"),translate_sid %>% dplyr::select(PATIENT_ID,PBP),by="PATIENT_ID") %>% dplyr::select(-PATIENT_ID),df) %>% dplyr::select(-dosage)
OS_C <- OS_C %>% select(-BIOPSY_SITE_TYPEOTHER)

OS_C <- make_multistate(OS_C)

vars_all <- names(OS_C)[grepl("AGE|GENDER|PROFILE_AFTER_TREATMENT|PANEL_VERSION|ANCESTRY_PC_NS|ANCESTRY_PC_ASH",names(OS_C)) & !grepl("\\.[1-3]$",names(OS_C))]
vars_1 <- names(OS_C)[grepl("CANCER_TYPE|START_YEAR|TREATMENT_CUM",names(OS_C)) & grepl("\\.1$",names(OS_C))]
vars_2 <- names(OS_C)[grepl("CANCER_TYPE|START_YEAR|TREATMENT_CUM|LINE_OF_TREATMENT|BIOPSY_SITE_TYPE",names(OS_C)) & grepl("\\.2$",names(OS_C))]
vars_3 <- names(OS_C)[grepl("CANCER_TYPE|START_YEAR|TREATMENT_CUM|LINE_OF_TREATMENT|BIOPSY_SITE_TYPE",names(OS_C)) & grepl("\\.3$",names(OS_C))]

vars <- c(vars_all,vars_1,vars_2,vars_3)

# only keep vars that are finite and pass prop-hazard
prophazard_pval_cutoff <- 0.01
formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars,collapse=" + "),sep="")
fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
notfinite_vars <- rownames(summary(fit)$conf.int)[!is.finite(summary(fit)$conf.int[,"upper .95"])]
vars <- vars[!vars %in% notfinite_vars]
formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars,collapse=" + "),sep="")
fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
strat_vars <- c()

while(length(vars[unlist(lapply(vars,function(x) grepl(x,paste(unlist((setDT(data.frame(suppressWarnings(cox.zph(fit))$table), keep.rownames = TRUE)[] %>% filter(p<prophazard_pval_cutoff))$rn),collapse=""))))] > 0)) {
    df <- data.frame(suppressWarnings(cox.zph(fit))$table)
    strat_vars <- c(strat_vars,vars[unlist(lapply(vars,function(x) grepl(x,paste(unlist((setDT(df, keep.rownames = TRUE)[] %>% filter(p<prophazard_pval_cutoff))$rn),collapse=""))))])
    if (length(strat_vars)==0) {
        break
    }
    vars <- vars[!vars %in% strat_vars]
    cont_strat_vars <- strat_vars[(sapply(OS_C[strat_vars], function(x) is.double(x) && length(unique(x))>2))]
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
        formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars,collapse=" + "),sep="")
    } else {
        formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars,collapse=" + ")," + strata(",paste(strat_vars,collapse=","),")",sep="")
    }
    fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
}

if (length(strat_vars)>0) {
    form <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + (1|MRN) + dosage.1 + dosage.2 + dosage.3 + ",paste(vars,collapse=" + ")," + strata(",paste(strat_vars,collapse=","),")",sep="")
} else {
    form <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + (1|MRN) + dosage.1 + dosage.2 + dosage.3 + ",paste(vars,collapse=" + "),sep="")
}                                                  

result_cancer_chrom <- data.frame()
len_of_loop <- length(rownames(allele_in))
pb <- txtProgressBar(min = 1, max = len_of_loop, style = 3)

for (i in 1:len_of_loop) {
    df <- setDT(data.frame(t(allele_in[i,])),keep.rownames = T)[] %>% filter(rn!="SNP")
    df <- df %>% rename(PBP=rn,dosage=contains("X")) %>% mutate(dosage=as.numeric(as.character(dosage)),PBP=gsub("^0_","",PBP))
    OS_df <- suppressWarnings(inner_join(OS_C,df,by="PBP"))
    OS_df <- OS_df %>% mutate(dosage.1 = ifelse(trans==1,dosage,0), dosage.2 = ifelse(trans==2,dosage,0), dosage.3 = ifelse(trans==3,dosage,0))
    fit <- coxme(as.formula(form),OS_df)
    haz <- fit$coefficients[["dosage.1"]]
    var <- diag(vcov(fit))[[1]]
    p_wald <- 2*pnorm(-abs(haz/sqrt(var)))

    temp <- data.frame(allele=allele_in[i,1],loghazard=haz,var=var,p_wald=p_wald)
    result_cancer_chrom <- rbind(result_cancer_chrom,temp)
    setTxtProgressBar(pb, i)    
}

fwrite(result_cancer_chrom,paste("result_coxme/",chr_number,".part",part,sep=""))
