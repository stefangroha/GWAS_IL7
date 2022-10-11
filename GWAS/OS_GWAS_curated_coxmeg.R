suppressPackageStartupMessages(library(mltools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(coxmeg))
suppressPackageStartupMessages(library(mbend))


# input 
args = commandArgs(trailingOnly=TRUE)
chr_number = args[[1]]
part = args[[2]]


# get sample names
sample_names <- c(t(suppressWarnings(read.table("/data/gusev/PROFILE/GWAS_Surv/IO_COHORT/KINSHIP/kinship.matrix.rel.id"))))
# load dosages
allele_in <- fread(paste("/data/gusev/PROFILE/GWAS_Surv/IO_COHORT/DOSAGES_split/PROF_2020_04.",chr_number,".QC.part",part,".traw",sep=""),data.table=F) %>% mutate(allele=paste(CHR,"_",POS,sep="")) %>% dplyr::select(-CHR,-POS,-SNP,-`(C)M`,-COUNTED,-ALT)
colnames(allele_in) <- gsub("^0_","",colnames(allele_in))
rownames(allele_in) <- allele_in$allele
allele_in <- allele_in %>% dplyr::select(-allele)
X <- t(allele_in)
X <- X[sample_names,]

# OS
pid2bl <- fread("/PHShome/szg99/immuno/Profile/OncDRS_data/ALL/GENOMIC_SPECIMEN.csv",fill=TRUE) %>% select(PATIENT_ID,SAMPLE_ACCESSION_NBR) %>% rename(BLNUM=SAMPLE_ACCESSION_NBR) %>% unique
bl2cbio <- fread("/PHShome/szg99/immuno/Profile/OncDRS_data/ALL/sample_id_to_cbio.txt") %>% select(SampleID,cBio_PatientID,cBio_SampleID) %>% rename(BLNUM=SampleID,CBIO_PATIENT=cBio_PatientID,CBIO_SAMPLE=cBio_SampleID) %>% unique
bl2cbio <- full_join(fread("/PHShome/szg99/immuno/Profile/OncDRS_data/LUNG/raw/PTID_CBIO.csv") %>% select(BLNUM,CBIO_PATIENT,CBIO_SAMPLE),bl2cbio,by=c("BLNUM","CBIO_PATIENT","CBIO_SAMPLE")) %>% distinct
pid2mrn <- fread("/PHShome/szg99/immuno/Profile/OncDRS_data/ALL/PT_INFO_STATUS_REGISTRATION.csv",fill=TRUE) %>% select(PATIENT_ID,DFCI_MRN) %>% rename(MRN=DFCI_MRN) %>% unique
pbp2bl <- fread("/data/gusev/PROFILE/2020_04/gusev-pseudonymize-output.txt") %>% rename(BLNUM=bl_number,PBP=pbp_number) %>% select(BLNUM,PBP)
temp <- right_join(bl2cbio,pid2bl,by="BLNUM") %>% group_by(PATIENT_ID) %>% summarise(CBIO_PATIENT=as.character(ifelse(all(is.na(CBIO_PATIENT)),NA,first(CBIO_PATIENT[!is.na(CBIO_PATIENT)]))))
translate_pid <- merge(temp,pid2mrn,by="PATIENT_ID") %>% select(PATIENT_ID,CBIO_PATIENT,MRN) %>% unique
translate_sid <- full_join(pbp2bl,right_join(bl2cbio,pid2bl,by="BLNUM"),by="BLNUM") %>% select(PATIENT_ID,BLNUM,PBP,CBIO_PATIENT,CBIO_SAMPLE) %>% unique

OS <- fread("../../OS/OS.csv",na.string="")
OS <- OS %>% mutate(LAST_DATE=as.Date(LAST_DATE),IO_START=as.Date(IO_START),IO_STOP=as.Date(IO_STOP))

OS_C <- OS# %>% filter(CANCER_TYPE==gsub("_"," ",cancer))
vars=c("CANCER_TYPE","START_YEAR","GENDER","PROFILE_AFTER_TREATMENT","TREATMENT_CUM","LINE_OF_TREATMENT","PANEL_VERSION","AGE_AT_TREATMENTSTART","ANCESTRY_EA_PC1","ANCESTRY_EA_PC2")
OS_C$ANCESTRY_EA_PC1 <- scale(OS_C$ANCESTRY_EA_PC1)
OS_C$ANCESTRY_EA_PC2 <- scale(OS_C$ANCESTRY_EA_PC2)
pbp2mrn <- inner_join(translate_pid %>% select(PATIENT_ID,MRN),translate_sid %>% select(PATIENT_ID,PBP),by="PATIENT_ID") %>% select(-PATIENT_ID)
OS_C <- inner_join(OS_C,pbp2mrn,by="MRN") %>% select(MRN,PBP,IO_START,tstart,tstop,event,one_of(vars)) %>% drop_na
mm <- model.matrix(as.formula(paste("~ ",paste(vars,collapse=" + "),sep="")), OS_C)  
OS_C <- OS_C %>% dplyr::select(-one_of(vars)) #%>% cbind(x)
OS_C <- cbind(OS_C,mm)
vars <- c(unlist(names(OS_C)[grepl(paste(vars,collapse="|"),names(OS_C))]))

# merge with first variant to remove samples which are NaN
df <- data.frame(PBP=names(allele_in),dosage=as.vector(as.matrix(allele_in)[1,]))
OS_C <- inner_join(OS_C,df) %>% dplyr::select(-dosage)
OS_C <- OS_C[match(sample_names,OS_C$PBP),]
row.names(OS_C) <- NULL
OS_C <- drop_na(OS_C)

temp <- merge(OS_C,fread("../curated_irAE_full.csv",na.strings="") %>% select(MRN,ADV_EVENT,ADV_EVENT_DT),by=c("MRN")) %>% mutate(irAE_time=as.integer(ifelse(ADV_EVENT,as.integer(difftime(ADV_EVENT_DT,IO_START,unit="days")),tstop)),irAE=as.numeric(ifelse(ADV_EVENT,1,0))) %>% select(-ADV_EVENT_DT,-ADV_EVENT)
temp <- temp %>% mutate(irAE_time=as.integer(ifelse(irAE_time==0,1,irAE_time)),tstop=as.integer(ifelse(tstop<=irAE_time,irAE_time,tstop)))
OS_C <- temp %>% dplyr::select(-tstart,-tstop,-event,-IO_START)

outcome <- as.matrix(OS_C[,c("irAE_time","irAE")])
rownames(outcome) <- OS_C$PBP
cov <- as.matrix(OS_C %>% select(-PBP,-MRN,-irAE_time,-irAE,-`(Intercept)`))
colnames(cov) <- gsub("/| |>|<|-","_",colnames(cov))
rownames(cov) <- OS_C$PBP

# genotype matrix
kinship_raw <- as.matrix(read.table("/data/gusev/PROFILE//GWAS_Surv/IO_COHORT/KINSHIP/kinship.matrix.rel"))
colnames(kinship_raw) <- sample_names
rownames(kinship_raw) <- colnames(kinship_raw)
kinship <- bend(kinship_raw)
kinship <- cov2cor(kinship)
kinship <- bend(kinship)

X <- X[intersect(sample_names,OS_C$PBP),]
kinship <- kinship[intersect(sample_names,OS_C$PBP),intersect(sample_names,OS_C$PBP)]
outcome <- outcome[intersect(sample_names,OS_C$PBP),]
cov <- cov[intersect(sample_names,OS_C$PBP),]

# fit procedure
fit <- coxmeg_m(X,outcome,kinship,"dense",cov=cov)
result <- cbind(allele=colnames(X),fit$summary) 

fwrite(result,paste("result_curated_coxmeg/",chr_number,"_part",part,sep=""))
