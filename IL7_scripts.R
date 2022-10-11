qqunif.plot<-function(pvalues, 
	should.thin=T, thin.obs.places=2, thin.exp.places=2, 
	xlab=expression(paste("Expected (",-log[10], " p-value)")),
	ylab=expression(paste("Observed (",-log[10], " p-value)")), 
	draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
	already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
	par.settings=list(superpose.symbol=list(pch=pch)), ...) {
	
	
	#error checking
	if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
	if(!(class(pvalues)=="numeric" || 
		(class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
		stop("pvalue vector is not numeric, can't draw plot")
	if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
	if (already.transformed==FALSE) {
		if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
	} else {
		if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
	}
	
	
	grp<-NULL
	n<-1
	exp.x<-c()
	if(is.list(pvalues)) {
		nn<-sapply(pvalues, length)
		rs<-cumsum(nn)
		re<-rs-nn+1
		n<-min(nn)
		if (!is.null(names(pvalues))) {
			grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
			names(pvalues)<-NULL
		} else {
			grp=factor(rep(1:length(pvalues), nn))
		}
		pvo<-pvalues
		pvalues<-numeric(sum(nn))
		exp.x<-numeric(sum(nn))
		for(i in 1:length(pvo)) {
			if (!already.transformed) {
				pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
				exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
			} else {
				pvalues[rs[i]:re[i]] <- pvo[[i]]
				exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
			}
		}
	} else {
		n <- length(pvalues)+1
		if (!already.transformed) {
			exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
			pvalues <- -log10(pvalues)
		} else {
			exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
		}
	}


	#this is a helper function to draw the confidence interval
	panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
		require(grid)
		conf.points = min(conf.points, n-1);
		mpts<-matrix(nrow=conf.points*2, ncol=2)
        	for(i in seq(from=1, to=conf.points)) {
            		mpts[i,1]<- -log10((i-.5)/n)
            		mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
            		mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
            		mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
        	}
        	grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
    	}

	#reduce number of points to plot
	if (should.thin==T) {
		if (!is.null(grp)) {
			thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
				exp.x = round(exp.x, thin.exp.places),
				grp=grp))
			grp = thin$grp
		} else {
			thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
				exp.x = round(exp.x, thin.exp.places)))
		}
		pvalues <- thin$pvalues
		exp.x <- thin$exp.x
	}
	gc()
	
	prepanel.qqunif= function(x,y,...) {
		A = list()
		A$xlim = range(x, y)*1.02
		A$xlim[1]=0
		A$ylim = A$xlim
		return(A)
	}

	#draw the plot
	xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
		prepanel=prepanel, scales=list(axs="i"), pch=pch,
		panel = function(x, y, ...) {
			if (draw.conf) {
				panel.qqconf(n, conf.points=conf.points, 
					conf.col=conf.col, conf.alpha=conf.alpha)
			};
			panel.xyplot(x,y, ...);
			panel.abline(0,1);
		}, par.settings=par.settings, ...
	)
}

library(gridExtra)
library(ttutils)

prophazard_pval_cutoff <- 0.01
prophazard_pval_cutoff_meta <- 0
vars_vec <- c("AGE_AT_TREATMENTSTART","GENDER","START_YEAR","PROFILE_AFTER_TREATMENT","TREATMENT_CUM","CONC_CHEMO","PANEL_VERSION","LINE_OF_TREATMENT")
vie <- list(trans1=c("TREATMENT_CUM","AGE_AT_TREATMENTSTART","GENDER"),trans2=c("AGE_AT_TREATMENTSTART","TREATMENT_CUM","CONC_CHEMO","LINE_OF_TREATMENT","START_YEAR","GENDER"),trans3=c("AGE_AT_TREATMENTSTART","TREATMENT_CUM","LINE_OF_TREATMENT","Tstart"))
via <- c("PANEL_VERSION","PROFILE_AFTER_TREATMENT","ANCESTRY_PC_ASH","ANCESTRY_PC_NS","ANCESTRY_ROT_PC1","ANCESTRY_ROT_PC2")


find_vars_fit <- function(vars_func,cov_vec,OS_C,strat_cancer=F) {
    nonstrat_vars <- c("ANCESTRY_ROT_PC1","ANCESTRY_ROT_PC2","ANCESTRY_PC_ASH","ANCESTRY_PC_NS","AGE_AT_TREATMENTSTART","TMB")
    nonstrat_vars <- c(nonstrat_vars,paste(nonstrat_vars,c(".1",".2",".3"),sep=""))
    
    for (cv in vars_func) {
        if (grepl("\\.[1-3]$",cv)) {
            tr <- gsub(".*\\.([1-3])$","\\1",cv)
            if (length(unique(OS_C[!is.na(OS_C[[cv]]) & OS_C[["trans"]]==tr,][[cv]]))<2) {
                vars_func <- vars_func[vars_func!=cv]
            }
        } else {
            if (length(unique(OS_C[!is.na(OS_C[[cv]]),][[cv]]))<2) {
               vars_func <- vars_func[vars_func!=cv] 
            }
        }
    } 
    if (strat_cancer) {
        formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans,CANCER_TYPE) + ",paste(vars_func,collapse=" + ")," + ",paste(cov_vec,collapse=" + "),sep="")
    } else {
        formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars_func,collapse=" + ")," + ",paste(cov_vec,collapse=" + "),sep="")
    }
    fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
#     if (strat_cancer) {
#         formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans,CANCER_TYPE) + ",paste(vars_func,collapse=" + ")," + ",paste(cov_vec,collapse=" + "),sep="")
#     } else {
#         formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars_func,collapse=" + ")," + ",paste(cov_vec,collapse=" + "),sep="")
#     }
#     fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
    strat_vars <- c()
    while(length(vars_func[unlist(lapply(vars_func,function(x) grepl(x,paste(unlist((setDT(data.frame(suppressWarnings(cox.zph(fit))$table), keep.rownames = TRUE)[] %>% filter(p<prophazard_pval_cutoff))$rn),collapse=""))))] > 0)) {
        df <- data.frame(suppressWarnings(cox.zph(fit))$table)
        strat_vars <- c(strat_vars,vars_func[unlist(lapply(vars_func,function(x) grepl(x,paste(unlist((setDT(df, keep.rownames = TRUE)[] %>% filter(p<prophazard_pval_cutoff))$rn),collapse=""))))])
        if (length(strat_vars)==0) {
            break
        }
        vars_func <- vars_func[!vars_func %in% strat_vars]  
        cont_strat_vars <- strat_vars[(sapply(OS_C[strat_vars], function(x) all(!isInteger(x))))]
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
            if (strat_cancer) {
                formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans,CANCER_TYPE) + ",paste(vars_func,collapse=" + ")," + ",paste(cov_vec,collapse=" + "),sep="")
            } else {
                formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars_func,collapse=" + ")," + ",paste(cov_vec,collapse=" + "),sep="")
            } 
        } else {
            if (strat_cancer) {
                formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans,CANCER_TYPE) + ",paste(vars_func,collapse=" + ")," + ",paste(cov_vec,collapse=" + ")," + strata(",paste(strat_vars,collapse=","),")",sep="")
            } else {
                formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars_func,collapse=" + ")," + ",paste(cov_vec,collapse=" + ")," + strata(",paste(strat_vars,collapse=","),")",sep="")
            }  
        }
        fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
    }
    vars_func <- c(vars_func,intersect(strat_vars,nonstrat_vars))
    nonfin_df <- data.frame(summary(fit)$conf.int)
    notfinite_vars <- (setDT(nonfin_df,keep.rownames = T)[] %>% filter_all(any_vars(is.double(.) & !is.finite(.))))$rn
    nonfin_df <- data.frame(summary(fit)$coefficients)
    notfinite_vars <- c(notfinite_vars,(setDT(nonfin_df,keep.rownames = T)[] %>% filter_all(any_vars(is.double(.) & !is.finite(.))))$rn)
#     notfinite_vars <- rownames(summary(fit)$conf.int)[!is.finite(summary(fit)$conf.int[,"upper .95"])]
#     notfinite_vars <- c(notfinite_vars,names(exp(fit$coefficients)[!is.finite(exp(fit$coefficients))]))
    vars_func <- vars_func[!vars_func %in% notfinite_vars]
    
    strat_vars <- strat_vars[!strat_vars %in% nonstrat_vars]
    if (length(strat_vars)==0) {
        if (strat_cancer) {
            formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans,CANCER_TYPE) + ",paste(vars_func,collapse=" + ")," + ",paste(cov_vec,collapse=" + "),sep="")
        } else {
            formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars_func,collapse=" + ")," + ",paste(cov_vec,collapse=" + "),sep="")
        } 
    } else {
        if (strat_cancer) {
            formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans,CANCER_TYPE) + ",paste(vars_func,collapse=" + ")," + ",paste(cov_vec,collapse=" + ")," + strata(",paste(strat_vars,collapse=","),")",sep="")
        } else {
            formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars_func,collapse=" + ")," + ",paste(cov_vec,collapse=" + ")," + strata(",paste(strat_vars,collapse=","),")",sep="")
        }  
    }
    fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
    return(fit)
}      


meta_survival_AE <- function(cov_vec_in,OS_temp,vars_in_each=vie,vars_in_all=via,vars_direct=c(),transition="1|2|3",plot=TRUE,name=c(),minx=c(),maxx=c(),cancer_variable="CANCER_TYPE",cancer_leaveout=c()) {
    temp <- OS_temp %>% select(MRN,trans,one_of(vars_in_all)) %>% drop_na
    x <- model.matrix(as.formula(paste("~ MRN + trans + ",paste(vars_in_all,collapse=" + "),sep="")), temp)
    OS_temp <- OS_temp %>% select(-one_of(vars_in_all)) #%>% cbind(x)
    OS_temp <- merge(OS_temp,x,by=c("MRN","trans"))

    vars_trans1 <- vars_in_each$trans1
    vars_trans1 <- names(OS_temp)[grepl(paste(vars_trans1,collapse="|"),names(OS_temp)) & grepl("\\.1$",names(OS_temp))]
    vars_trans2 <- vars_in_each$trans2
    vars_trans2 <- names(OS_temp)[grepl(paste(vars_trans2,collapse="|"),names(OS_temp)) & grepl("\\.2$",names(OS_temp))]
    vars_trans3 <- vars_in_each$trans3
    vars_trans3 <- names(OS_temp)[grepl(paste(vars_trans3,collapse="|"),names(OS_temp)) & grepl("\\.3$",names(OS_temp))]

    cov_vec <- names(OS_temp)[grepl(paste("^",paste(cov_vec_in,collapse="\\.[1-3]$|^"),"\\.[1-3]$",sep=""),names(OS_temp)) & grepl(paste("\\.",transition,"$",sep=""),names(OS_temp))]
    vars_trans1 <- vars_trans1[!grepl(paste(cov_vec,collapse="|"),vars_trans1)]
    vars_trans2 <- vars_trans2[!grepl(paste(cov_vec,collapse="|"),vars_trans2)]
    vars_trans3 <- vars_trans3[!grepl(paste(cov_vec,collapse="|"),vars_trans3)]
    
    vars_in_all <- names(OS_temp)[grepl(paste(vars_in_all,collapse="|"),names(OS_temp))& !grepl("\\.[1-3]$",names(OS_temp))]
    vars_in_all <- vars_in_all[!grepl(paste(cov_vec,collapse="|"),vars_in_all)]
    
    vars <- c(vars_trans1,vars_trans2,vars_trans3,vars_in_all,vars_direct)
    
    OS_temp <- OS_temp %>% select(MRN,matches(cancer_variable),Tstart,Tstop,status,trans,one_of(vars),one_of(cov_vec)) %>% drop_na %>% filter(grepl(transition,as.character(trans)))
    
    df_hazard <- data.frame(cancer=c(),variables=c(),log_hr=c(),log_hr_sd=c(),pval=c(),n=c(),pval_prophaz=c())
    if (length(cancer_leaveout)>0) {
        cancers <- unique(OS_temp[,cancer_variable])[!grepl(paste(cancer_leaveout,collapse="|"),unique(OS_temp[,cancer_variable]))]
    } else {
        cancers <- unique(OS_temp[,cancer_variable])
    }
    for (c in cancers) {
#         print(c)
        tryCatch({
            OS_C <- OS_temp %>% filter(get(cancer_variable)==c)
    #         OS_C <- OS_temp %>% filter(CANCER_TYPE==c) # %>% filter(trans==1)
            if (sum(OS_C$status)>0) {
                # make formula
                fit <- find_vars_fit(vars,cov_vec,OS_C)
                vars_meta <- unlist(lapply(cov_vec,function(x) rownames(summary(fit)$coefficients)[grepl(x,rownames(summary(fit)$coefficients))]))
                temp <- data.table(cancer=c,variables=vars_meta,log_hr=summary(fit)$coefficients[vars_meta,c("coef")],log_hr_sd=summary(fit)$coefficients[vars_meta,c("se(coef)")],pval=summary(fit)$coefficients[vars_meta,c("Pr(>|z|)")])
                temp$n <- dim(OS_C %>% distinct(MRN))[[1]]
                temp$pval_prophaz <- suppressWarnings(cox.zph(fit))$table[vars_meta,"p"]
                df_hazard <- rbind(df_hazard,temp)
            }
        }, 
            error=function (e) {
                print(paste("Error:",e))
                print(paste("Cancer",c,"was unsuccesful"))}
        )
    }
    df_hazard_ret <- df_hazard                              
    df_hazard <- df_hazard %>% drop_na(pval_prophaz) %>% filter(pval_prophaz>prophazard_pval_cutoff_meta) %>% filter(!is.na(log_hr)) %>% filter(!(abs(log_hr_sd/log_hr)>100 & abs(log_hr)>0.01))

    df_hazard$p <- df_hazard$pval
    df_hazard$p_prophaz <- df_hazard$pval_prophaz
#     df_hazard$p <- signif(df_hazard$pval,2)
#     df_hazard$p_prophaz <- signif(df_hazard$pval_prophaz,2)
    res <- df_hazard %>% group_by(variables) %>% do(result=metagen(log_hr,log_hr_sd,data=.,studlab=paste(cancer),comb.random=F,comb.fixed=T,method.tau = "SJ",hakn = TRUE,prediction=F,sm="SMD"))
#     if(plot) {
#         if (length(res$result)>0) {
#             for (i in 1:length(res$result)) {
#      forest(res$result[[i]],leftcols=c("studlab","TE","seTE","n","p","p_prophaz"),print.tau2=F,print.pval.Q=F,print.I2=F,digits=4)
#                     grid.text(res$variables[[i]], 0.5,0.95, gp=gpar(cex=1.5))
#                     grid.text(paste("p=",signif(res$result[[i]]$pval.fixed,2),sep=""), 0.3,0.225, gp=gpar(cex=1))
#                     grid.text(paste("p=",signif(res$result[[i]]$pval.random,2),sep=""), 0.3,0.175, gp=gpar(cex=1))
#             }
#         }
#     }
    if(plot) {
        if (length(res$result)>0) {
            for (i in 1:length(res$result)) {
                res$result[[i]]$n_print <- as.character(as.integer(res$result[[i]]$data$n))
                res$result[[i]]$pval_print <- as.character(formatC(signif(as.double(res$result[[i]]$pval),2),digits=2,format="fg", flag="#"))
                res$result[[i]]$TE_print <- as.character(formatC(signif(as.double(res$result[[i]]$TE),2),digits=2,format="fg", flag="#"))
#                 print(names(res$result[[i]]$df.Q))
                fp <- grid.grabExpr(print(forest(res$result[[i]],colgap=unit(0,"points"),comb.random=F,leftcols=c("studlab","TE_print","n_print","pval_print"),leftlabs=c("Cancer Type","logHR","N","p"),rightcols=c("ci"),ref=0,print.tau2=F,print.pval.Q=F,print.I2=F,smlab="logHR",xlab="logHR",xlim=c(ifelse(length(minx)==0,ifelse(min(res$result[[i]]$TE) < -1,1.2,ifelse(min(res$result[[i]]$TE)>0,0.8,2))*min(res$result[[i]]$TE),minx[[i]]),ifelse(length(maxx)==0,ifelse(max(res$result[[i]]$TE)>1,1.2,ifelse(max(res$result[[i]]$TE)<0,0.8,2))*max(res$result[[i]]$TE),maxx[[i]])),text.fixed=paste("Fixed effects model:    logHR=",signif(res$result[[i]]$TE.fixed,2),", p=",signif(res$result[[i]]$pval.fixed,2),sep=""))))
#                 fp <- grid.grabExpr(print(forest(res$result[[i]],fontsize=24,plotwidth=unit(250, "points"),squaresize=0.8,colgap=unit(0,"points"),comb.random=F,leftcols=c("studlab","TE_print","n_print","pval_print"),leftlabs=c("Cancer Type","logHR","N","p"),rightcols=c("ci"),ref=0,print.tau2=F,print.pval.Q=F,print.I2=F,smlab="logHR",xlab="logHR",xlim=c(ifelse(length(min)==0,ifelse(min(res$result[[i]]$TE) < -1,1.2,ifelse(min(res$result[[i]]$TE)>0,0.8,2))*min(res$result[[i]]$TE),min[[i]]),ifelse(length(max)==0,ifelse(max(res$result[[i]]$TE)>1,1.2,ifelse(max(res$result[[i]]$TE)<0,0.8,2))*max(res$result[[i]]$TE),max[[i]])),text.fixed=paste("Fixed effects model:    logHR=",signif(res$result[[i]]$TE.fixed,2),", p=",signif(res$result[[i]]$pval.fixed,2),sep=""))))
                title <- textGrob(ifelse(length(name)==0,res$variables[[i]],name[[i]]), gp=gpar(fontface="bold"))
#                 png(paste(res$variables[[i]],".png",sep=""),width=1124,height=560,pointsize = 28)
                grid.arrange(top=title, fp)
#                 dev.off()
            }
        }
    }
    return(list(res,df_hazard_ret))
}   
                                       
selectcancer_survival_AE <- function(cancer,cov_vec_in,OS_temp,vars_in_each=vie, vars_in_all=via,vars_direct=c(),transition="1|2|3") {
    temp <- OS_temp %>% select(MRN,trans,one_of(vars_in_all)) %>% drop_na
    x <- model.matrix(as.formula(paste("~ MRN + trans + ",paste(vars_in_all,collapse=" + "),sep="")), temp)
    OS_temp <- OS_temp %>% select(-one_of(vars_in_all)) #%>% cbind(x)
    OS_temp <- merge(OS_temp,x,by=c("MRN","trans"))

    vars_trans1 <- vars_in_each$trans1
    vars_trans1 <- names(OS_temp)[grepl(paste(vars_trans1,collapse="|"),names(OS_temp)) & grepl("\\.1$",names(OS_temp))]
    vars_trans2 <- vars_in_each$trans2
    vars_trans2 <- names(OS_temp)[grepl(paste(vars_trans2,collapse="|"),names(OS_temp)) & grepl("\\.2$",names(OS_temp))]
    vars_trans3 <- vars_in_each$trans3
    vars_trans3 <- names(OS_temp)[grepl(paste(vars_trans3,collapse="|"),names(OS_temp)) & grepl("\\.3$",names(OS_temp))]

    cov_vec <- names(OS_temp)[grepl(paste("^",paste(cov_vec_in,collapse="\\.[1-3]$|^"),"\\.[1-3]$",sep=""),names(OS_temp)) & grepl(paste("\\.",transition,"$",sep=""),names(OS_temp))]
    vars_trans1 <- vars_trans1[!grepl(paste(cov_vec,collapse="|"),vars_trans1)]
    vars_trans2 <- vars_trans2[!grepl(paste(cov_vec,collapse="|"),vars_trans2)]
    vars_trans3 <- vars_trans3[!grepl(paste(cov_vec,collapse="|"),vars_trans3)]
    
    vars_in_all <- names(OS_temp)[grepl(paste(vars_in_all,collapse="|"),names(OS_temp))& !grepl("\\.[1-3]$",names(OS_temp))]
    vars_in_all <- vars_in_all[!grepl(paste(cov_vec,collapse="|"),vars_in_all)]
    
    vars <- c(vars_trans1,vars_trans2,vars_trans3,vars_in_all,vars_direct)
    
    OS_temp <- OS_temp %>% select(CANCER_TYPE,Tstart,Tstop,status,trans,one_of(vars),one_of(cov_vec)) %>% drop_na %>% filter(grepl(transition,as.character(trans)))
    
    OS_C <- OS_temp %>% filter(CANCER_TYPE==cancer)
    if (sum(OS_C$status)>0) {
        # make formula
        fit <- find_vars_fit(vars,cov_vec,OS_C)
    } else {
        fit <- NA
    }
    return(fit)
}


meta_survival_onehot <- function(cov_vec,OS_temp,vars_in=vars_vec,plot=TRUE,name=c(),minx=c(),maxx=c()) {
    df_hazard <- data.frame(cancer=c(),variables=c(),log_hr=c(),log_hr_sd=c(),pval=c(),n=c(),pval_prophaz=c())
    for (c in unique(OS_temp$CANCER_TYPE)) {
        tryCatch({
            OS_C <- OS_temp %>% filter(CANCER_TYPE==c) %>% select(MRN,tstart,tstop,event,one_of(vars_in),one_of(cov_vec)) %>% drop_na
            if (length(levels(factor(OS_C[,cov_vec])))>1) {
                if (sum(OS_C$event)>0) {
                    vars <- vars_in
                    vars_cont <- vars[vars %in% names(OS_C)[unlist(sapply(as_tibble(OS_C), is.double))]]
                    vars <- vars[!vars %in% vars_cont]
                    x <- apply(model.matrix(as.formula(paste("~ ",paste(vars,collapse=" + "),sep="")), OS_C),2,as.integer)         
                    OS_C <- OS_C %>% select(-one_of(vars))
                    OS_C <- cbind(OS_C,x)
                    vars <- c(unlist(names(OS_C)[grepl(paste(vars,collapse="|"),names(OS_C))]),vars_cont)
                    formula <- paste("Surv(tstart,tstop,event) ~ ",paste(vars,collapse=" + ")," + ",paste(cov_vec,collapse=" + "),sep="")
                    fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
                    nonfinite_vars <- names(fit$coefficients)[unlist(is.na(fit$coefficients))]
                    vars <- vars[!vars %in% nonfinite_vars]
                    formula <- paste("Surv(tstart,tstop,event) ~ ",paste(vars,collapse=" + ")," + ",paste(cov_vec,collapse=" + "),sep="")
                    fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
                    strat_vars <- c()
                    while(length(vars[unlist(lapply(vars,function(x) grepl(x,paste(unlist((setDT(data.frame(suppressWarnings(cox.zph(fit))$table), keep.rownames = TRUE)[] %>% filter(p<prophazard_pval_cutoff))$rn),collapse=""))))] > 0)) {
                        df <- data.frame(suppressWarnings(cox.zph(fit))$table)
                        strat_vars <- c(strat_vars,vars[unlist(lapply(vars,function(x) grepl(x,paste(unlist((setDT(df, keep.rownames = TRUE)[] %>% filter(p<prophazard_pval_cutoff))$rn),collapse=""))))])
                        if (length(strat_vars)==0) {
                            break
                        }
                        vars <- vars[!vars %in% strat_vars]
                        cont_strat_vars <- strat_vars[(sapply(OS_C[strat_vars], function(x) all(!isInteger(x))))]
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
                            formula <- paste("Surv(tstart,tstop,event) ~ ",paste(vars,collapse=" + ")," + ",paste(cov_vec,collapse=" + "),sep="")
                        } else {
                            formula <- paste("Surv(tstart,tstop,event) ~ ",paste(vars,collapse=" + ")," + ",paste(cov_vec,collapse=" + ")," + strata(",paste(strat_vars,collapse=","),")",sep="")
                        }                                         
                        fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
                    }   
    #                 print(formula)
    #                 print(fit)
    #                 print(cox.zph(fit))
                    vars <- unlist(lapply(cov_vec,function(x) rownames(summary(fit)$coefficients)[grepl(x,rownames(summary(fit)$coefficients))]))
                    temp <- data.table(cancer=c,variables=vars,log_hr=summary(fit)$coefficients[vars,c("coef")],log_hr_sd=summary(fit)$coefficients[vars,c("se(coef)")],pval=summary(fit)$coefficients[vars,c("Pr(>|z|)")])
                    temp$n <- dim(OS_C %>% distinct(MRN))[[1]]
                    temp$pval_prophaz <- suppressWarnings(cox.zph(fit))$table[vars,"p"]
                    df_hazard <- rbind(df_hazard,temp)
                }
            }
        }, 
            error=function (e) {
                print(paste("Error:",e))
                print(paste("Cancer",c,"was unsuccesful"))}
        )
    }
    if(all(dim(df_hazard)==c(0,0))) {
        return(list(NA,NA))
    }
    df_hazard_ret <- df_hazard                             
    df_hazard <- df_hazard %>% drop_na(pval_prophaz) %>% filter(pval_prophaz>prophazard_pval_cutoff_meta) %>% filter(!is.na(log_hr)) %>% filter(!(abs(log_hr_sd/log_hr)>100 & abs(log_hr)>0.01))
#     df_hazard$p <- signif(df_hazard$pval,2)
#     df_hazard$p_prophaz <- signif(df_hazard$pval_prophaz,2)
    res <- df_hazard %>% group_by(variables) %>% do(result=suppressWarnings(metagen(log_hr,log_hr_sd,data=.,studlab=paste(cancer),comb.random=T,comb.fixed=T,method.tau = "SJ",hakn = TRUE,prediction=F,sm="SMD")))
    if(plot) {
        if (length(res$result)>0) {
            for (i in 1:length(res$result)) {
                res$result[[i]]$n_print <- as.character(as.integer(res$result[[i]]$data$n))
                res$result[[i]]$pval_print <- as.character(formatC(signif(as.double(res$result[[i]]$pval),2),digits=2,format="fg", flag="#"))
                res$result[[i]]$TE_print <- as.character(formatC(signif(as.double(res$result[[i]]$TE),2),digits=2,format="fg", flag="#"))
#                 print(names(res$result[[i]]$df.Q))
                fp <- grid.grabExpr(print(forest(res$result[[i]],colgap=unit(0,"points"),comb.random=F,leftcols=c("studlab","TE_print","n_print","pval_print"),leftlabs=c("Cancer Type","logHR","N","p"),rightcols=c("ci"),ref=0,print.tau2=F,print.pval.Q=F,print.I2=F,smlab="logHR",xlab="logHR",xlim=c(ifelse(length(minx)==0,ifelse(min(res$result[[i]]$TE) < -1,1.2,ifelse(min(res$result[[i]]$TE)>0,0.8,2))*min(res$result[[i]]$TE),minx[[i]]),ifelse(length(maxx)==0,ifelse(max(res$result[[i]]$TE)>1,1.2,ifelse(max(res$result[[i]]$TE)<0,0.8,2))*max(res$result[[i]]$TE),maxx[[i]])),text.fixed=paste("Fixed effects model:    logHR=",signif(res$result[[i]]$TE.fixed,2),", p=",signif(res$result[[i]]$pval.fixed,2),sep=""))))
#                 fp <- grid.grabExpr(print(forest(res$result[[i]],fontsize=24,plotwidth=unit(250, "points"),squaresize=0.8,colgap=unit(0,"points"),comb.random=F,leftcols=c("studlab","TE_print","n_print","pval_print"),leftlabs=c("Cancer Type","logHR","N","p"),rightcols=c("ci"),ref=0,print.tau2=F,print.pval.Q=F,print.I2=F,smlab="logHR",xlab="logHR",xlim=c(ifelse(length(min)==0,ifelse(min(res$result[[i]]$TE) < -1,1.2,ifelse(min(res$result[[i]]$TE)>0,0.8,2))*min(res$result[[i]]$TE),min[[i]]),ifelse(length(max)==0,ifelse(max(res$result[[i]]$TE)>1,1.2,ifelse(max(res$result[[i]]$TE)<0,0.8,2))*max(res$result[[i]]$TE),max[[i]])),text.fixed=paste("Fixed effects model:    logHR=",signif(res$result[[i]]$TE.fixed,2),", p=",signif(res$result[[i]]$pval.fixed,2),sep=""))))
                title <- textGrob(ifelse(length(name)==0,res$variables[[i]],name[[i]]), gp=gpar(fontface="bold"))
#                 png(paste(res$variables[[i]],".png",sep=""),width=1124,height=560,pointsize = 28)
                grid.arrange(top=title, fp)
#                 dev.off()
            }
        }
    }
    return(list(res,df_hazard_ret))
}
                                   
                                   
selectcancer_survival_onehot <- function(cancer,cov_vec,OS_temp,vars_in=vars_vec) {
    OS_C <- OS_temp %>% filter(CANCER_TYPE==cancer)
    OS_C <- OS_C %>% select(tstart,tstop,event,one_of(vars_in),one_of(cov_vec)) %>% drop_na
    
    vars <- vars_in
    vars_cont <- vars[vars %in% names(OS_C)[unlist(sapply(as_tibble(OS_C), is.double))]]
    vars <- vars[!vars %in% vars_cont]
    x <- apply(model.matrix(as.formula(paste("~ ",paste(vars,collapse=" + "),sep="")), OS_C),2,as.integer)         
    OS_C <- OS_C %>% select(-one_of(vars))
    OS_C <- cbind(OS_C,x)
    vars <- c(unlist(names(OS_C)[grepl(paste(vars,collapse="|"),names(OS_C))]),vars_cont)
    formula <- paste("Surv(tstart,tstop,event) ~ ",paste(vars,collapse=" + ")," + ",paste(cov_vec,collapse=" + "),sep="")
    fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
    nonfinite_vars <- names(fit$coefficients)[unlist(is.na(fit$coefficients))]
    vars <- vars[!vars %in% nonfinite_vars]
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
            formula <- paste("Surv(tstart,tstop,event) ~ ",paste(vars,collapse=" + ")," + ",paste(cov_vec,collapse=" + "),sep="")
        } else {
            formula <- paste("Surv(tstart,tstop,event) ~ ",paste(vars,collapse=" + ")," + ",paste(cov_vec,collapse=" + ")," + strata(",paste(strat_vars,collapse=","),")",sep="")
        }                                         
        fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
    }   
    return(fit)
}    
  
                                          
meta_survival_interaction <- function(cov_vec,inter_var,OS_temp,vars_in=vars_vec,cancer_exclude=c(),plot=TRUE) {
    vars_in <- vars_in[!vars_in %in% cov_vec]
    df_hazard <- data.frame(cancer=c(),variables=c(),log_hr=c(),log_hr_sd=c(),pval=c(),n=c(),pval_prophaz=c())
    for (c in unique(OS_temp$CANCER_TYPE)[!unique(OS_temp$CANCER_TYPE) %in% cancer_exclude]) {
        tryCatch({
            OS_C <- OS_temp %>% filter(CANCER_TYPE==c) %>% select(tstart,tstop,event,one_of(vars_in),one_of(cov_vec)) %>% drop_na
            if (length(levels(factor(OS_C[,cov_vec])))>1) {
                if (sum(OS_C$event)>0) {
                    vars <- vars_in
                    vars_cont <- vars[vars %in% names(OS_C)[unlist(sapply(as_tibble(OS_C), is.double))]]
#                     vars_cont <- vars[vars %in% names(OS_C)[unlist(sapply(as_tibble(OS_C), function(x) all(!isInteger(x))))]]
                    vars <- vars[!vars %in% vars_cont]
                    x <- apply(model.matrix(as.formula(paste("~ ",paste(vars,collapse=" + "),sep="")), OS_C),2,as.integer)         
                    OS_C <- OS_C %>% select(-one_of(vars))
                    OS_C <- cbind(OS_C,x)
                    vars <- c(unlist(names(OS_C)[grepl(paste(vars,collapse="|"),names(OS_C))]),vars_cont)
                    formula <- paste("Surv(tstart,tstop,event) ~ ",paste(vars,collapse=" + ")," + ",paste(cov_vec,collapse="*"),sep="")
                    fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
                    nonfinite_vars <- names(fit$coefficients)[unlist(is.na(fit$coefficients))]
                    vars <- vars[!vars %in% nonfinite_vars]
                    formula <- paste("Surv(tstart,tstop,event) ~ ",paste(vars,collapse=" + ")," + ",paste(cov_vec,collapse="*"),sep="")
                    fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
                    strat_vars <- c()
                    while(length(vars[unlist(lapply(vars,function(x) grepl(x,paste(unlist((setDT(data.frame(suppressWarnings(cox.zph(fit))$table), keep.rownames = TRUE)[] %>% filter(p<prophazard_pval_cutoff))$rn),collapse=""))))] > 0)) {
                        df <- data.frame(suppressWarnings(cox.zph(fit))$table)
                        strat_vars <- c(strat_vars,vars[unlist(lapply(vars,function(x) grepl(x,paste(unlist((setDT(df, keep.rownames = TRUE)[] %>% filter(p<prophazard_pval_cutoff))$rn),collapse=""))))])
                        if (length(strat_vars)==0) {
                            break
                        }
                        vars <- vars[!vars %in% strat_vars]
                        cont_strat_vars <- strat_vars[(sapply(OS_C[strat_vars], function(x) all(!isInteger(x))))]
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
                            formula <- paste("Surv(tstart,tstop,event) ~ ",paste(vars,collapse=" + ")," + ",paste(cov_vec,collapse="*"),sep="")
                        } else {
                            formula <- paste("Surv(tstart,tstop,event) ~ ",paste(vars,collapse=" + ")," + ",paste(cov_vec,collapse="*")," + strata(",paste(strat_vars,collapse=","),")",sep="")
                        }                               
                        fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
                    } 
#                     print(c)
#                     print(formula)
#                     print(fit)
#                     print(cox.zph(fit))
                    vars <- unique(unlist(lapply(cov_vec,function(x) rownames(summary(fit)$coefficients)[grepl(x,rownames(summary(fit)$coefficients))])))
                    vars_interaction <- vars[grepl(":",vars)]
                    for (var in vars_interaction) {
                        vars <- unlist(strsplit(var,":"))
                        var_nottreatment <- vars[!grepl(inter_var,vars)]
                        var_treatment <- vars[grepl(inter_var,vars)]
                        beta <- fit$coefficients[c(var_nottreatment,var_treatment,var)]
                        beta <- c(beta,list(beta_diff1=sum(fit$coefficients[c(var_treatment,var)]),beta_diff2=sum(fit$coefficients[c(var_nottreatment,var)])))
                        names(beta) <- c(paste("beta_1_",var_nottreatment,sep=""),paste("beta_2_",var_nottreatment,sep=""),paste("beta_3_",var_nottreatment,sep=""),paste("beta_3T_",var_nottreatment,sep=""),paste("beta_3m_",var_nottreatment,sep=""))
                        variance_nottreatment <- fit$var[names(fit$coefficients)==var_nottreatment,names(fit$coefficients)==var_nottreatment]
                        variance_treatment <- fit$var[names(fit$coefficients)==var_treatment,names(fit$coefficients)==var_treatment]
                        variance_inter <- fit$var[names(fit$coefficients)==var,names(fit$coefficients)==var]
                        variance_diff1 <- variance_treatment + variance_inter + 2*fit$var[names(fit$coefficients)==var_treatment,names(fit$coefficients)==var]
                        variance_diff2 <- variance_nottreatment + variance_inter + 2*fit$var[names(fit$coefficients)==var_nottreatment,names(fit$coefficients)==var]
                        variance <- list(variance_nottreatment,variance_treatment,variance_inter,variance_diff1,variance_diff2)
                        names(variance) <- gsub("beta_","var_",names(beta))
                        pval <- 2*pnorm(-abs(unlist(beta)/sqrt(unlist(variance))))
                        names(pval) <- gsub("beta_","pval_",names(beta))
                        l <- cbind(beta,variance,pval)
                        temp <- data.frame(matrix(unlist(t(l)), nrow=dim(l)[1], byrow=T))
                        names(temp) <- c("beta","variance","pval")
                        temp[,"variables"] <- gsub("beta_","",names(beta))
                        temp$cancer = c
                        temp$n <- summary(fit)$n
                        df_hazard <- rbind(df_hazard,temp)
                    }
                }
            }
        }, 
            error=function (e) {
                print(paste("Error:",e))
                print(paste("Cancer",c,"was unsuccesful"))}
        )
    }
    if(all(dim(df_hazard)==c(0,0))) {
        return(list(NA,NA))
    }
    df_hazard_ret <- df_hazard
    df_hazard$sd <- sqrt(df_hazard$variance)
    df_hazard <- df_hazard %>% filter(!is.na(beta)) %>% filter(!(abs(sd/beta)>100 & abs(beta)>0.01))
    df_hazard$p <- signif(df_hazard$pval,2)
#     df_hazard$p_prophaz <- signif(df_hazard$pval_prophaz,2)
#     df_hazard <- df_hazard %>% rowwise %>% mutate(variables = paste(sort(unlist(str_split(variables,":"))),collapse=":")) %>% ungroup
    res <- df_hazard %>% group_by(variables) %>% do(result=suppressWarnings(metagen(beta,sd,data=.,studlab=paste(cancer),comb.random=T,comb.fixed=T,method.tau = "SJ",hakn = TRUE,prediction=F,sm="SMD")))
    if(plot) {
        if (length(res$result)>0) {
            for (i in 1:length(res$result)) {
              forest(res$result[[i]],leftcols=c("studlab","TE","seTE","n","p"),print.tau2=F,print.pval.Q=F,print.I2=F)
                grid.text(res$variables[[i]], 0.5,0.95, gp=gpar(cex=1.5))
                grid.text(paste("p=",signif(res$result[[i]]$pval.fixed,2),sep=""), 0.3,0.225, gp=gpar(cex=1))
                grid.text(paste("p=",signif(res$result[[i]]$pval.random,2),sep=""), 0.3,0.175, gp=gpar(cex=1))
            }
        }
    }
    return(list(res,df_hazard_ret))
}                                          
                                          
last_before_treatment <- function(vars,dates,io_start) {
    dates_before <- dates[as.integer(difftime(dates,io_start,unit="days"))< 7]#-14]
#     dates_after <- dates[as.integer(difftime(dates,io_start,unit="days"))>=0]
    if (length(dates_before)>0) {
        return(first(vars[which(dates==max(dates_before)),]))
#     } else if (length(dates_after)>0) {
#         return(first(vars[which(dates==min(dates_after)),]))
    } else {
        return(NA)
    }
}

lastbeforeorfirstafter_treatment  <- function(vars,dates,io_start) {
    dates_before <- dates[as.integer(difftime(dates,io_start,unit="days"))<0]#-14]
    dates_after <- dates[as.integer(difftime(dates,io_start,unit="days"))>=0]
    if (length(dates_before)>0) {
        return(first(vars[which(dates==max(dates_before)),]))
    } else if (length(dates_after)>0) {
        return(first(vars[which(dates==min(dates_after)),]))
    } else {
        return(NA)
    }
}

find_vars_fit_int <- function(vars_func,cov_vec,OS_C,strat_cancer=F) {
    nonstrat_vars <- c("ANCESTRY_ROT_PC1","ANCESTRY_ROT_PC2","ANCESTRY_PC_ASH","ANCESTRY_PC_NS","AGE_AT_TREATMENTSTART","TMB")
    nonstrat_vars <- c(nonstrat_vars,paste(nonstrat_vars,c(".1",".2",".3"),sep=""))
    
    for (cv in vars_func) {
        if (grepl("\\.[1-3]$",cv)) {
            tr <- gsub(".*\\.([1-3])$","\\1",cv)
            if (length(unique(OS_C[!is.na(OS_C[[cv]]) & OS_C[["trans"]]==tr,][[cv]]))<2) {
                vars_func <- vars_func[vars_func!=cv]
            }
        } else {
            if (length(unique(OS_C[!is.na(OS_C[[cv]]),][[cv]]))<2) {
               vars_func <- vars_func[vars_func!=cv] 
            }
        }
    }
    if (strat_cancer) {
        formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans,CANCER_TYPE) + ",paste(vars_func,collapse=" + ")," + ",paste(paste(cov_vec[grepl("\\.1$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.2$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.3$",cov_vec)],collapse="*"),sep=" + "),sep="")
    } else {
        formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars_func,collapse=" + ")," + ",paste(paste(cov_vec[grepl("\\.1$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.2$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.3$",cov_vec)],collapse="*"),sep=" + "),sep="")
    }
    fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
#     notfinite_vars <- rownames(summary(fit)$conf.int)[!is.finite(summary(fit)$conf.int[,"upper .95"])]
#     notfinite_vars <- c(notfinite_vars,names(exp(fit$coefficients)[!is.finite(exp(fit$coefficients))]))
#     vars_func <- vars_func[!vars_func %in% notfinite_vars]
#     if (strat_cancer) {
#         formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans,CANCER_TYPE) + ",paste(vars_func,collapse=" + ")," + ",paste(paste(cov_vec[grepl("\\.1$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.2$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.3$",cov_vec)],collapse="*"),sep=" + "),sep="")
#     } else {
#         formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars_func,collapse=" + ")," + ",paste(paste(cov_vec[grepl("\\.1$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.2$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.3$",cov_vec)],collapse="*"),sep=" + "),sep="")
#     }
#     fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
    strat_vars <- c()
    while(length(vars_func[unlist(lapply(vars_func,function(x) grepl(x,paste(unlist((setDT(data.frame(suppressWarnings(cox.zph(fit))$table), keep.rownames = TRUE)[] %>% filter(p<prophazard_pval_cutoff))$rn),collapse=""))))] > 0)) {
        df <- data.frame(suppressWarnings(cox.zph(fit))$table)
        strat_vars <- c(strat_vars,vars_func[unlist(lapply(vars_func,function(x) grepl(x,paste(unlist((setDT(df, keep.rownames = TRUE)[] %>% filter(p<prophazard_pval_cutoff))$rn),collapse=""))))])
        if (length(strat_vars)==0) {
            break
        }
        vars_func <- vars_func[!vars_func %in% strat_vars]
                                                           
        cont_strat_vars <- strat_vars[(sapply(OS_C[strat_vars], function(x) all(!isInteger(x))))]
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
            if (strat_cancer) {
                formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans,CANCER_TYPE) + ",paste(vars_func,collapse=" + ")," + ",paste(paste(cov_vec[grepl("\\.1$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.2$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.3$",cov_vec)],collapse="*"),sep=" + "),sep="")
            } else {
                formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars_func,collapse=" + ")," + ",paste(paste(cov_vec[grepl("\\.1$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.2$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.3$",cov_vec)],collapse="*"),sep=" + "),sep="")
            } 
        } else {
            if (strat_cancer) {
                formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans,CANCER_TYPE) + ",paste(vars_func,collapse=" + ")," + ",paste(paste(cov_vec[grepl("\\.1$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.2$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.3$",cov_vec)],collapse="*"),sep=" + ")," + strata(",paste(strat_vars,collapse=","),")",sep="")
            } else {
                formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars_func,collapse=" + ")," + ",paste(paste(cov_vec[grepl("\\.1$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.2$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.3$",cov_vec)],collapse="*"),sep=" + ")," + strata(",paste(strat_vars,collapse=","),")",sep="")
            }  
        }
        fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
    }
    vars_func <- c(vars_func,intersect(strat_vars,nonstrat_vars))
    nonfin_df <- data.frame(summary(fit)$conf.int)
    notfinite_vars <- (setDT(nonfin_df,keep.rownames = T)[] %>% filter_all(any_vars(is.double(.) & !is.finite(.))))$rn
    nonfin_df <- data.frame(summary(fit)$coefficients)
    notfinite_vars <- c(notfinite_vars,(setDT(nonfin_df,keep.rownames = T)[] %>% filter_all(any_vars(is.double(.) & !is.finite(.))))$rn)
#     notfinite_vars <- rownames(summary(fit)$conf.int)[!is.finite(summary(fit)$conf.int[,"upper .95"])]
#     notfinite_vars <- c(notfinite_vars,names(exp(fit$coefficients)[!is.finite(exp(fit$coefficients))]))
    vars_func <- vars_func[!vars_func %in% notfinite_vars]                                                       
                                                           
    strat_vars <- strat_vars[!strat_vars %in% nonstrat_vars]
    if (length(strat_vars)==0) {
        if (strat_cancer) {
            formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans,CANCER_TYPE) + ",paste(vars_func,collapse=" + ")," + ",paste(paste(cov_vec[grepl("\\.1$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.2$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.3$",cov_vec)],collapse="*"),sep=" + "),sep="")
        } else {
            formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars_func,collapse=" + ")," + ",paste(paste(cov_vec[grepl("\\.1$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.2$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.3$",cov_vec)],collapse="*"),sep=" + "),sep="")
        } 
    } else {
        if (strat_cancer) {
            formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans,CANCER_TYPE) + ",paste(vars_func,collapse=" + ")," + ",paste(paste(cov_vec[grepl("\\.1$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.2$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.3$",cov_vec)],collapse="*"),sep=" + ")," + strata(",paste(strat_vars,collapse=","),")",sep="")
        } else {
            formula <- paste("Surv(Tstart,Tstop,status) ~ strata(trans) + ",paste(vars_func,collapse=" + ")," + ",paste(paste(cov_vec[grepl("\\.1$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.2$",cov_vec)],collapse="*"),paste(cov_vec[grepl("\\.3$",cov_vec)],collapse="*"),sep=" + ")," + strata(",paste(strat_vars,collapse=","),")",sep="")
        }  
    }
    fit <- suppressWarnings(coxph(as.formula(formula),data=OS_C))
    return(fit)
} 

meta_survival_interaction_AE <- function(cov_vec_in,OS_temp,interaction_var="rs7011565_norm",transition="1|2|3",vars_in_each=vie,vars_in_all=via,vars_direct=c(),cancer_variable="CANCER_TYPE",cancer_leaveout=c(),plot=TRUE) {
    temp <- OS_temp %>% select(MRN,trans,one_of(vars_in_all)) %>% drop_na
    x <- model.matrix(as.formula(paste("~ MRN + trans + ",paste(vars_in_all,collapse=" + "),sep="")), temp)
    OS_temp <- OS_temp %>% select(-one_of(vars_in_all)) #%>% cbind(x)
    OS_temp <- merge(OS_temp,x,by=c("MRN","trans"))

    vars_trans1 <- vars_in_each$trans1
    vars_trans1 <- names(OS_temp)[grepl(paste(vars_trans1,collapse="|"),names(OS_temp)) & grepl("\\.1$",names(OS_temp))]
    vars_trans2 <- vars_in_each$trans2
    vars_trans2 <- names(OS_temp)[grepl(paste(vars_trans2,collapse="|"),names(OS_temp)) & grepl("\\.2$",names(OS_temp))]
    vars_trans3 <- vars_in_each$trans3
    vars_trans3 <- names(OS_temp)[grepl(paste(vars_trans3,collapse="|"),names(OS_temp)) & grepl("\\.3$",names(OS_temp))]

    cov_vec <- names(OS_temp)[grepl(paste("^",paste(cov_vec_in,collapse="\\.[1-3]$|^"),"\\.[1-3]$",sep=""),names(OS_temp)) & grepl(paste("\\.",transition,"$",sep=""),names(OS_temp))]
    vars_trans1 <- vars_trans1[!grepl(paste(cov_vec,collapse="|"),vars_trans1)]
    vars_trans2 <- vars_trans2[!grepl(paste(cov_vec,collapse="|"),vars_trans2)]
    vars_trans3 <- vars_trans3[!grepl(paste(cov_vec,collapse="|"),vars_trans3)]
    
    vars_in_all <- names(OS_temp)[grepl(paste(vars_in_all,collapse="|"),names(OS_temp))& !grepl("\\.[1-3]$",names(OS_temp))]
    vars_in_all <- vars_in_all[!grepl(paste(cov_vec,collapse="|"),vars_in_all)]
    
    vars_include <- c(vars_trans1,vars_trans2,vars_trans3,vars_in_all,vars_direct)
    
    OS_temp <- OS_temp %>% select(matches(cancer_variable),Tstart,Tstop,status,trans,one_of(vars_include),one_of(cov_vec)) %>% drop_na %>% filter(grepl(transition,as.character(trans)))
    
    df_hazard <- data.frame(cancer=c(),variables=c(),log_hr=c(),log_hr_sd=c(),pval=c(),n=c(),pval_prophaz=c())
    if (length(cancer_leaveout)>0) {
        cancers <- unique(OS_temp[,cancer_variable])[!grepl(paste(cancer_leaveout,collapse="|"),unique(OS_temp[,cancer_variable]))]
    } else {
        cancers <- unique(OS_temp[,cancer_variable])
    }
    for (c in cancers) {
        tryCatch({
            OS_C <- OS_temp %>% filter(get(cancer_variable)==c)
    #         OS_C <- OS_temp %>% filter(CANCER_TYPE==c) # %>% filter(trans==1)
            if (sum(OS_C$status)>0) {
                # make formula
                fit <- find_vars_fit_int(vars_include,cov_vec,OS_C)
                                           
                vars <- unique(unlist(lapply(cov_vec,function(x) rownames(summary(fit)$coefficients)[grepl(x,rownames(summary(fit)$coefficients))])))
                vars_interaction <- vars[grepl(":",vars)]
                for (var in vars_interaction) {
                    vars <- unlist(strsplit(var,":"))
                    var_nottreatment <- vars[!grepl(interaction_var,vars)]
                    var_treatment <- vars[grepl(interaction_var,vars)]
                    beta <- fit$coefficients[c(var_nottreatment,var_treatment,var)]
                    beta <- c(beta,list(beta_diff1=sum(fit$coefficients[c(var_treatment,var)]),beta_diff2=sum(fit$coefficients[c(var_nottreatment,var)])))
                    names(beta) <- c(paste("beta_1_",var_nottreatment,sep=""),paste("beta_2_",var_nottreatment,sep=""),paste("beta_3_",var_nottreatment,sep=""),paste("beta_3T_",var_nottreatment,sep=""),paste("beta_3m_",var_nottreatment,sep=""))
                    variance_nottreatment <- fit$var[names(fit$coefficients)==var_nottreatment,names(fit$coefficients)==var_nottreatment]
                    variance_treatment <- fit$var[names(fit$coefficients)==var_treatment,names(fit$coefficients)==var_treatment]
                    variance_inter <- fit$var[names(fit$coefficients)==var,names(fit$coefficients)==var]
                    variance_diff1 <- variance_treatment + variance_inter + 2*fit$var[names(fit$coefficients)==var_treatment,names(fit$coefficients)==var]
                    variance_diff2 <- variance_nottreatment + variance_inter + 2*fit$var[names(fit$coefficients)==var_nottreatment,names(fit$coefficients)==var]
                    variance <- list(variance_nottreatment,variance_treatment,variance_inter,variance_diff1,variance_diff2)
                    names(variance) <- gsub("beta_","var_",names(beta))
                    pval <- 2*pnorm(-abs(unlist(beta)/sqrt(unlist(variance))))
                    names(pval) <- gsub("beta_","pval_",names(beta))
                    l <- cbind(beta,variance,pval)
                    temp <- data.frame(matrix(unlist(t(l)), nrow=dim(l)[1], byrow=T))
                    names(temp) <- c("beta","variance","pval")
                    temp[,"variables"] <- gsub("beta_","",names(beta))
                    temp$cancer = c
                    temp$n <- summary(fit)$n
                    df_hazard <- rbind(df_hazard,temp)
                }                           
            }
        }, 
            error=function (e) {
                print(paste("Error:",e))
                print(paste("Cancer",c,"was unsuccesful"))}
        )
    }
    if(all(dim(df_hazard)==c(0,0))) {
        return(list(NA,NA))
    }
    df_hazard_ret <- df_hazard
#     print(df_hazard)
    df_hazard$sd <- sqrt(df_hazard$variance)
    df_hazard <- df_hazard %>% filter(!is.na(beta)) %>% filter(!(abs(sd/beta)>100 & abs(beta)>0.01))
    df_hazard$p <- signif(df_hazard$pval,2)
#     df_hazard$p_prophaz <- signif(df_hazard$pval_prophaz,2)
#     df_hazard <- df_hazard %>% rowwise %>% mutate(variables = paste(sort(unlist(str_split(variables,":"))),collapse=":")) %>% ungroup
    res <- df_hazard %>% group_by(variables) %>% do(result=suppressWarnings(metagen(beta,sd,data=.,studlab=paste(cancer),comb.random=T,comb.fixed=T,method.tau = "SJ",hakn = TRUE,prediction=F,sm="SMD")))
    if(plot) {
        if (length(res$result)>0) {
            for (i in 1:length(res$result)) {
              forest(res$result[[i]],leftcols=c("studlab","TE","seTE","n","p"),print.tau2=F,print.pval.Q=F,print.I2=F)
                grid.text(res$variables[[i]], 0.5,0.95, gp=gpar(cex=1.5))
                grid.text(paste("p=",signif(res$result[[i]]$pval.fixed,2),sep=""), 0.3,0.225, gp=gpar(cex=1))
                grid.text(paste("p=",signif(res$result[[i]]$pval.random,2),sep=""), 0.3,0.175, gp=gpar(cex=1))
            }
        }
    }
    return(list(res,df_hazard_ret))
}   
                                      
meta_survival_AE_stratcancer <- function(cov_vec_in,OS_temp,vars_in_each=vie,vars_in_all=via,vars_direct=c(),transition="1|2|3",plot=TRUE,name=c(),minx=c(),maxx=c(),cancer_variable="CANCER_TYPE",cancer_leaveout=c()) {
    temp <- OS_temp %>% select(MRN,trans,one_of(vars_in_all)) %>% drop_na
    x <- model.matrix(as.formula(paste("~ MRN + trans + ",paste(vars_in_all,collapse=" + "),sep="")), temp)
    OS_temp <- OS_temp %>% select(-one_of(vars_in_all)) #%>% cbind(x)
    OS_temp <- merge(OS_temp,x,by=c("MRN","trans"))

    vars_trans1 <- vars_in_each$trans1
    vars_trans1 <- names(OS_temp)[grepl(paste(vars_trans1,collapse="|"),names(OS_temp)) & grepl("\\.1$",names(OS_temp))]
    vars_trans2 <- vars_in_each$trans2
    vars_trans2 <- names(OS_temp)[grepl(paste(vars_trans2,collapse="|"),names(OS_temp)) & grepl("\\.2$",names(OS_temp))]
    vars_trans3 <- vars_in_each$trans3
    vars_trans3 <- names(OS_temp)[grepl(paste(vars_trans3,collapse="|"),names(OS_temp)) & grepl("\\.3$",names(OS_temp))]

    cov_vec <- names(OS_temp)[grepl(paste("^",paste(cov_vec_in,collapse="\\.[1-3]$|^"),"\\.[1-3]$",sep=""),names(OS_temp)) & grepl(paste("\\.",transition,"$",sep=""),names(OS_temp))]
    vars_trans1 <- vars_trans1[!grepl(paste(cov_vec,collapse="|"),vars_trans1)]
    vars_trans2 <- vars_trans2[!grepl(paste(cov_vec,collapse="|"),vars_trans2)]
    vars_trans3 <- vars_trans3[!grepl(paste(cov_vec,collapse="|"),vars_trans3)]
    
    vars_in_all <- names(OS_temp)[grepl(paste(vars_in_all,collapse="|"),names(OS_temp))& !grepl("\\.[1-3]$",names(OS_temp))]
    vars_in_all <- vars_in_all[!grepl(paste(cov_vec,collapse="|"),vars_in_all)]
    
    vars <- c(vars_trans1,vars_trans2,vars_trans3,vars_in_all,vars_direct)
    
    OS_temp <- OS_temp %>% select(MRN,CANCER_TYPE,matches(cancer_variable),Tstart,Tstop,status,trans,one_of(vars),one_of(cov_vec)) %>% drop_na %>% filter(grepl(transition,as.character(trans)))
    
    df_hazard <- data.frame(cancer=c(),variables=c(),log_hr=c(),log_hr_sd=c(),pval=c(),n=c(),pval_prophaz=c())
    if (length(cancer_leaveout)>0) {
        cancers <- unique(OS_temp[,cancer_variable])[!grepl(paste(cancer_leaveout,collapse="|"),unique(OS_temp[,cancer_variable]))]
    } else {
        cancers <- unique(OS_temp[,cancer_variable])
    }
    for (c in cancers) {
        tryCatch({
            OS_C <- OS_temp %>% filter(get(cancer_variable)==c)
    #         OS_C <- OS_temp %>% filter(CANCER_TYPE==c) # %>% filter(trans==1)
            if (sum(OS_C$status)>0) {
                # make formula
                fit <- find_vars_fit(vars,cov_vec,OS_C,strat_cancer=T)
                vars_meta <- unlist(lapply(cov_vec,function(x) rownames(summary(fit)$coefficients)[grepl(x,rownames(summary(fit)$coefficients))]))
                temp <- data.table(cancer=c,variables=vars_meta,log_hr=summary(fit)$coefficients[vars_meta,c("coef")],log_hr_sd=summary(fit)$coefficients[vars_meta,c("se(coef)")],pval=summary(fit)$coefficients[vars_meta,c("Pr(>|z|)")])
                temp$n <- dim(OS_C %>% distinct(MRN))[[1]]
                temp$pval_prophaz <- suppressWarnings(cox.zph(fit))$table[vars_meta,"p"]
                df_hazard <- rbind(df_hazard,temp)
            }
        }, 
            error=function (e) {
                print(paste("Error:",e))
                print(paste("Cancer",c,"was unsuccesful"))}
        )
    }
    df_hazard_ret <- df_hazard                              
    df_hazard <- df_hazard %>% drop_na(pval_prophaz) %>% filter(pval_prophaz>prophazard_pval_cutoff_meta) %>% filter(!is.na(log_hr)) %>% filter(!(abs(log_hr_sd/log_hr)>100 & abs(log_hr)>0.01))

    df_hazard$p <- df_hazard$pval
    df_hazard$p_prophaz <- df_hazard$pval_prophaz
#     df_hazard$p <- signif(df_hazard$pval,2)
#     df_hazard$p_prophaz <- signif(df_hazard$pval_prophaz,2)
    res <- df_hazard %>% group_by(variables) %>% do(result=metagen(log_hr,log_hr_sd,data=.,studlab=paste(cancer),comb.random=F,comb.fixed=T,method.tau = "SJ",hakn = TRUE,prediction=F,sm="SMD"))
#     if(plot) {
#         if (length(res$result)>0) {
#             for (i in 1:length(res$result)) {
#      forest(res$result[[i]],leftcols=c("studlab","TE","seTE","n","p","p_prophaz"),print.tau2=F,print.pval.Q=F,print.I2=F,digits=4)
#                     grid.text(res$variables[[i]], 0.5,0.95, gp=gpar(cex=1.5))
#                     grid.text(paste("p=",signif(res$result[[i]]$pval.fixed,2),sep=""), 0.3,0.225, gp=gpar(cex=1))
#                     grid.text(paste("p=",signif(res$result[[i]]$pval.random,2),sep=""), 0.3,0.175, gp=gpar(cex=1))
#             }
#         }
#     }
    if(plot) {
        if (length(res$result)>0) {
            for (i in 1:length(res$result)) {
                res$result[[i]]$n_print <- as.character(as.integer(res$result[[i]]$data$n))
                res$result[[i]]$pval_print <- as.character(formatC(signif(as.double(res$result[[i]]$pval),2),digits=2,format="fg", flag="#"))
                res$result[[i]]$TE_print <- as.character(formatC(signif(as.double(res$result[[i]]$TE),2),digits=2,format="fg", flag="#"))
#                 print(names(res$result[[i]]$df.Q))
                fp <- grid.grabExpr(print(forest(res$result[[i]],colgap=unit(0,"points"),comb.random=F,leftcols=c("studlab","TE_print","n_print","pval_print"),leftlabs=c("Cancer Type","logHR","N","p"),rightcols=c("ci"),ref=0,print.tau2=F,print.pval.Q=F,print.I2=F,smlab="logHR",xlab="logHR",xlim=c(ifelse(length(minx)==0,ifelse(min(res$result[[i]]$TE) < -1,1.2,ifelse(min(res$result[[i]]$TE)>0,0.8,2))*min(res$result[[i]]$TE),minx[[i]]),ifelse(length(maxx)==0,ifelse(max(res$result[[i]]$TE)>1,1.2,ifelse(max(res$result[[i]]$TE)<0,0.8,2))*max(res$result[[i]]$TE),maxx[[i]])),text.fixed=paste("Fixed effects model:    logHR=",signif(res$result[[i]]$TE.fixed,2),", p=",signif(res$result[[i]]$pval.fixed,2),sep=""))))
#                 fp <- grid.grabExpr(print(forest(res$result[[i]],fontsize=24,plotwidth=unit(250, "points"),squaresize=0.8,colgap=unit(0,"points"),comb.random=F,leftcols=c("studlab","TE_print","n_print","pval_print"),leftlabs=c("Cancer Type","logHR","N","p"),rightcols=c("ci"),ref=0,print.tau2=F,print.pval.Q=F,print.I2=F,smlab="logHR",xlab="logHR",xlim=c(ifelse(length(min)==0,ifelse(min(res$result[[i]]$TE) < -1,1.2,ifelse(min(res$result[[i]]$TE)>0,0.8,2))*min(res$result[[i]]$TE),min[[i]]),ifelse(length(max)==0,ifelse(max(res$result[[i]]$TE)>1,1.2,ifelse(max(res$result[[i]]$TE)<0,0.8,2))*max(res$result[[i]]$TE),max[[i]])),text.fixed=paste("Fixed effects model:    logHR=",signif(res$result[[i]]$TE.fixed,2),", p=",signif(res$result[[i]]$pval.fixed,2),sep=""))))
                title <- textGrob(ifelse(length(name)==0,res$variables[[i]],name[[i]]), gp=gpar(fontface="bold"))
#                 png(paste(res$variables[[i]],".png",sep=""),width=1124,height=560,pointsize = 28)
                grid.arrange(top=title, fp)
#                 dev.off()
            }
        }
    }
    return(list(res,df_hazard_ret))
}                                      
                                          
make_multistate <- function(OS_temp) {
    if ("CBIO_PATIENT" %in% names(OS_temp)) {
        OS_temp <- OS_temp %>% select(-CBIO_PATIENT)
    }
    temp <- merge(OS_temp,fread("adverse_events_truenegcorrected.csv",na.strings="") %>% select(MRN,ADV_EVENT,ADV_EVENT_DT) %>% mutate(ADV_EVENT_DT=as.Date(ADV_EVENT_DT)),by=c("MRN")) %>% mutate(irAE_time=as.integer(ifelse(ADV_EVENT,as.integer(difftime(ADV_EVENT_DT,IO_START,unit="days")),tstop)),irAE=as.numeric(ifelse(ADV_EVENT,1,0))) %>% select(-ADV_EVENT_DT,-ADV_EVENT)
#     temp <- merge(OS_temp,fread("predcur.csv",na.strings="") %>% select(MRN,ADV_EVENT,ADV_EVENT_DT) %>% mutate(ADV_EVENT_DT=as.Date(ADV_EVENT_DT)),by=c("MRN")) %>% mutate(irAE_time=as.integer(ifelse(ADV_EVENT,as.integer(difftime(ADV_EVENT_DT,IO_START,unit="days")),tstop)),irAE=as.numeric(ifelse(ADV_EVENT,1,0))) %>% select(-ADV_EVENT_DT,-ADV_EVENT)
    temp <- temp %>% mutate(irAE_time=as.integer(ifelse(irAE_time==0,1,irAE_time)),tstop=as.integer(ifelse(tstop<=irAE_time,irAE_time,tstop)))
    temp <- temp %>% drop_na(irAE_time)
    tmat <- transMat(x = list(c(2,3), c(3), c()), names = c("Treatment", "AE", "Death"))
    OS_mult <- msprep(data = temp, trans = tmat, time = c(NA, "irAE_time", "tstop"), status = c(NA, "irAE", "event"), id="MRN")
    OS_mult <- merge(OS_mult,temp,by="MRN") %>% rowwise %>% mutate(Tstart=max(tstart,Tstart),time=Tstop-Tstart) %>% select(-irAE_time,-tstart,-tstop,-event,-irAE,-id)
    OS_mult <- OS_mult %>% filter(Tstop>Tstart)

    attributes(OS_mult)$class <- c("msdata","data.frame")
    attributes(OS_mult)$trans <- tmat

    print(events(OS_mult))
    OS_mult <- OS_mult %>% select_if(function(x) ifelse(is.factor(x),length(levels(x))>=2,ifelse(is.character(x),length(levels(as.factor(x)))>=2,TRUE)))
    covs <- names(OS_mult %>% select(-MRN,-from,-to,-trans,-Tstart,-Tstop,-time,-status,-IO_START,-IO_STOP,-AGE_AT_TREATMENTSTART_BINNED,-LAST_DATE,-PRIMARY_CANCER_DIAGNOSIS,-CANCER_TYPE,-othertreatment_after_io,-othertreatment_while_io,-BIOPSY_SITE_TYPE,-MET_BEFORE_TREATMENTSTART,-before_2018))
    OS_mult <- expand.covs(OS_mult, covs, longnames = T)
    return(OS_mult)
}
                                     
make_multistate_noleft <- function(OS_temp) {
    if ("CBIO_PATIENT" %in% names(OS_temp)) {
        OS_temp <- OS_temp %>% select(-CBIO_PATIENT)
    }
    OS_temp$tstart <- 0
    temp <- merge(OS_temp,fread("adverse_events_truenegcorrected.csv",na.strings="") %>% select(MRN,ADV_EVENT,ADV_EVENT_DT) %>% mutate(ADV_EVENT_DT=as.Date(ADV_EVENT_DT)),by=c("MRN")) %>% mutate(irAE_time=as.integer(ifelse(ADV_EVENT,as.integer(difftime(ADV_EVENT_DT,IO_START,unit="days")),tstop)),irAE=as.numeric(ifelse(ADV_EVENT,1,0))) %>% select(-ADV_EVENT_DT,-ADV_EVENT)
#     temp <- merge(OS_temp,fread("predcur.csv",na.strings="") %>% select(MRN,ADV_EVENT,ADV_EVENT_DT) %>% mutate(ADV_EVENT_DT=as.Date(ADV_EVENT_DT)),by=c("MRN")) %>% mutate(irAE_time=as.integer(ifelse(ADV_EVENT,as.integer(difftime(ADV_EVENT_DT,IO_START,unit="days")),tstop)),irAE=as.numeric(ifelse(ADV_EVENT,1,0))) %>% select(-ADV_EVENT_DT,-ADV_EVENT)
    temp <- temp %>% mutate(irAE_time=as.integer(ifelse(irAE_time==0,1,irAE_time)),tstop=as.integer(ifelse(tstop<=irAE_time,irAE_time,tstop)))
    temp <- temp %>% drop_na(irAE_time)
    tmat <- transMat(x = list(c(2,3), c(3), c()), names = c("Treatment", "AE", "Death"))
    OS_mult <- msprep(data = temp, trans = tmat, time = c(NA, "irAE_time", "tstop"), status = c(NA, "irAE", "event"), id="MRN")
    OS_mult <- merge(OS_mult,temp,by="MRN") %>% rowwise %>% mutate(Tstart=max(tstart,Tstart),time=Tstop-Tstart) %>% select(-irAE_time,-tstart,-tstop,-event,-irAE,-id)
    OS_mult <- OS_mult %>% filter(Tstop>Tstart)

    attributes(OS_mult)$class <- c("msdata","data.frame")
    attributes(OS_mult)$trans <- tmat

    print(events(OS_mult))
    OS_mult <- OS_mult %>% select_if(function(x) ifelse(is.factor(x),length(levels(x))>=2,ifelse(is.character(x),length(levels(as.factor(x)))>=2,TRUE)))
    covs <- names(OS_mult %>% select(-MRN,-from,-to,-trans,-Tstart,-Tstop,-time,-status,-IO_START,-IO_STOP,-AGE_AT_TREATMENTSTART_BINNED,-LAST_DATE,-PRIMARY_CANCER_DIAGNOSIS,-CANCER_TYPE,-othertreatment_after_io,-othertreatment_while_io,-BIOPSY_SITE_TYPE,-MET_BEFORE_TREATMENTSTART,-before_2018))
    OS_mult <- expand.covs(OS_mult, covs, longnames = T)
    return(OS_mult)
}                                     
                                     
make_multistate_cur <- function(OS_temp) {
    if ("CBIO_PATIENT" %in% names(OS_temp)) {
        OS_temp <- OS_temp %>% select(-CBIO_PATIENT)
    }
    temp <- merge(OS_temp,fread("curated_irAE_full.csv",na.strings="") %>% select(MRN,ADV_EVENT,ADV_EVENT_DT) %>% mutate(ADV_EVENT_DT=as.Date(ADV_EVENT_DT)),by=c("MRN")) %>% mutate(irAE_time=as.integer(ifelse(ADV_EVENT,as.integer(difftime(ADV_EVENT_DT,IO_START,unit="days")),tstop)),irAE=as.numeric(ifelse(ADV_EVENT,1,0))) %>% select(-ADV_EVENT_DT,-ADV_EVENT)
    temp <- temp %>% mutate(irAE_time=as.integer(ifelse(irAE_time==0,1,irAE_time)),tstop=as.integer(ifelse(tstop<=irAE_time,irAE_time,tstop)))
    temp <- temp %>% drop_na(irAE_time)
    tmat <- transMat(x = list(c(2,3), c(3), c()), names = c("Treatment", "AE", "Death"))
    OS_mult <- msprep(data = temp, trans = tmat, time = c(NA, "irAE_time", "tstop"), status = c(NA, "irAE", "event"), id="MRN")
    OS_mult <- merge(OS_mult,temp,by="MRN") %>% rowwise %>% mutate(Tstart=max(tstart,Tstart),time=Tstop-Tstart) %>% select(-irAE_time,-tstart,-tstop,-event,-irAE,-id)
    OS_mult <- OS_mult %>% filter(Tstop>Tstart)

    attributes(OS_mult)$class <- c("msdata","data.frame")
    attributes(OS_mult)$trans <- tmat

    print(events(OS_mult))
    OS_mult <- OS_mult %>% select_if(function(x) ifelse(is.factor(x),length(levels(x))>=2,ifelse(is.character(x),length(levels(as.factor(x)))>=2,TRUE)))
    covs <- names(OS_mult %>% select(-MRN,-from,-to,-trans,-Tstart,-Tstop,-time,-status,-IO_START,-IO_STOP,-AGE_AT_TREATMENTSTART_BINNED,-LAST_DATE,-PRIMARY_CANCER_DIAGNOSIS,-CANCER_TYPE,-othertreatment_after_io,-othertreatment_while_io,-BIOPSY_SITE_TYPE,-MET_BEFORE_TREATMENTSTART,-before_2018))
    OS_mult <- expand.covs(OS_mult, covs, longnames = T)
    return(OS_mult)
}
                                     
make_multistate_trueneg <- function(OS_temp) {
    if ("CBIO_PATIENT" %in% names(OS_temp)) {
        OS_temp <- OS_temp %>% select(-CBIO_PATIENT)
    }
    temp <- merge(OS_temp,fread("adverse_events_trueneg.csv",na.strings="") %>% select(MRN,true_negative,ADV_EVENT_DT) %>% rename(ADV_EVENT=true_negative) %>% mutate(ADV_EVENT_DT=as.Date(ADV_EVENT_DT)),by=c("MRN")) %>% mutate(irAE_time=as.integer(ifelse(ADV_EVENT,as.integer(difftime(ADV_EVENT_DT,IO_START,unit="days")),tstop)),irAE=as.numeric(ifelse(ADV_EVENT,1,0))) %>% select(-ADV_EVENT_DT,-ADV_EVENT)
#     temp <- merge(OS_temp,fread("predcur.csv",na.strings="") %>% select(MRN,ADV_EVENT,ADV_EVENT_DT) %>% mutate(ADV_EVENT_DT=as.Date(ADV_EVENT_DT)),by=c("MRN")) %>% mutate(irAE_time=as.integer(ifelse(ADV_EVENT,as.integer(difftime(ADV_EVENT_DT,IO_START,unit="days")),tstop)),irAE=as.numeric(ifelse(ADV_EVENT,1,0))) %>% select(-ADV_EVENT_DT,-ADV_EVENT)
    temp <- temp %>% mutate(irAE_time=as.integer(ifelse(irAE_time==0,1,irAE_time)),tstop=as.integer(ifelse(tstop<=irAE_time,irAE_time,tstop)))
    temp <- temp %>% drop_na(irAE_time)
    tmat <- transMat(x = list(c(2,3), c(3), c()), names = c("Treatment", "AE", "Death"))
    OS_mult <- msprep(data = temp, trans = tmat, time = c(NA, "irAE_time", "tstop"), status = c(NA, "irAE", "event"), id="MRN")
    OS_mult <- merge(OS_mult,temp,by="MRN") %>% rowwise %>% mutate(Tstart=max(tstart,Tstart),time=Tstop-Tstart) %>% select(-irAE_time,-tstart,-tstop,-event,-irAE,-id)
    OS_mult <- OS_mult %>% filter(Tstop>Tstart)

    attributes(OS_mult)$class <- c("msdata","data.frame")
    attributes(OS_mult)$trans <- tmat

    print(events(OS_mult))
    OS_mult <- OS_mult %>% select_if(function(x) ifelse(is.factor(x),length(levels(x))>=2,ifelse(is.character(x),length(levels(as.factor(x)))>=2,TRUE)))
    covs <- names(OS_mult %>% select(-MRN,-from,-to,-trans,-Tstart,-Tstop,-time,-status,-IO_START,-IO_STOP,-AGE_AT_TREATMENTSTART_BINNED,-LAST_DATE,-PRIMARY_CANCER_DIAGNOSIS,-CANCER_TYPE,-othertreatment_after_io,-othertreatment_while_io,-BIOPSY_SITE_TYPE,-MET_BEFORE_TREATMENTSTART,-before_2018))
    OS_mult <- expand.covs(OS_mult, covs, longnames = T)
    return(OS_mult)
}       
                                     
make_multistate_type <- function(OS_temp,type) {
    if ("CBIO_PATIENT" %in% names(OS_temp)) {
        OS_temp <- OS_temp %>% select(-CBIO_PATIENT)
    }
    temp <- inner_join(OS_temp,fread("adverse_events_truenegcorrected_type.csv",na.strings="")  %>% filter(is.na(TYPE) | grepl(type,tolower(TYPE))) %>% select(MRN,ADV_EVENT,ADV_EVENT_DT) %>% mutate(ADV_EVENT_DT=as.Date(ADV_EVENT_DT)),by=c("MRN")) %>% mutate(irAE_time=as.integer(ifelse(ADV_EVENT,as.integer(difftime(ADV_EVENT_DT,IO_START,unit="days")),tstop)),irAE=as.numeric(ifelse(ADV_EVENT,1,0))) %>% select(-ADV_EVENT_DT,-ADV_EVENT)
#     temp <- merge(OS_temp,fread("predcur.csv",na.strings="") %>% select(MRN,ADV_EVENT,ADV_EVENT_DT) %>% mutate(ADV_EVENT_DT=as.Date(ADV_EVENT_DT)),by=c("MRN")) %>% mutate(irAE_time=as.integer(ifelse(ADV_EVENT,as.integer(difftime(ADV_EVENT_DT,IO_START,unit="days")),tstop)),irAE=as.numeric(ifelse(ADV_EVENT,1,0))) %>% select(-ADV_EVENT_DT,-ADV_EVENT)
    temp <- temp %>% mutate(irAE_time=as.integer(ifelse(irAE_time==0,1,irAE_time)),tstop=as.integer(ifelse(tstop<=irAE_time,irAE_time,tstop)))
    temp <- temp %>% drop_na(irAE_time)
    tmat <- transMat(x = list(c(2,3), c(3), c()), names = c("Treatment", "AE", "Death"))
    OS_mult <- msprep(data = temp, trans = tmat, time = c(NA, "irAE_time", "tstop"), status = c(NA, "irAE", "event"), id="MRN")
    OS_mult <- merge(OS_mult,temp,by="MRN") %>% rowwise %>% mutate(Tstart=max(tstart,Tstart),time=Tstop-Tstart) %>% select(-irAE_time,-tstart,-tstop,-event,-irAE,-id)
    OS_mult <- OS_mult %>% filter(Tstop>Tstart)

    attributes(OS_mult)$class <- c("msdata","data.frame")
    attributes(OS_mult)$trans <- tmat

    print(events(OS_mult))
    OS_mult <- OS_mult %>% select_if(function(x) ifelse(is.factor(x),length(levels(x))>=2,ifelse(is.character(x),length(levels(as.factor(x)))>=2,TRUE)))
    covs <- names(OS_mult %>% select(-MRN,-from,-to,-trans,-Tstart,-Tstop,-time,-status,-IO_START,-IO_STOP,-AGE_AT_TREATMENTSTART_BINNED,-LAST_DATE,-PRIMARY_CANCER_DIAGNOSIS,-CANCER_TYPE,-othertreatment_after_io,-othertreatment_while_io,-BIOPSY_SITE_TYPE,-MET_BEFORE_TREATMENTSTART,-before_2018))
    OS_mult <- expand.covs(OS_mult, covs, longnames = T)
    return(OS_mult)
}                                        