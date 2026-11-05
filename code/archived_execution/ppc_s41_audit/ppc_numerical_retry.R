# 88h — 공통 사후예측 점검·분산 분해(88e 정신질환판 일반화, 2026-10-21 유방암부터). E 변형 민감도는 질환별 E 파일이 있을 때만 별도.
#  사용(01_스크립트/): TAG=BRCA ST_VER=std Rscript 88h_ppc_varpart.R    산출: results_<VER>/S42_ppc.csv, S42_pit_hist.csv, S42_variance_partition.csv
suppressPackageStartupMessages({library(INLA);library(sf);library(spdep);library(jsonlite)})
TAG<-Sys.getenv("TAG");stopifnot(nzchar(TAG));VER<-Sys.getenv("ST_VER",unset="std");SUF<-if(nzchar(VER))paste0("_",VER) else ""
ROOT<-normalizePath(Sys.getenv("FINAL_ROOT"),mustWork=TRUE)
SHP<-"${HOME}/Library/CloudStorage/[REDACTED_EMAIL]/내 드라이브/S.K/기타/R_Studio/FBD_DATA_ZIP/final.shp"
OUTD<-file.path(ROOT,"04_분석",TAG);SRCRES<-file.path(OUTD,paste0("results",SUF));RES<-Sys.getenv("PPC_OUT");dir.create(RES,recursive=TRUE,showWarnings=FALSE)
inla.setOption(num.threads="1:1")
LOG<-function(...)cat(sprintf("[88h %s] ",format(Sys.time(),"%H:%M:%S")),sprintf(...),"\n",sep="")
S<-fromJSON(file.path(OUTD,sprintf("config_%s%s.json",TAG,SUF)),simplifyVector=TRUE)
INFO<-fromJSON(file.path(OUTD,sprintf("panel_INFO_%s%s.json",TAG,SUF)))
DD<-read.csv(file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);DD$fbd_region<-enc2utf8(DD$fbd_region)
HB<-readLines(file.path(OUTD,sprintf("HM_best_%s%s.txt",TAG,SUF)));kv<-function(k)sub(paste0("^",k,"="),"",grep(paste0("^",k,"="),HB,value=TRUE)[1])
MAINFAM<-kv("main_family");BEST<-kv("best_set");FORCED<-S$forced;SETS<-lapply(S$sets_after_vif,as.character);HYPV<-SETS[[BEST]];MAINV<-unique(c(HYPV,FORCED))
ROB<-read.csv(file.path(SRCRES,"S7_multiverse_summary.csv"),stringsAsFactors=FALSE);rob_flag<-tolower(trimws(as.character(ROB$robust)));stopifnot(all(rob_flag %in% c("true","false")));ROBV<-ROB$cov[rob_flag=="true"];stopifnot(length(ROBV)==4L,!anyNA(ROBV));writeLines(ROBV,file.path(RES,"actual_robust_variables.txt"))
LOG("robust %d: %s",length(ROBV),paste(ROBV,collapse=";"))
M3<-read.csv(file.path(SRCRES,"S3_main_IRR.csv"),stringsAsFactors=FALSE);msign<-setNames(sign(log(M3$IRR)),M3$cov)

G<-st_read(SHP,quiet=TRUE);G$region<-enc2utf8(G$region);n<-nrow(G)
GF<-file.path(ROOT,"input_graph.graph")
g<-inla.read.graph(GF);stopifnot(g$n==n)
stopifnot(digest::digest(file=GF,algo="sha256")==digest::digest(file=file.path(SRCRES,"nb_BRCA.graph"),algo="sha256"))
stopifnot(length(MAINV)==24L,setequal(MAINV,M3$cov),MAINFAM=="nbinomial")
file.copy(GF,file.path(RES,"actual_input_graph.graph"),overwrite=FALSE)
writeLines(MAINV,file.path(RES,"actual_MAINV.txt"))
DD<-DD[order(DD$year,match(DD$fbd_region,G$region)),];DD$area<-match(DD$fbd_region,G$region);Y0<-min(DD$year);DD$time<-DD$year-Y0+1;DD$yiid<-DD$time;DD$id<-seq_len(nrow(DD))
pcU<-function(u=0.5,a=0.01,pu=0.5,pa=2/3)list(prec=list(prior="pc.prec",param=c(u,a)),phi=list(prior="pc",param=c(pu,pa)))
pc<-pcU();pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
SPf<-function(gr="g",hyper="pc",extra="")sprintf("f(area,model='bym2',graph=%s,scale.model=TRUE,constr=TRUE,hyper=%s%s)",gr,hyper,extra)
TI<-"f(time,model='rw1',scale.model=TRUE,hyper=pcp)";IT<-"f(yiid,model='iid',hyper=pcp)"
bt<-function(v)sprintf("`%s`",v)
fit<-function(vars,data=DD,sp=SPf(),extra=character(0),robust=FALSE,env=parent.frame()){
  rhs<-c("1",bt(vars),extra,sp,TI,IT);rhs<-rhs[nzchar(rhs)];form<-as.formula(paste("cases ~",paste(rhs,collapse="+")),env=env)
  ci<-if(robust)list(control.vb=list(emergency=50),int.strategy="eb",strategy="gaussian") else list(control.vb=list(emergency=50))
  inla(form,family=MAINFAM,data=data,E=E,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE,cpo=FALSE,config=FALSE),control.inla=ci,verbose=FALSE)}
safe<-function(...){tryCatch(fit(...),error=function(e){message("  fit failed: ",conditionMessage(e)," -> robust");tryCatch(fit(...,robust=TRUE),error=function(e2){message("  robust failed: ",conditionMessage(e2));NULL})})}
irr<-function(m,v){r<-m$summary.fixed[v,,drop=FALSE];pg<-1-inla.pmarginal(0,m$marginals.fixed[[v]])
  data.frame(cov=v,IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),P_gt1=round(pg,3),cred=(r$`0.025quant`>0|r$`0.975quant`<0),stringsAsFactors=FALSE)}
set.seed(20261020)
## (Q3) 주모형 재적합: CPO/PIT + 사후표본
rhs<-c("1",bt(MAINV),SPf(),TI,IT);form<-as.formula(paste("cases ~",paste(rhs,collapse="+")))
initial_fit<-readRDS(Sys.getenv("PPC_INITIAL_RDS"))
stopifnot(isTRUE(initial_fit$ok),all(is.finite(initial_fit$mode$theta)))
m<-inla(form,family=MAINFAM,data=DD,E=E,control.fixed=CF,control.compute=list(cpo=TRUE,waic=TRUE,config=TRUE),
        control.mode=list(result=initial_fit,restart=TRUE),
        control.inla=list(control.vb=list(emergency=50),tolerance=1e-6,tolerance.f=1e-6,tolerance.g=1e-6,tolerance.x=1e-6),
        control.predictor=list(compute=TRUE,link=1),verbose=FALSE)
LOG("main refit WAIC %.1f; cpo fail %d",m$waic$waic,sum(m$cpo$failure>0,na.rm=TRUE))
saveRDS(m,file.path(RES,"model_configTRUE.rds"))
diag<-list(ok=m$ok,mode_status=m$mode$mode.status,misc_mode_status=m$misc$mode.status,internal_warnings=m$misc$warnings,
 fixed_finite=all(is.finite(as.matrix(m$summary.fixed))),hyper_finite=all(is.finite(as.matrix(m$summary.hyperpar))),
 waic=m$waic$waic,waic_finite=is.finite(m$waic$waic),strategy="default_integration_strict_optimization",initialization=Sys.getenv("PPC_INITIAL_RDS"),tolerance=1e-6,config=TRUE,threads="1:1")
write_json(diag,file.path(RES,"PPC_numerical_diagnostics.json"),auto_unbox=TRUE,pretty=TRUE,null="null",digits=NA)
LOG("numerical diagnostics: %s",toJSON(diag,auto_unbox=TRUE,null="null",digits=NA))
stopifnot(diag$fixed_finite,diag$hyper_finite,diag$waic_finite)
if(!is.null(diag$ok))stopifnot(isTRUE(diag$ok))
for(x in list(diag$mode_status,diag$misc_mode_status))if(length(x))stopifnot(all(x==0))
writeLines(if(length(unlist(diag$internal_warnings))) "requires_warning_review" else "no_internal_warnings_recorded",file.path(RES,"PPC_diagnostic_review_status.txt"))
write.csv(DD,file.path(RES,"indexed_model_input.csv"),row.names=FALSE)
write.csv(data.frame(variable=rownames(m$summary.fixed),m$summary.fixed,check.names=FALSE),file.path(RES,"beta_unrounded.csv"),row.names=FALSE)
write.csv(data.frame(year=DD$year,fbd_region=DD$fbd_region,waic=m$waic$local.waic),file.path(RES,"pointwise_waic.csv"),row.names=FALSE)
pit<-m$cpo$pit;stopifnot(all(is.finite(pit)),all(pit>=0 & pit<=1)); br<-seq(0,1,0.1); h<-hist(pit,breaks=br,plot=FALSE)
write.csv(data.frame(bin_lo=br[-11],bin_hi=br[-1],count=h$counts,share=round(h$counts/sum(h$counts),4)),file.path(RES,"S42_pit_hist.csv"),row.names=FALSE)
NS<-as.integer(Sys.getenv("NS","200"));S<-inla.posterior.sample(NS,m,num.threads="1:1",seed=20261104)
saveRDS(S,file.path(RES,"joint_posterior_draws_improved.rds"))
n<-nrow(DD);hn<-names(S[[1]]$hyperpar);szn<-grep("size",hn,value=TRUE);stopifnot(length(szn)==1L)
Y<-matrix(NA,n,NS)
for(s in seq_len(NS)){
 ids<-grep("^Predictor:",rownames(S[[s]]$latent));ord<-as.integer(sub("^Predictor:","",rownames(S[[s]]$latent)[ids]))
 stopifnot(length(ids)==n,setequal(ord,seq_len(n)))
 eta<-S[[s]]$latent[ids[match(seq_len(n),ord)],1];sz<-S[[s]]$hyperpar[szn]
 stopifnot(is.finite(sz),sz>0,all(is.finite(eta)))
 Y[,s]<-rnbinom(n,size=sz,mu=DD$E*exp(eta))
}
saveRDS(list(predictive_counts=Y,keys=DD[,c("year","fbd_region")],seed=20261104,NS=NS),file.path(RES,"predictive_draws.rds"))

lo<-apply(Y,1,quantile,0.025);hi<-apply(Y,1,quantile,0.975);lo50<-apply(Y,1,quantile,0.25);hi50<-apply(Y,1,quantile,0.75)
cov95<-mean(DD$cases>=lo&DD$cases<=hi);cov50<-mean(DD$cases>=lo50&DD$cases<=hi50)
ppp_var<-mean(apply(Y,2,function(y)var(log((y+0.5)/DD$E)))>=var(log((DD$cases+0.5)/DD$E)))
ksp<-suppressWarnings(ks.test(pit,"punif")$p.value)
ppc<-data.frame(n_obs=n,n_samples=NS,coverage95=cov95,coverage50=cov50,ppp_variance_logSIR=ppp_var,
 pit_ks_p_uncalibrated=ksp,pit_mean=mean(pit),pit_sd=sd(pit),LCPO=-mean(log(m$cpo$cpo)),cpo_failed=sum(m$cpo$failure>0,na.rm=TRUE),WAIC=m$waic$waic)
write.csv(ppc,file.path(RES,"S42_ppc.csv"),row.names=FALSE)
write.csv(data.frame(year=DD$year,fbd_region=DD$fbd_region,cases=DD$cases,E=DD$E,pred_lo95=lo,pred_hi95=hi,pred_lo50=lo50,pred_hi50=hi,pit=pit,cpo=m$cpo$cpo,cpo_failure=m$cpo$failure),file.path(RES,"PPC_pointwise.csv"),row.names=FALSE)
# Descriptive posterior-mean component variation; NOT R2 or explained variance.
X<-model.matrix(as.formula(paste("~",paste(bt(MAINV),collapse="+"))),DD)
b<-setNames(m$summary.fixed$mean,rownames(m$summary.fixed));colnames(X)<-gsub("`","",colnames(X))
fx<-rowSums(sweep(X[,-1,drop=FALSE],2,b[colnames(X)[-1]],"*"))
sp<-m$summary.random$area$mean[DD$area]
tt<-m$summary.random$time$mean[DD$time]+m$summary.random$yiid$mean[DD$yiid]
components<-data.frame(fixed_covariates=fx,spatial_combined=sp,year_effects=tt)
v<-sapply(components,var);total<-var(m$summary.linear.predictor$mean)
write.csv(data.frame(component=names(v),variance=v,descriptive_share_of_component_variance_sum=v/sum(v),linear_predictor_variance=total,interpretation="Posterior-mean descriptive variation, not variance explained; covariances nonzero"),file.path(RES,"S42_variance_partition.csv"),row.names=FALSE)
write.csv(cov(components),file.path(RES,"component_covariance.csv"))
delta<-max(abs(rowSums(components)+b["(Intercept)"]-m$summary.linear.predictor$mean))
write_json(list(n=nrow(DD),NS=NS,main_variables=MAINV,family=MAINFAM,config=TRUE,threads="1:1",graph=GF,WAIC=m$waic$waic,component_linear_predictor_max_error=delta,structured_unstructured_split="not reported: avoids plug-in hyperparameter decomposition",inla_version=as.character(packageVersion("INLA"))),file.path(RES,"PPC_method_metadata.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
# Marginal-summary means need not add exactly to the separately approximated predictor mean.
# Keep the discrepancy explicit; do not call this a variance-explained decomposition.
residual<-m$summary.linear.predictor$mean-(rowSums(components)+b["(Intercept)"])
write.csv(data.frame(year=DD$year,fbd_region=DD$fbd_region,component_sum=rowSums(components)+b["(Intercept)"],linear_predictor_mean=m$summary.linear.predictor$mean,approximation_residual=residual),file.path(RES,"component_reconstruction_audit.csv"),row.names=FALSE)
stopifnot(is.finite(delta),all(is.finite(as.matrix(components))))
write_json(list(numerical_warnings=diag$internal_warnings,component_max_abs_residual=delta,original_tolerance=1e-5,component_check_passed=delta<1e-5,scientific_validation=if(length(unlist(diag$internal_warnings))||delta>=1e-5) "requires_review" else "diagnostic_checks_passed",note="Computational completion does not imply scientific validation; residual retained without assigning it to a substantive component"),file.path(RES,"PPC_review_status.json"),auto_unbox=TRUE,pretty=TRUE,null="null",digits=NA)
writeLines(capture.output(sessionInfo()),file.path(RES,"sessionInfo.txt"))
writeLines("PPC complete",file.path(RES,"COMPLETE"))
LOG("done")
