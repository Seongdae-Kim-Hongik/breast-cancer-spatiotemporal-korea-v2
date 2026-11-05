# 88h — 공통 사후예측 점검·분산 분해(88e 정신질환판 일반화, 2026-10-21 유방암부터). E 변형 민감도는 질환별 E 파일이 있을 때만 별도.
#  사용(01_스크립트/): TAG=BRCA ST_VER=std Rscript 88h_ppc_varpart.R    산출: results_<VER>/S42_ppc.csv, S42_pit_hist.csv, S42_variance_partition.csv
suppressPackageStartupMessages({library(INLA);library(sf);library(spdep);library(jsonlite)})
TAG<-Sys.getenv("TAG");stopifnot(nzchar(TAG));VER<-Sys.getenv("ST_VER",unset="std");SUF<-if(nzchar(VER))paste0("_",VER) else ""
ROOT<-normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE))),".."))
SHP<-"${HOME}/Library/CloudStorage/[REDACTED_EMAIL]/내 드라이브/S.K/기타/R_Studio/FBD_DATA_ZIP/final.shp"
OUTD<-file.path(ROOT,"04_분석",TAG);RES<-file.path(OUTD,paste0("results",SUF))
inla.setOption(num.threads="1:1")
LOG<-function(...)cat(sprintf("[88h %s] ",format(Sys.time(),"%H:%M:%S")),sprintf(...),"\n",sep="")
S<-fromJSON(file.path(OUTD,sprintf("config_%s%s.json",TAG,SUF)),simplifyVector=TRUE)
INFO<-fromJSON(file.path(OUTD,sprintf("panel_INFO_%s%s.json",TAG,SUF)))
DD<-read.csv(file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);DD$fbd_region<-enc2utf8(DD$fbd_region)
HB<-readLines(file.path(OUTD,sprintf("HM_best_%s%s.txt",TAG,SUF)));kv<-function(k)sub(paste0("^",k,"="),"",grep(paste0("^",k,"="),HB,value=TRUE)[1])
MAINFAM<-kv("main_family");BEST<-kv("best_set");FORCED<-S$forced;SETS<-lapply(S$sets_after_vif,as.character);HYPV<-SETS[[BEST]];MAINV<-unique(c(HYPV,FORCED))
ROB<-read.csv(file.path(RES,"S7_multiverse_summary.csv"),stringsAsFactors=FALSE);ROBV<-ROB$cov[ROB$robust]
LOG("robust %d: %s",length(ROBV),paste(ROBV,collapse=";"))
M3<-read.csv(file.path(RES,"S3_main_IRR.csv"),stringsAsFactors=FALSE);msign<-setNames(sign(log(M3$IRR)),M3$cov)

G<-st_read(SHP,quiet=TRUE);G$region<-enc2utf8(G$region);n<-nrow(G)
g<-inla.read.graph(file.path(RES,sprintf("nb_%s.graph",TAG)));stopifnot(g$n==n)
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
m<-inla(form,family=MAINFAM,data=DD,E=E,control.fixed=CF,control.compute=list(cpo=TRUE,waic=TRUE,config=TRUE),
        control.inla=list(control.vb=list(emergency=50)),control.predictor=list(compute=TRUE),verbose=FALSE)
LOG("main refit WAIC %.1f; cpo fail %d",m$waic$waic,sum(m$cpo$failure>0,na.rm=TRUE))
pit<-m$cpo$pit; br<-seq(0,1,0.1); h<-hist(pit,breaks=br,plot=FALSE)
write.csv(data.frame(bin_lo=br[-11],bin_hi=br[-1],count=h$counts,share=round(h$counts/sum(h$counts),4)),file.path(RES,"S42_pit_hist.csv"),row.names=FALSE)
NS<-as.integer(Sys.getenv("NS","200"));S<-inla.posterior.sample(NS,m,num.threads="1:1")
n<-nrow(DD);hn<-names(S[[1]]$hyperpar);szn<-grep("size",hn,value=TRUE)[1]
Y<-matrix(NA,n,NS)
for(s in seq_len(NS)){eta<-S[[s]]$latent[grep("^Predictor:",rownames(S[[s]]$latent)),1];sz<-S[[s]]$hyperpar[szn];Y[,s]<-if(is.na(szn))rpois(n,DD$E*exp(eta)) else rnbinom(n,size=sz,mu=DD$E*exp(eta))}
lo<-apply(Y,1,quantile,0.025);hi<-apply(Y,1,quantile,0.975);lo50<-apply(Y,1,quantile,0.25);hi50<-apply(Y,1,quantile,0.75)
cov95<-mean(DD$cases>=lo&DD$cases<=hi);cov50<-mean(DD$cases>=lo50&DD$cases<=hi50)
ppp_var<-mean(apply(Y,2,function(y)var(log((y+0.5)/DD$E)))>=var(log((DD$cases+0.5)/DD$E)))
ksp<-suppressWarnings(ks.test(pit,"punif")$p.value)
ppc<-data.frame(n_obs=n,n_samples=NS,coverage95=round(cov95,4),coverage50=round(cov50,4),ppp_variance_logSIR=round(ppp_var,3),
  pit_ks_p=signif(ksp,3),pit_mean=round(mean(pit,na.rm=TRUE),3),pit_sd=round(sd(pit,na.rm=TRUE),3),LCPO=round(-mean(log(m$cpo$cpo),na.rm=TRUE),4))
if(NS>=200){write.csv(ppc,file.path(RES,"S42_ppc.csv"),row.names=FALSE)};print(ppc)
## 분산 분해(선형예측자, 사후평균 기준)
X<-model.matrix(as.formula(paste("~",paste(bt(MAINV),collapse="+"))),DD);b<-m$summary.fixed$mean;names(b)<-rownames(m$summary.fixed)
colnames(X)<-gsub("`","",colnames(X));fx<-as.vector(X[,-1,drop=FALSE]%*%b[colnames(X)[-1]])
ra<-m$summary.random$area$mean;comb<-ra[1:g$n][DD$area];ustar<-ra[(g$n+1):(2*g$n)][DD$area]
phi<-m$summary.hyperpar[grep("Phi",rownames(m$summary.hyperpar)),"mean"];tau<-m$summary.hyperpar[grep("Precision for area",rownames(m$summary.hyperpar)),"mean"]
str<-sqrt(phi/tau)*ustar;unstr<-comb-str
tt<-m$summary.random$time$mean[DD$time]+m$summary.random$yiid$mean[DD$yiid]
fxh<-as.vector(X[,bt2<-colnames(X)[-1][colnames(X)[-1]%in%HYPV],drop=FALSE]%*%b[bt2]);fxf<-fx-fxh
comp<-list(fixed_covariates=fx,fixed_hypothesis=fxh,fixed_forced_controls=fxf,spatial_structured=str,spatial_unstructured=unstr,year_effects=tt)
v<-sapply(comp,var);vsum<-sum(v[c("fixed_covariates","spatial_structured","spatial_unstructured","year_effects")]);tot<-var(m$summary.linear.predictor$mean)
vp<-data.frame(component=names(v),variance=signif(v,4),share_of_sum=round(v/vsum,3),share_of_linear_predictor_var=round(v/tot,3))
vp<-rbind(vp,data.frame(component="linear_predictor_total",variance=signif(tot,4),share_of_sum=NA,share_of_linear_predictor_var=1))
write.csv(vp,file.path(RES,"S42_variance_partition.csv"),row.names=FALSE);print(vp)
LOG("done")
