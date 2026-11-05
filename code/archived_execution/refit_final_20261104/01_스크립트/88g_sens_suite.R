# 88g — 공간 트랙 질환 공통 민감도 묶음(88c·88f 정신질환판 일반화, 2026-10-21 유방암 2호부터).
#  사용(01_스크립트/): TAG=BRCA ST_VER=std GROUP=1..6 [NCO=GCA,INFECT] [FOCAL=walking_practice_rate] [SMOOTH=sex_ratio] [DROPS='{"name":["v1","v2"]}'] Rscript 88g_sens_suite.R
#  입력 = 88g_sens_prep.py 산출 S41_panel_*. 산출 = results_<VER>/S41_part_<analysis>.csv(공변량별 IRR + 적합 요약), S41_hot_<analysis>.csv(핫스팟)
#  GROUP 1 = 우도·시간·평활·Spatial+·지연·극단·이동평균·비교결과 / 2 = 시도 빼기 / 3 = 연도 빼기·짝홀·2024 조건부 검증 / 4 = 대체 잡음 20회
#        5 = 같은 기간 비지연 비교 / 6 = 이원 분해(강건 집합)·변수 제외(DROPS)
suppressPackageStartupMessages({library(INLA);library(sf);library(spdep);library(jsonlite)})
TAG<-Sys.getenv("TAG");stopifnot(nzchar(TAG));VER<-Sys.getenv("ST_VER",unset="std");SUF<-if(nzchar(VER))paste0("_",VER) else ""
ROOT<-normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE))),".."))
SHP<-"${HOME}/Library/CloudStorage/[REDACTED_EMAIL]/내 드라이브/S.K/기타/R_Studio/FBD_DATA_ZIP/final.shp"
OUTD<-file.path(ROOT,"04_분석",TAG);RES<-file.path(OUTD,paste0("results",SUF))
inla.setOption(num.threads="1:1")
LOG<-function(...)cat(sprintf("[88g %s %s] ",TAG,format(Sys.time(),"%H:%M:%S")),sprintf(...),"\n",sep="")
S<-fromJSON(file.path(OUTD,sprintf("config_%s%s.json",TAG,SUF)),simplifyVector=TRUE)
INFO<-fromJSON(file.path(OUTD,sprintf("panel_INFO_%s%s.json",TAG,SUF)))
DD<-read.csv(file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);DD$fbd_region<-enc2utf8(DD$fbd_region)
HB<-readLines(file.path(OUTD,sprintf("HM_best_%s%s.txt",TAG,SUF)));kv<-function(k)sub(paste0("^",k,"="),"",grep(paste0("^",k,"="),HB,value=TRUE)[1])
MAINFAM<-kv("main_family");BEST<-kv("best_set");FORCED<-S$forced;SETS<-lapply(S$sets_after_vif,as.character);HYPV<-SETS[[BEST]];MAINV<-unique(c(HYPV,FORCED))
ROB<-read.csv(file.path(RES,"S7_multiverse_summary.csv"),stringsAsFactors=FALSE);stopifnot(is.logical(ROB$robust),!anyNA(ROB$robust));ROBV<-ROB$cov[ROB$robust];stopifnot(!anyNA(ROBV),all(ROBV%in%MAINV));NR<-length(ROBV)
M3<-read.csv(file.path(RES,"S3_main_IRR.csv"),stringsAsFactors=FALSE);MF<-read.csv(file.path(RES,"refit_main_fixed_unrounded.csv"));msign<-setNames(sign(MF$beta),MF$cov)
FOCAL<-Sys.getenv("FOCAL",unset=if(NR)ROBV[1] else HYPV[1])
LOG("family %s, set %s (%d hyp + %d forced), robust %d (%s), focal %s",MAINFAM,BEST,length(HYPV),length(FORCED),NR,paste(ROBV,collapse=";"),FOCAL)

G<-st_read(SHP,quiet=TRUE);G$region<-enc2utf8(G$region);n<-nrow(G)
g<-inla.read.graph(file.path(RES,sprintf("nb_%s.graph",TAG)));stopifnot(g$n==n)
DD<-DD[order(DD$year,match(DD$fbd_region,G$region)),];DD$area<-match(DD$fbd_region,G$region);Y0<-min(DD$year);DD$time<-DD$year-Y0+1;DD$yiid<-DD$time;DD$id<-seq_len(nrow(DD))
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
SPf<-function(gr="g")sprintf("f(area,model='bym2',graph=%s,scale.model=TRUE,constr=TRUE,hyper=pc)",gr)
TI<-"f(time,n=nt,values=seq_len(nt),model='rw1',scale.model=TRUE,hyper=pcp)";IT<-"f(yiid,model='iid',hyper=pcp)"
bt<-function(v)sprintf("`%s`",v)
irr<-function(m,v){r<-m$summary.fixed[v,,drop=FALSE];pg<-1-inla.pmarginal(0,m$marginals.fixed[[v]])
  data.frame(cov=v,IRR=exp(r$mean),lo=exp(r$`0.025quant`),hi=exp(r$`0.975quant`),P_gt1=pg,beta=r$mean,cred=(r$`0.025quant`>0|r$`0.975quant`<0),stringsAsFactors=FALSE)}
M6<-read.csv(file.path(RES,"S6_spatial_RR.csv"),stringsAsFactors=FALSE);hot0<-M6$fbd_region[M6$class=="hotspot"]
rdp<-function(f){x<-read.csv(file.path(RES,f),fileEncoding="UTF-8",stringsAsFactors=FALSE);x$fbd_region<-enc2utf8(x$fbd_region);x}
mrg<-function(d,X,cols){d<-d[,setdiff(names(d),cols)];d<-merge(d,X[,c("year","fbd_region",cols)],by=c("year","fbd_region"),all.x=TRUE,sort=FALSE);d[order(d$id),]}
fitX<-function(vars,data,fam=MAINFAM,ti=TI,extra=character(0),sp=SPf(),cfg=FALSE){
  nt<-max(data$time);stopifnot(nt>=2,all(data$time>=1))
  rhs<-c("1",bt(vars),extra,sp,ti,IT);rhs<-rhs[nzchar(rhs)];form<-as.formula(paste("cases ~",paste(rhs,collapse="+")))
  strategy<-"default"
  m<-tryCatch(inla(form,family=fam,data=data,E=E,control.fixed=CF,control.compute=list(waic=TRUE,config=cfg),control.predictor=list(compute=TRUE,link=1),
    control.inla=list(control.vb=list(emergency=50)),verbose=FALSE),error=function(e){message("fail: ",conditionMessage(e));NULL})
  if(is.null(m)){strategy<-"eb_gaussian_retry";m<-tryCatch(inla(form,family=fam,data=data,E=E,control.fixed=CF,control.compute=list(waic=TRUE,config=cfg),control.predictor=list(compute=TRUE,link=1),
    control.inla=list(control.vb=list(emergency=50),int.strategy="eb",strategy="gaussian"),verbose=FALSE),error=function(e){message("robust fail: ",conditionMessage(e));NULL})}
  if(!is.null(m))attr(m,"numerical_strategy")<-strategy
  m}
summ<-function(an,item,m,data,note="",vars=MAINV){
  if(is.null(m)){writeLines(toJSON(list(analysis=an,success=FALSE,reason="Both default and EB/Gaussian attempts failed; see group log"),auto_unbox=TRUE),file.path(RES,paste0("S41_status_",an,".json")));LOG("%s FAILED",an);return(invisible())}
  diag<-list(analysis=an,success=TRUE,strategy=attr(m,"numerical_strategy"),n_obs=nrow(data),
    mode_status=m$mode$mode.status,misc_mode_status=m$misc$mode.status,ok=m$ok,internal_warnings=m$misc$warnings,
    fixed_finite=all(is.finite(as.matrix(m$summary.fixed))),hyper_finite=all(is.finite(as.matrix(m$summary.hyperpar))),
    waic=m$waic$waic,waic_finite=is.finite(m$waic$waic))
  writeLines(toJSON(diag,auto_unbox=TRUE,null="null"),file.path(RES,paste0("S41_status_",an,".json")))
  r<-do.call(rbind,lapply(intersect(vars,rownames(m$summary.fixed)),function(v)irr(m,v)))
  kept<-sum(sapply(ROBV,function(v){x<-r[r$cov==v,];nrow(x)==1&&x$cred&&sign(log(x$IRR))==msign[v]}))
  hj<-NA;nh<-NA
  if(!is.null(m$summary.random$area)){ua<-sort(unique(data$area));p1<-sapply(seq_along(ua),function(k)1-inla.pmarginal(0,m$marginals.random$area[[ua[k]]]))
    hot<-G$region[ua][p1>0.95];nh<-length(hot);h0<-intersect(hot0,G$region[ua]);hj<-round(length(intersect(hot,h0))/length(union(hot,h0)),3)
    write.csv(data.frame(fbd_region=G$region[ua],RR=round(exp(m$summary.random$area$mean[ua]),4),P_gt1=round(p1,4)),file.path(RES,sprintf("S41_hot_%s.csv",an)),row.names=FALSE,fileEncoding="UTF-8")}
  w<-r[r$cov==FOCAL,]
  out<-data.frame(analysis=an,item=item,r,WAIC=round(m$waic$waic,1),n_obs=nrow(data),years=sprintf("%d-%d",min(data$year),max(data$year)),robust_kept=kept,n_robust=NR,n_hot=nh,hot_jaccard_vs_main=hj,note=note,stringsAsFactors=FALSE)
  write.csv(out,file.path(RES,sprintf("S41_part_%s.csv",an)),row.names=FALSE,fileEncoding="UTF-8")
  LOG("%-22s kept %d/%d hotJ %s %s %s",an,kept,NR,hj,FOCAL,if(nrow(w))sprintf("%.3f (%.3f-%.3f)%s",w$IRR,w$lo,w$hi,if(w$cred)"*" else "") else "-")}
dn<-function(an)file.exists(file.path(RES,sprintf("S41_part_%s.csv",an)))
reidx<-function(d){d$time<-d$year-min(d$year)+1L;d$yiid<-d$time;d}
SX<-rdp("S41_panel_sido.csv");DD<-mrg(DD,SX,c("sido","sido_year"))
GROUP<-as.integer(Sys.getenv("GROUP","1"))
if(GROUP==1){
  if(!dn("C1_poisson_lognormal")){summ("C1_poisson_lognormal","likelihood",fitX(MAINV,DD,fam="poisson",extra="f(id,model='iid',hyper=pcp)"),DD,"Poisson + observation-level iid")}
  if(!dn("C2_rw2_year")){summ("C2_rw2_year","time",fitX(MAINV,DD,ti="f(time,n=nt,values=seq_len(nt),model='rw2',scale.model=TRUE,hyper=pcp)"),DD,"RW2 year effect")}
  if(!dn("C3_ar1_year")){summ("C3_ar1_year","time",fitX(MAINV,DD,ti="f(time,n=nt,values=seq_len(nt),model='ar1',hyper=list(prec=list(prior='pc.prec',param=c(0.5,0.01))))"),DD,"AR1 year effect")}
  if(!dn("C4_sido_year_iid")){summ("C4_sido_year_iid","time",fitX(MAINV,DD,extra="f(sido_year,model='iid',hyper=pcp)"),DD,"+ province-by-year iid shocks")}
  SM<-intersect(strsplit(Sys.getenv("SMOOTH",unset="sex_ratio"),",")[[1]],MAINV)
  if(length(SM)&&!dn("B3_smooth_forced")){d<-DD;ex<-character(0);for(v in SM){gv<-paste0("g_",v);d[[gv]]<-inla.group(d[[v]],n=20);ex<-c(ex,sprintf("f(%s,model='rw2',scale.model=TRUE,hyper=pcp)",gv))}
    summ("B3_smooth_forced","standardization",fitX(setdiff(MAINV,SM),d,extra=ex),d,paste("RW2 smooth:",paste(SM,collapse="+")),vars=setdiff(MAINV,SM))}
  for(nm in c("splus_cc2_K23","lag1","lag2","extreme","ma3")){an<-c(splus_cc2_K23="D1_spatial_plus_cc2_K23",lag1="E1_lag1",lag2="E2_lag2",extreme="F3_extreme_R7",ma3="G1_kchs_ma3")[[nm]]
    if(nm=="splus_cc2_K23"){
      md<-jsonlite::fromJSON(file.path(RES,"S41_splus_cc2_K23_definition.json"))
      stopifnot(BEST==md$main_set,setequal(MAINV,md$covariates),
        digest::digest(file=file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),algo="sha256")==md$panel_sha256,
        digest::digest(file=file.path(RES,sprintf("nb_%s.graph",TAG)),algo="sha256")==md$graph_sha256,
        digest::digest(file=file.path(RES,"S41_panel_splus_cc2_K23.csv"),algo="sha256")==md$output_sha256)
    }
    if(dn(an))next;X<-rdp(sprintf("S41_panel_%s.csv",nm));cols<-intersect(names(X)[-(1:2)],MAINV);d<-mrg(DD,X,cols);d<-d[complete.cases(d[,cols]),];if(nm%in%c("lag1","lag2"))d<-reidx(d)
    summ(an,"robustness",fitX(MAINV,d),d,nm)}
  for(t2 in strsplit(Sys.getenv("NCO",unset=""),",")[[1]]){if(!nzchar(t2))next;an<-paste0("H1_nco_",t2);if(dn(an))next;X<-rdp(sprintf("S41_panel_nco_%s.csv",t2));d<-DD;k<-match(paste(d$year,d$fbd_region),paste(X$year,X$fbd_region))
    d$cases<-X[[paste0("cases_",t2)]][k];d$E<-X[[paste0("E_",t2)]][k];d<-d[complete.cases(d[,c("cases","E")]),];summ(an,"comparison_outcome",fitX(MAINV,d),d,paste("outcome =",t2))}
}
if(GROUP==2){for(s in sort(unique(DD$sido))){an<-paste0("F1_drop_sido_",s);if(dn(an))next;d<-DD[DD$sido!=s,];summ(an,"leave_one_province",fitX(MAINV,d),d,paste("without province",s))}}
if(GROUP==3){
  for(y in sort(unique(DD$year))){an<-paste0("F2_drop_year_",y);if(dn(an))next;d<-reidx(DD[DD$year!=y,]);summ(an,"leave_one_year",fitX(MAINV,d),d,paste("without",y))}
  for(pp in c("odd","even")){an<-paste0("H4_",pp,"_years");if(dn(an))next;d<-reidx(DD[(DD$year%%2==1)==(pp=="odd"),]);summ(an,"split_half",fitX(MAINV,d),d,paste(pp,"years"))}
  if(!file.exists(file.path(RES,"S41_oos2024.csv"))){d<-DD;obs<-d$cases;d$cases[d$year==2024]<-NA;m<-fitX(MAINV,d,cfg=TRUE)
    ds<-if(is.null(m))list(analysis="conditional_2024",success=FALSE,reason="Fit returned NULL; see group3 log") else list(analysis="conditional_2024",success=TRUE,strategy=attr(m,"numerical_strategy"),n_rows=nrow(d),n_observed=sum(!is.na(d$cases)),n_heldout=sum(is.na(d$cases)),ok=m$ok,mode_status=m$mode$mode.status,misc_mode_status=m$misc$mode.status,internal_warnings=m$misc$warnings,fixed_finite=all(is.finite(as.matrix(m$summary.fixed))),hyper_finite=all(is.finite(as.matrix(m$summary.hyperpar))),waic_finite=is.finite(m$waic$waic),offset_note="held-out-year E and covariates retained; conditional check, not prospective forecast")
    writeLines(toJSON(ds,auto_unbox=TRUE,null="null"),file.path(RES,"S41_status_conditional_2024.json"))
    if(is.null(m)||!isTRUE(m$ok)||!isTRUE(m$mode$mode.status==0)||!isTRUE(ds$fixed_finite)||!isTRUE(ds$hyper_finite))stop("Conditional2024 fit failed numerical checks; no posterior sampling")
    set.seed(202611042L)
    SS<-inla.posterior.sample(200,m,num.threads="1:1");k<-which(DD$year==2024);hn<-names(SS[[1]]$hyperpar);szn<-grep("size",hn,value=TRUE)[1]
    Y<-sapply(SS,function(s){eta<-s$latent[grep("^Predictor:",rownames(s$latent)),1][k];if(is.na(szn))rpois(length(k),DD$E[k]*exp(eta)) else rnbinom(length(k),size=s$hyperpar[szn],mu=DD$E[k]*exp(eta))})
    lo<-apply(Y,1,quantile,.025);hi<-apply(Y,1,quantile,.975);med<-apply(Y,1,median)
    res<-data.frame(n=length(k),coverage95=round(mean(obs[k]>=lo&obs[k]<=hi),4),median_abs_pct_error=round(median(abs(med-obs[k])/obs[k])*100,2),spearman_obs_pred=round(cor(obs[k],med,method="spearman"),4))
    write.csv(res,file.path(RES,"S41_oos2024.csv"),row.names=FALSE);LOG("2024 conditional check: cov95 %.3f MdAPE %.2f%% rho %.3f",res$coverage95,res$median_abs_pct_error,res$spearman_obs_pred)}
}
if(GROUP==4){for(r_ in 1:20){an<-sprintf("G2_noise_r%02d",r_);if(dn(an))next;f<-sprintf("S41_panel_noise_r%02d.csv",r_);if(!file.exists(file.path(RES,f))){LOG("no imputed cells, skip noise");break}
  X<-rdp(f);cols<-intersect(names(X)[-(1:2)],MAINV);d<-mrg(DD,X,cols);summ(an,"imputation_noise",fitX(MAINV,d),d,"noise on imputed cells")}}
if(GROUP==5){for(y0 in c(2012,2013)){an<-sprintf("E0_nolag_from%d",y0);if(dn(an))next;d<-reidx(DD[DD$year>=y0,]);summ(an,"lag_comparator",fitX(MAINV,d),d,sprintf("unlagged, same years as lag analysis (from %d)",y0))}}
if(GROUP==6){
  if(NR&&!dn("bw_twoway")){d<-DD;BW<-character(0)
    for(v in ROBV){mi<-ave(d[[v]],d$area,FUN=mean);mt<-ave(d[[v]],d$year,FUN=mean);d[[paste0(v,"__between")]]<-mi;d[[paste0(v,"__within")]]<-d[[v]]-mi-mt+mean(d[[v]]);BW<-c(BW,paste0(v,c("__between","__within")))}
    summ("bw_twoway","between_within",fitX(c(setdiff(MAINV,ROBV),BW),d),d,"robust set: district mean + two-way centred deviation",vars=c(setdiff(MAINV,ROBV),BW))}
  DR<-Sys.getenv("DROPS",unset="");if(nzchar(DR)){DL<-fromJSON(DR);for(nm in names(DL)){an<-paste0("drop_",nm);if(dn(an))next;v<-setdiff(MAINV,DL[[nm]]);summ(an,"drop_covariates",fitX(v,DD),DD,paste("without",paste(DL[[nm]],collapse="+")),vars=v)}}
}
LOG("GROUP %d done",GROUP)
