# 81 — 질환 시공간 결정요인 INLA 적합 (시군구×연도 패널). 사용: Rscript 81_fit_inla.R <TAG>   [ST_VER=v2 → config/panel *_v2 입력, results_v2/, HM_*_v2]
#  2026-10-19 반영: 미관측·연도별 결측 시군구 Moran NA 처리(moran_obs, 관측 지역 부분그래프) — _PATCH_NOTES.md
#  모형: cases ~ offset(log E) + forced + set + f(area, bym2, PC) + f(time, rw1) + f(yiid, iid)  [+ f(id, iid) Type I 민감도]
#  단계(각 CSV 저장, 있으면 건너뜀): S1 우도 선택 → S2 세트 비교·선정(rules) → S3 주모형·음성대조·공간 RR·잔차 Moran
#          → S4 창별 재적합 + 창×변수 상호작용 → S5 COVID 전후 → S6 연도별 Moran → S7 다중우주 → S9 Type I 민감도
suppressPackageStartupMessages({library(INLA);library(sf);library(spdep);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE);TAG<-args[1];stopifnot(!is.na(TAG))
ROOT<-normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE))),".."))
SHP<-"${HOME}/Library/CloudStorage/[REDACTED_EMAIL]/내 드라이브/S.K/기타/R_Studio/FBD_DATA_ZIP/final.shp"
VER<-Sys.getenv("ST_VER",unset="");SUF<-if(nzchar(VER))paste0("_",VER) else ""   # 버전 분기
OUTD<-file.path(ROOT,"04_분석",TAG);RES<-Sys.getenv("AUDIT_OUT");dir.create(RES,recursive=TRUE,showWarnings=FALSE)
NT<-as.integer(Sys.getenv("ST_THREADS",unset="1"));inla.setOption(num.threads=sprintf("%d:1",NT))   # 2026-10-19 실측: 2스레드면 같은 모형 WAIC가 ±50 요동(phi 0.84~0.99) → 1스레드(±4)로 고정, 세트 비교 재현성 확보
S<-fromJSON(file.path(OUTD,sprintf("config_%s%s.json",TAG,SUF)),simplifyVector=TRUE);RU<-S$rules   # 80이 정규화해 저장한 설정(sets json 원본은 _st_config가 해석)
INFO<-fromJSON(file.path(OUTD,sprintf("panel_INFO_%s%s.json",TAG,SUF)))
META<-read.csv(file.path(OUTD,sprintf("panel_META_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE)
LAB<-setNames(META$label_en,META$variable);KOLAB<-setNames(META$label_ko,META$variable)
DD<-read.csv(file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE)
if(nzchar(Sys.getenv("AUDIT_PANEL")))DD<-read.csv(Sys.getenv("AUDIT_PANEL"),stringsAsFactors=FALSE)
DD$fbd_region<-enc2utf8(DD$fbd_region)
FORCED<-S$forced;NEG<-if(is.null(S$negative_controls))character(0) else S$negative_controls
SETS<-S$sets_after_vif;SETS<-lapply(SETS,function(v)as.character(v));SENS<-if(is.null(S$sensitivity_sets))character(0) else S$sensitivity_sets;SENSV<-if(is.null(S$sensitivity_only_vars))character(0) else S$sensitivity_only_vars
CAND<-setdiff(unique(unlist(SETS[!names(SETS)%in%SENS])),c(FORCED,SENSV))   # 다중우주·선정 후보 풀 = 사전지정 세트의 변수(민감도 전용 변수 제외)
stopifnot(all(c(FORCED,CAND,SENSV,NEG)%in%names(DD)))
WTOL<-RU$waic_tol;CPF<-RU$cpo_fail_frac
NMV<-RU$multiverse_n;if(nzchar(Sys.getenv("ST_MV_N")))NMV<-as.integer(Sys.getenv("ST_MV_N"))
WIN<-INFO$windows;POSTY<-RU$covid_post_from   # 창은 80이 관측 연도에 맞춰 자른 것(panel_INFO)
LOG<-function(...)cat(sprintf("[81 %s] ",format(Sys.time(),"%H:%M:%S")),sprintf(...),"\n",sep="")

## ── 이웃행렬: queen 인접, 고립 폴리곤은 최근접 1개 연결(대칭화) ──
G<-st_read(SHP,quiet=TRUE);G$region<-enc2utf8(G$region)
stopifnot(all(unique(DD$fbd_region)%in%G$region));OBSR<-G$region%in%unique(DD$fbd_region)   # 미관측 시군구(outcome 전 연도 결측)는 그래프에 남기고 모형에서는 자료 없음
GF<-file.path(RES,sprintf("nb_%s.graph",TAG));NBF<-file.path(RES,sprintf("nb_%s_summary.csv",TAG))
nb<-poly2nb(G,queen=TRUE);isl<-which(card(nb)==0)
if(length(isl)){cen<-st_coordinates(st_centroid(st_geometry(G)));kn<-knearneigh(cen,k=1)$nn
  for(i in isl){j<-kn[i,1];nb[[i]]<-sort(unique(c(nb[[i]][nb[[i]]>0],j)));nb[[j]]<-sort(unique(c(nb[[j]][nb[[j]]>0],i)))}}
stopifnot(all(card(nb)>0));nb2INLA(GF,nb);g<-inla.read.graph(GF);lw<-nb2listw(nb,style="W",zero.policy=TRUE)
nbo<-subset(nb,OBSR);lwo<-nb2listw(nbo,style="W",zero.policy=TRUE);GO<-G$region[OBSR]   # Moran용: 관측 시군구만
moran_obs<-function(x,regions){k<-!is.na(x);if(sum(k)<10)return(list(statistic=NA,p.value=NA));kk<-G$region%in%regions[k];nbk<-subset(nb,kk);lwk<-nb2listw(nbk,style="W",zero.policy=TRUE);moran.mc(x[k][match(G$region[kk],regions[k])],lwk,999,zero.policy=TRUE)}   # METS 패치: NA 지역 제외 부분그래프 Moran
write.csv(data.frame(n_regions=length(nb),n_links=sum(card(nb)),n_subgraphs=n.comp.nb(nb)$nc,islands_linked=length(isl),island_regions=paste(G$region[isl],collapse=";"),mean_neighbors=round(mean(card(nb)),2)),NBF,row.names=FALSE)
LOG("graph: %d regions, %d islands linked (%s)",length(nb),length(isl),paste(G$region[isl],collapse=","))

## ── 데이터 색인 ──
DD$area<-match(DD$fbd_region,G$region);Y0<-min(DD$year);DD$time<-DD$year-Y0+1;DD$yiid<-DD$time;DD$id<-seq_len(nrow(DD))
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
SP<-"f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)";TI<-"f(time,model='rw1',scale.model=TRUE,hyper=pcp)";IT<-"f(yiid,model='iid',hyper=pcp)";ST1<-"f(id,model='iid',hyper=pcp)"
bt<-function(v)sprintf("`%s`",v)
fitf<-function(vars,fam,data=DD,cpo=TRUE,typeI=FALSE,extra=character(0),re=TRUE,ecol="E",robust=FALSE){
  rhs<-c("1",bt(vars),extra,if(re)c(SP,TI,IT),if(typeI)ST1);form<-as.formula(paste("cases ~",paste(rhs,collapse="+")))
  ci<-if(robust)list(control.vb=list(emergency=50),int.strategy="eb",strategy="gaussian") else list(control.vb=list(emergency=50))
  inla(form,family=fam,data=data,E=data[[ecol]],control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE,cpo=cpo,config=FALSE),control.inla=ci,verbose=FALSE)}
safefit<-function(...)tryCatch(fitf(...),error=function(e){message("  ✗ fit failed: ",conditionMessage(e)," → retry robust(eb/gaussian)");tryCatch(fitf(...,robust=TRUE),error=function(e2){message("  ✗ retry failed: ",conditionMessage(e2));NULL})})
irr<-function(m,v){r<-m$summary.fixed[v,,drop=FALSE];data.frame(cov=v,IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),
  P_gt1=round(1-INLA::inla.pmarginal(0,m$marginals.fixed[[v]]),3),cred=(r$`0.025quant`>0|r$`0.975quant`<0),stringsAsFactors=FALSE)}
fitstat<-function(m,fam){cp<-m$cpo$cpo;data.frame(family=fam,WAIC=round(m$waic$waic,1),DIC=round(m$dic$dic,1),LCPO=round(-mean(log(cp[is.finite(cp)&cp>0])),4),cpo_fail=sum(m$cpo$failure>0,na.rm=TRUE),
  phi=round(if("Phi for area"%in%rownames(m$summary.hyperpar))m$summary.hyperpar["Phi for area","mean"] else NA,3),
  size=round(if(any(grepl("size",rownames(m$summary.hyperpar))))m$summary.hyperpar[grep("size",rownames(m$summary.hyperpar))[1],"mean"] else NA,3))}
FAMS<-c("poisson","nbinomial")
timed<-function(label,expr){t0<-Sys.time();r<-expr;LOG("%s · %.1f min",label,as.numeric(difftime(Sys.time(),t0,units="mins")));r}


if(nzchar(Sys.getenv("AUDIT_GRAPH")))g<-inla.read.graph(Sys.getenv("AUDIT_GRAPH"))
MAINV<-unique(c(SETS[[RU$main_set]],FORCED))
mode<-Sys.getenv("AUDIT_MODE",unset="main")
if(mode=="post2014")DD<-DD[DD$year>=2014,]
if(mode=="khIV"){
 DD$space_st<-DD$area;DD$time_st<-DD$time
 form_extra<-"f(space_st,model='besag',graph=g,scale.model=TRUE,constr=TRUE,hyper=pcp,group=time_st,control.group=list(model='rw1',scale.model=TRUE))"
}else form_extra<-character(0)
set.seed(20261104)
m<-fitf(MAINV,"nbinomial",extra=form_extra)
saveRDS(m,file.path(RES,"model.rds"))
write.csv(do.call(rbind,lapply(MAINV,function(v)irr(m,v))),file.path(RES,"IRR.csv"),row.names=FALSE)
write.csv(fitstat(m,"nbinomial"),file.path(RES,"fit.csv"),row.names=FALSE)
write.csv(data.frame(year=DD$year,fbd_region=DD$fbd_region,waic=m$waic$local.waic),file.path(RES,"pointwise_waic.csv"),row.names=FALSE)
write.csv(m$summary.hyperpar,file.path(RES,"hyper.csv"))
LOG("AUDIT DONE %s",mode)
