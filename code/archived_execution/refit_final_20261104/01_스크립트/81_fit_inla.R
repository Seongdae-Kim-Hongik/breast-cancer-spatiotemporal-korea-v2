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
OUTD<-file.path(ROOT,"04_분석",TAG);RES<-file.path(OUTD,paste0("results",SUF));dir.create(RES,showWarnings=FALSE)
NT<-as.integer(Sys.getenv("ST_THREADS",unset="1"));inla.setOption(num.threads=sprintf("%d:1",NT))   # 2026-10-19 실측: 2스레드면 같은 모형 WAIC가 ±50 요동(phi 0.84~0.99) → 1스레드(±4)로 고정, 세트 비교 재현성 확보
S<-fromJSON(file.path(OUTD,sprintf("config_%s%s.json",TAG,SUF)),simplifyVector=TRUE);RU<-S$rules   # 80이 정규화해 저장한 설정(sets json 원본은 _st_config가 해석)
INFO<-fromJSON(file.path(OUTD,sprintf("panel_INFO_%s%s.json",TAG,SUF)))
META<-read.csv(file.path(OUTD,sprintf("panel_META_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE)
LAB<-setNames(META$label_en,META$variable);KOLAB<-setNames(META$label_ko,META$variable)
DD<-read.csv(file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE)
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
if(nzchar(Sys.getenv("BRCA_GRAPH"))){
 supplied<-inla.read.graph(Sys.getenv("BRCA_GRAPH"));stopifnot(supplied$n==nrow(G))
 for(i in seq_along(nb))nb[[i]]<-as.integer(sort(supplied$nbs[[i]]))
 attr(nb,"sym")<-TRUE
}
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

## ── S1 우도 선택(강제통제만) ──
F1<-file.path(RES,"S1_likelihood.csv")
if(!file.exists(F1)){lk<-NULL;for(fam in FAMS){m<-timed(sprintf("S1 %s",fam),safefit(FORCED,fam));if(is.null(m))next;lk<-rbind(lk,cbind(set="REF_forced_only",nvar=length(FORCED),n_hyp=0,fitstat(m,fam)))}
  write.csv(lk,F1,row.names=FALSE)}
lk<-read.csv(F1);MAINFAM<-lk$family[which.min(lk$WAIC)];LOG("S1 main family = %s (%s)",MAINFAM,paste(sprintf("%s %.1f",lk$family,lk$WAIC),collapse=" vs "))

## ── S2 세트 비교(주 우도) ──
CF1<-file.path(OUTD,sprintf("HM_setcompare_%s%s.csv",TAG,SUF));CF2<-file.path(OUTD,sprintf("HM_setcompare_IRR_%s%s.csv",TAG,SUF))
cmp<-if(file.exists(CF1))read.csv(CF1,stringsAsFactors=FALSE) else NULL;est<-if(file.exists(CF2))read.csv(CF2,stringsAsFactors=FALSE) else NULL
if(is.null(cmp)||!("REF_forced_only"%in%cmp$set)){rf<-lk[lk$family==MAINFAM,];cmp<-rbind(cmp,data.frame(set="REF_forced_only",family=MAINFAM,nvar=length(FORCED),n_hyp=0,WAIC=rf$WAIC,DIC=rf$DIC,LCPO=rf$LCPO,cpo_fail=rf$cpo_fail,minutes=NA,sensitivity=FALSE))}
if(!"sensitivity"%in%names(cmp))cmp$sensitivity<-cmp$set%in%SENS
for(sn in names(SETS)){if(any(cmp$set==sn))next;cv<-unique(c(SETS[[sn]],FORCED));t0<-Sys.time();m<-safefit(cv,MAINFAM);if(is.null(m))next;fs<-fitstat(m,MAINFAM)
  cmp<-rbind(cmp,data.frame(set=sn,family=MAINFAM,nvar=length(cv),n_hyp=length(SETS[[sn]]),WAIC=fs$WAIC,DIC=fs$DIC,LCPO=fs$LCPO,cpo_fail=fs$cpo_fail,minutes=round(as.numeric(difftime(Sys.time(),t0,units="mins")),1),sensitivity=sn%in%SENS))
  for(v in cv)est<-rbind(est,data.frame(set=sn,family=MAINFAM,irr(m,v),forced=v%in%FORCED))
  LOG("S2 %-26s nvar %2d WAIC %9.1f LCPO %.4f cpo_fail %d · %.1f min",sn,length(cv),fs$WAIC,fs$LCPO,fs$cpo_fail,as.numeric(difftime(Sys.time(),t0,units="mins")))
  write.csv(cmp,CF1,row.names=FALSE);write.csv(est,CF2,row.names=FALSE)}
# 선정 규칙(sets json rules): 주 우도에서 WAIC 최소; 허용폭 내면 가설 변수 수 최소; 동률이면 LCPO(CPO 실패 <cpo_fail_frac일 때만)
mc<-cmp[!grepl("^REF_",cmp$set)&!cmp$sensitivity,];mc<-mc[order(mc$WAIC),];best<-mc$set[1]   # 민감도 세트(S1b)는 선정 후보에서 제외
near<-mc[mc$WAIC-mc$WAIC[1]<=WTOL,];if(nrow(near)>1){near<-near[near$n_hyp==min(near$n_hyp),];near<-if(all(near$cpo_fail<CPF*nrow(DD)))near[order(near$LCPO,near$WAIC),] else near[order(near$WAIC),];best<-near$set[1]}
RULEBEST<-best;MSET<-RU$main_set;if(!is.null(MSET)&&length(MSET)==1&&!is.na(MSET)&&nzchar(MSET)){stopifnot(MSET%in%cmp$set);best<-MSET}
writeLines(c(sprintf("tag=%s",TAG),sprintf("main_family=%s",MAINFAM),sprintf("best_set=%s",best),sprintf("rule_best_set=%s",RULEBEST),
  sprintf("main_set_mode=%s",if(best==RULEBEST&&(is.null(MSET)||length(MSET)!=1||is.na(MSET)))"rule" else "json rules main_set (author decision); rule-selected set kept as sensitivity"),
  sprintf("family_rule=lower WAIC of REF_forced_only (%s)",paste(sprintf("%s %.1f",lk$family,lk$WAIC),collapse=" vs ")),
  sprintf("set_rule=min WAIC in main family; within %g WAIC choose fewest hypothesis variables; ties by lower LCPO (only if CPO failures <%g%% of cells) [rules from sets json]",WTOL,100*CPF),
  sprintf("ranking=%s",paste(sprintf("%s(%.1f, LCPO %.4f, n_hyp %d)",mc$set,mc$WAIC,mc$LCPO,mc$n_hyp),collapse=" > ")),
  sprintf("ref_forced_only_WAIC=%.1f",cmp$WAIC[cmp$set=="REF_forced_only"]),sprintf("sensitivity_sets_not_eligible=%s",paste(SENS,collapse=";"))),file.path(OUTD,sprintf("HM_best_%s%s.txt",TAG,SUF)))
LOG("S2 best set = %s (rule %s)",best,RULEBEST)
MAINV<-unique(c(SETS[[best]],FORCED));HYPV<-SETS[[best]]

## ── S3 주모형 + 음성대조 + 공간 RR + 잔차 Moran ──
F3<-file.path(RES,"S3_main_IRR.csv");F3f<-file.path(RES,"S3_main_fit.csv");F3n<-file.path(RES,"S3_negcontrol.csv");F6<-file.path(RES,"S6_spatial_RR.csv");F8<-file.path(RES,"S8_residual_moran.csv");F3h<-file.path(RES,"S3_main_hyper.csv");F3y<-file.path(RES,"S3_year_effect.csv")
if(!all(file.exists(F3,F3f,F3n,F6,F8))){
  m<-timed("S3 main",safefit(MAINV,MAINFAM));stopifnot(!is.null(m));saveRDS(m,file.path(RES,"refit_main.rds"));write.csv(data.frame(year=DD$year,fbd_region=DD$fbd_region,waic=m$waic$local.waic),file.path(RES,"refit_main_pointwise_waic.csv"),row.names=FALSE)
  ir<-do.call(rbind,lapply(MAINV,function(v)irr(m,v)));ir$label_en<-LAB[ir$cov];ir$label_ko<-KOLAB[ir$cov];ir$role<-ifelse(ir$cov%in%FORCED,"forced","hypothesis");ir$set<-best;ir$family<-MAINFAM
  write.csv(ir,F3,row.names=FALSE)
  hp<-m$summary.hyperpar;write.csv(data.frame(hyperpar=rownames(hp),round(hp[,c("mean","0.025quant","0.975quant")],4)),F3h,row.names=FALSE)
  # 공간 RR(bym2 총효과 exp) · 초과확률
  n<-length(nb);ur<-m$summary.random$area[1:n,];pep<-sapply(m$marginals.random$area[1:n],function(mg)1-inla.pmarginal(0,mg))
  rr<-data.frame(area=1:n,fbd_region=G$region,sigungu_cd=G$sigungu_cd,RR=round(exp(ur$mean),3),lo=round(exp(ur$`0.025quant`),3),hi=round(exp(ur$`0.975quant`),3),P_gt1=round(pep,3))
  rr$class<-ifelse(rr$P_gt1>0.95,"hotspot",ifelse(rr$P_gt1<0.05,"coldspot","neutral"));rr$observed<-OBSR;rr$class[!OBSR]<-"unobserved"
  agg<-aggregate(cbind(cases,E)~fbd_region,DD,sum);rr$SIR_obs<-round((agg$cases/agg$E)[match(rr$fbd_region,agg$fbd_region)],3)
  write.csv(rr,F6,row.names=FALSE)
  ye<-m$summary.random$time;write.csv(data.frame(year=Y0+ye$ID-1,rw1_mean=round(ye$mean,4),lo=round(ye$`0.025quant`,4),hi=round(ye$`0.975quant`,4)),F3y,row.names=FALSE)
  # 잔차 Moran(적합값 기준 Pearson, 시군구 평균)
  mu<-m$summary.fitted.values$mean*DD$E;sz<-if(MAINFAM=="nbinomial")m$summary.hyperpar[grep("size",rownames(m$summary.hyperpar))[1],"mean"] else Inf
  pr<-(DD$cases-mu)/sqrt(mu+mu^2/sz);pra<-tapply(pr,DD$area,mean);set.seed(1)
  mr<-moran_obs(as.numeric(pra[as.character(which(OBSR))]),GO);mo<-moran_obs(log(agg$cases/agg$E)[match(GO,agg$fbd_region)],GO)
  write.csv(data.frame(quantity=c("observed_logSIR_pooled","pearson_residual_area_mean"),moran_I=round(c(mo$statistic,mr$statistic),4),p_perm=round(c(mo$p.value,mr$p.value),4),nperm=999),F8,row.names=FALSE)
  fs<-fitstat(m,MAINFAM)
  write.csv(data.frame(fs,set=best,nvar=length(MAINV),n_hyp=length(HYPV),n_obs=nrow(DD),n_regions=sum(OBSR),n_years=length(unique(DD$year)),total_cases=sum(DD$cases),
    n_hotspot=sum(rr$class=="hotspot"),n_coldspot=sum(rr$class=="coldspot"),n_unobserved=sum(!OBSR),RR_min=min(rr$RR[OBSR]),RR_max=max(rr$RR[OBSR]),n_cred_hyp=sum(ir$cred[ir$role=="hypothesis"]),
    cred_hyp=paste(ir$cov[ir$cred&ir$role=="hypothesis"],collapse=";"),resid_moran_I=round(mr$statistic,4),resid_moran_p=round(mr$p.value,4),obs_moran_I=round(mo$statistic,4),obs_moran_p=round(mo$p.value,4)),F3f,row.names=FALSE)
  # 음성대조: 주모형 + 음성대조 1개씩
  nc<-NULL;for(v in NEG){mm<-timed(sprintf("S3 negcontrol %s",v),safefit(c(MAINV,v),MAINFAM));if(is.null(mm))next;r<-irr(mm,v);r$label_en<-LAB[v];r$label_ko<-KOLAB[v];r$WAIC<-round(mm$waic$waic,1);nc<-rbind(nc,r)}
  write.csv(if(is.null(nc))data.frame(cov=character(0))else nc,F3n,row.names=FALSE)
  LOG("S3 done: credible hyp = %s; hotspots %d coldspots %d",paste(ir$cov[ir$cred&ir$role=="hypothesis"],collapse=","),sum(rr$class=="hotspot"),sum(rr$class=="coldspot"))
  rm(m)}

## ── S4 창별 재적합 + 창×변수 상호작용(pooled, 겹치지 않는 창) ──
F4<-file.path(RES,"S4_window_IRR.csv");F4i<-file.path(RES,"S4_interaction.csv");F4w<-file.path(RES,"S4_window_fit.csv")
if(!all(file.exists(F4,F4i,F4w))){wr<-NULL;wf<-NULL
  for(wn in names(WIN)){yy<-WIN[[wn]];d<-DD[DD$year>=yy[1]&DD$year<=yy[2],];d$time<-d$year-min(d$year)+1;d$yiid<-d$time;d$id<-seq_len(nrow(d))
    m<-timed(sprintf("S4 window %s",wn),safefit(MAINV,MAINFAM,data=d));if(is.null(m))next
    r<-do.call(rbind,lapply(MAINV,function(v)irr(m,v)));r$window<-wn;r$years<-sprintf("%d-%d",yy[1],yy[2]);r$label_en<-LAB[r$cov];r$label_ko<-KOLAB[r$cov];r$role<-ifelse(r$cov%in%FORCED,"forced","hypothesis");wr<-rbind(wr,r)
    wf<-rbind(wf,data.frame(window=wn,years=sprintf("%d-%d",yy[1],yy[2]),n_obs=nrow(d),n_years=length(unique(d$year)),cases=sum(d$cases),WAIC=round(m$waic$waic,1),n_cred_hyp=sum(r$cred[r$role=="hypothesis"])))}
  write.csv(wr,F4,row.names=FALSE);write.csv(wf,F4w,row.names=FALSE)
  d<-DD;wl<-names(WIN);d$W<-factor(d$window,levels=wl);ex<-character(0)
  for(k in 2:length(wl)){dn<-sprintf("W%d",k);d[[dn]]<-as.integer(d$W==wl[k]);ex<-c(ex,dn);for(v in HYPV){nm<-sprintf("%s_x_W%d",v,k);d[[nm]]<-d[[v]]*d[[dn]];ex<-c(ex,bt(nm))}}
  m<-timed("S4 interaction pooled",safefit(MAINV,MAINFAM,data=d,extra=ex));stopifnot(!is.null(m))
  ir<-NULL;for(k in 2:length(wl))for(v in HYPV){nm<-sprintf("%s_x_W%d",v,k);r<-irr(m,nm);r$cov<-v;r$window<-wl[k];r$ref_window<-wl[1];r$label_en<-LAB[v];r$label_ko<-KOLAB[v];names(r)[names(r)=="IRR"]<-"IRR_ratio";ir<-rbind(ir,r)}
  ir$WAIC<-round(m$waic$waic,1);write.csv(ir,F4i,row.names=FALSE);LOG("S4 done: interaction credible = %d/%d",sum(ir$cred),nrow(ir));rm(m)}

## ── S5 COVID 전후(post 더미 + 변수×post) ──
F5<-file.path(RES,"S5_covid.csv")
if(!file.exists(F5)){d<-DD;ex<-"post";for(v in HYPV){nm<-sprintf("%s_x_post",v);d[[nm]]<-d[[v]]*d$post;ex<-c(ex,bt(nm))}
  m<-timed("S5 covid",safefit(MAINV,MAINFAM,data=d,extra=ex));stopifnot(!is.null(m))
  r<-irr(m,"post");r$term<-"post_main";r$cov<-"post";r$label_en<-sprintf("COVID period (%d+) main effect",POSTY);r$label_ko<-sprintf("COVID 시기(%d~) 주효과",POSTY);out<-r
  for(v in HYPV){r<-irr(m,sprintf("%s_x_post",v));r$term<-"var_x_post";r$cov<-v;r$label_en<-LAB[v];r$label_ko<-KOLAB[v];out<-rbind(out,r)}
  for(v in HYPV){r<-irr(m,v);r$term<-"var_pre";r$cov<-v;r$label_en<-LAB[v];r$label_ko<-KOLAB[v];out<-rbind(out,r)}
  out$WAIC<-round(m$waic$waic,1);out$post_from<-POSTY;write.csv(out,F5,row.names=FALSE);LOG("S5 done");rm(m)}

## ── S6 연도별 Moran's I(관측 log SIR) ──
F6m<-file.path(RES,"S6_moran_year.csv")
if(!file.exists(F6m)){set.seed(2);my<-do.call(rbind,lapply(sort(unique(DD$year)),function(y){d<-DD[DD$year==y,];x<-log((d$cases+0.5)/d$E)[match(GO,d$fbd_region)];mm<-moran_obs(x,GO)
  data.frame(year=y,moran_I=round(mm$statistic,4),p_perm=round(mm$p.value,4),n_regions=nrow(d),crude_rate_per100k=round(sum(d$cases)/sum(d$population)*1e5,1))}));write.csv(my,F6m,row.names=FALSE);LOG("S6 done: Moran max %.3f (%d)",max(my$moran_I),my$year[which.max(my$moran_I)])}

if(Sys.getenv("ST_ONLY_CORE")=="1"){
 main_saved<-readRDS(file.path(RES,"refit_main.rds"));rr<-main_saved$summary.fixed
 write.csv(data.frame(cov=rownames(rr),beta=rr$mean,sd=rr$sd,lo=rr$`0.025quant`,hi=rr$`0.975quant`),file.path(RES,"refit_main_fixed_unrounded.csv"),row.names=FALSE)
 LOG("CORE ONLY DONE — multiverse deferred to isolated workers")
 quit(save="no",status=0)
}

## ── S7 다중우주: forced 고정 + 후보 변수 무작위 부분집합 ──
F7s<-file.path(RES,"S7_multiverse_specs.csv");F7r<-file.path(RES,"S7_multiverse_results.csv");F7<-file.path(RES,"S7_multiverse_summary.csv")
if(!file.exists(F7)){set.seed(20261019);POOL<-CAND;inla.setOption(inla.timeout=600);LOG("S7 inla.timeout=600s per spec (failed/timeout specs are skipped and logged)")
  if(!file.exists(F7s)){specs<-list();keys<-character(0);while(length(specs)<NMV){inc<-POOL[runif(length(POOL))<0.5];if(length(inc)<1)next;k<-paste(sort(inc),collapse="|");if(length(POOL)>=5&&k%in%keys)next;specs[[length(specs)+1]]<-inc;keys<-c(keys,k)}
    write.csv(data.frame(spec=seq_along(specs),nvar=sapply(specs,length),vars=sapply(specs,paste,collapse="|")),F7s,row.names=FALSE)}
  sp<-read.csv(F7s,stringsAsFactors=FALSE);mv<-if(file.exists(F7r))read.csv(F7r,stringsAsFactors=FALSE) else NULL;t0<-Sys.time()
  for(i in sp$spec){if(!is.null(mv)&&i%in%mv$spec)next;vv<-strsplit(sp$vars[sp$spec==i],"\\|")[[1]];m<-safefit(c(vv,FORCED),MAINFAM,cpo=FALSE);if(is.null(m)){LOG("S7 spec %d skipped (fit failed/timeout)",i);next}
    r<-do.call(rbind,lapply(vv,function(v)irr(m,v)));r$spec<-i;r$nvar<-length(vv);mv<-rbind(mv,r);if(i%%10==0){write.csv(mv,F7r,row.names=FALSE);LOG("S7 spec %d/%d · %.1f min elapsed",i,nrow(sp),as.numeric(difftime(Sys.time(),t0,units="mins")))}}
  write.csv(mv,F7r,row.names=FALSE);inla.setOption(inla.timeout=0)
  main<-read.csv(F3,stringsAsFactors=FALSE);sgn<-setNames(sign(log(main$IRR)),main$cov)
  sm<-do.call(rbind,lapply(POOL,function(v){x<-mv[mv$cov==v,];if(!nrow(x))return(NULL);ref<-if(v%in%names(sgn))sgn[v] else sign(median(log(x$IRR)))
    data.frame(cov=v,label_en=LAB[v],label_ko=KOLAB[v],in_main_set=v%in%HYPV,n_specs=nrow(x),cred_frac=round(mean(x$cred),3),sign_agree=round(mean(sign(log(x$IRR))==ref),3),
      IRR_median=round(median(x$IRR),3),IRR_min=round(min(x$IRR),3),IRR_max=round(max(x$IRR),3),robust=(mean(x$cred)>=0.75&mean(sign(log(x$IRR))==ref)>=0.95),provisional=nrow(x)<50)}))
  sm$n_specs_total<-nrow(sp);write.csv(sm,F7,row.names=FALSE);LOG("S7 done: robust = %s",paste(sm$cov[sm$robust],collapse=","))}

## ── S9 Type I 공간×시간 iid 민감도 ──
F9<-file.path(RES,"S9_typeI_IRR.csv")
if(!file.exists(F9)){m<-timed("S9 Type I",safefit(MAINV,MAINFAM,typeI=TRUE));if(!is.null(m)){r<-do.call(rbind,lapply(MAINV,function(v)irr(m,v)));r$label_en<-LAB[r$cov];r$role<-ifelse(r$cov%in%FORCED,"forced","hypothesis");r$WAIC<-round(m$waic$waic,1);r$model<-"BYM2+RW1+IID+TypeI(iid space-time)";write.csv(r,F9,row.names=FALSE);LOG("S9 done")}}
## ── S10 분모 민감도(E_sens: json offset_sensitivity, 예 METS population) ──
F10<-file.path(RES,"S10_offset_sensitivity_IRR.csv")
if("E_sens"%in%names(DD)&&!file.exists(F10)){m<-timed("S10 offset sensitivity",safefit(MAINV,MAINFAM,ecol="E_sens"));if(!is.null(m)){r<-do.call(rbind,lapply(MAINV,function(v)irr(m,v)));r$label_en<-LAB[r$cov];r$role<-ifelse(r$cov%in%FORCED,"forced","hypothesis");r$WAIC<-round(m$waic$waic,1);r$offset<-INFO$offset_sensitivity;write.csv(r,F10,row.names=FALSE);LOG("S10 done")}}
LOG("ALL DONE %s%s",TAG,SUF)
