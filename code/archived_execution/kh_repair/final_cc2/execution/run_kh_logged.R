# 88k — (2026-10-21 공통판: TAG·FOCAL 환경변수, S20_panel_extra 없으면 건너뜀, 산출 S43_; 쓰는 그룹 nb·struct·kh·combine) 원본 88 — 정신질환(MENTAL, 주분석 판 std) 리뷰어 대비 민감도·보충 분석. 81 의 우도(NB)·강제통제·주 세트(S1_narrative_full)·prior·제약을 그대로 재사용.
#  사용(01_스크립트/ 에서): Rscript 88_mental_sensitivity.R <GROUP>   GROUP = nb | struct | kh | conf | period | maps | combine
#  입력(읽기만): 04_분석/MENTAL/{config,panel,panel_INFO,HM_best}_MENTAL_std.*, results_std/nb_MENTAL.graph, results_std/S9·S10, results_std/S20_panel_extra.csv(88a)
#  산출(새 이름만): results_std/S43_part_<analysis>.csv(분석별) → combine: results_std/S43_sensitivity_IRR.csv(긴 형식) · S43_sensitivity_fit.csv
#                   maps: results_std/S20_obs_fit_resid.csv, figs_std/S20_obs_fit_resid_pooled.png · _windows.png
#  분석 목록(그룹): nb = main_queen(주, 재적합) · knn3·4·5 | struct = ICAR(besag)·iid·no_spatial·PC prior 완화 2종
#    kh = Knorr-Held Type II·III·IV(generic0, Kronecker 구조 + 합 0 제약) | conf = RSR·도로포장률 추가·도시성 5분위·도시성 rw2·관외진료 비율·의원 밀도·둘 다
#    period = 2020-21 제외·2015년 시작·경계변경 4곳 제외·boundary_flag 전체 제외·치매 제외 E(E_std_nondem)
suppressPackageStartupMessages({library(INLA);library(sf);library(spdep);library(jsonlite);library(Matrix)})
args<-commandArgs(trailingOnly=TRUE);GROUP<-args[1];stopifnot(!is.na(GROUP))
# Repaired KH-only drop-in. Fitting requires a deliberate later execution gate.
if (!GROUP %in% c("kh","kh4x")) stop("This reviewed patch covers kh / kh4x only")
if (Sys.getenv("KH_ALLOW_FITS") != "YES_AFTER_REVIEW") stop("No fits authorized by default")
KH_HELPER_DIR <- Sys.getenv("KH_HELPER_DIR")
if (!nzchar(KH_HELPER_DIR)) stop("Set KH_HELPER_DIR to reviewed kh_repair directory")
source(file.path(KH_HELPER_DIR,"kh_constraints.R"))

TAG<-Sys.getenv("TAG",unset="MENTAL");FOCAL<-Sys.getenv("FOCAL",unset="walking_practice_rate");VER<-Sys.getenv("ST_VER",unset="std");SUF<-paste0("_",VER)   # _v 사본(2026-10-19): 판 인자화, 원본 88*.R 은 std 고정 그대로
ROOT<-Sys.getenv("KH_EXTRA_ROOT")
if (!nzchar(ROOT)) stop("Set KH_EXTRA_ROOT explicitly to isolated extra_prepared root")
ROOT<-normalizePath(ROOT,mustWork=TRUE)
if (TAG!="BRCA" || VER!="std") stop("Patch validated for corrected BRCA std only")
SHP<-"${HOME}/Library/CloudStorage/[REDACTED_EMAIL]/내 드라이브/S.K/기타/R_Studio/FBD_DATA_ZIP/final.shp"
OUTD<-file.path(ROOT,"04_분석",TAG);RES<-file.path(OUTD,paste0("results",SUF));FIG<-file.path(OUTD,paste0("figs",SUF))
inla.setOption(num.threads="1:1")   # 1스레드 고정(플레이북 §4)
LOG<-function(...)cat(sprintf("[88 %s %s] ",GROUP,format(Sys.time(),"%H:%M:%S")),sprintf(...),"\n",sep="")
S<-fromJSON(file.path(OUTD,sprintf("config_%s%s.json",TAG,SUF)),simplifyVector=TRUE)
INFO<-fromJSON(file.path(OUTD,sprintf("panel_INFO_%s%s.json",TAG,SUF)))
META<-read.csv(file.path(OUTD,sprintf("panel_META_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);LAB<-setNames(META$label_en,META$variable)
DD<-read.csv(file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);DD$fbd_region<-enc2utf8(DD$fbd_region)
HB<-readLines(file.path(OUTD,sprintf("HM_best_%s%s.txt",TAG,SUF)));kv<-function(k)sub(paste0("^",k,"="),"",grep(paste0("^",k,"="),HB,value=TRUE)[1])
MAINFAM<-kv("main_family");BEST<-kv("best_set");FORCED<-S$forced;SETS<-lapply(S$sets_after_vif,as.character);HYPV<-SETS[[BEST]];MAINV<-unique(c(HYPV,FORCED))
if(file.exists(file.path(RES,"S20_panel_extra.csv"))){EX<-read.csv(file.path(RES,"S20_panel_extra.csv"),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);EX$fbd_region<-enc2utf8(EX$fbd_region)
DD<-merge(DD,EX,by=c("year","fbd_region"),all.x=TRUE,sort=FALSE);stopifnot(nrow(DD)==INFO$n_rows)}
WIN<-INFO$windows

## Frozen, independently reviewed primary graph override (no re-derived queen graph).
G<-st_read(SHP,quiet=TRUE);G$region<-enc2utf8(G$region);n<-nrow(G)
KH_GRAPH_FILE<-Sys.getenv("KH_GRAPH_FILE",unset=file.path(ROOT,"input_graph.graph"))
g<-inla.read.graph(KH_GRAPH_FILE)
kh_assert_dependencies(g,KH_GRAPH_FILE,G$region,
  file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),KH_HELPER_DIR)
DD<-DD[order(DD$year,match(DD$fbd_region,G$region)),];DD$area<-match(DD$fbd_region,G$region);Y0<-min(DD$year);DD$time<-DD$year-Y0+1;DD$yiid<-DD$time;DD$id<-seq_len(nrow(DD))
pcU<-function(u=0.5,a=0.01,pu=0.5,pa=2/3)list(prec=list(prior="pc.prec",param=c(u,a)),phi=list(prior="pc",param=c(pu,pa)))
pc<-pcU();pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
SPf<-function(gr="g",hyper="pc",extra="")sprintf("f(area,model='bym2',graph=%s,scale.model=TRUE,constr=TRUE,hyper=%s%s)",gr,hyper,extra)
TI<-"f(time,model='rw1',scale.model=TRUE,hyper=pcp)";IT<-"f(yiid,model='iid',hyper=pcp)"
bt<-function(v)sprintf("`%s`",v)
fit<-function(vars,data=DD,sp=SPf(),extra=character(0),ecol="E",robust=FALSE,env=parent.frame()){
  rhs<-c("1",bt(vars),extra,sp,TI,IT);rhs<-rhs[nzchar(rhs)];form<-as.formula(paste("cases ~",paste(rhs,collapse="+")),env=env)
  ci<-if(robust)list(control.vb=list(emergency=50),int.strategy="eb",strategy="gaussian") else list(control.vb=list(emergency=50))
  data$EOFF<-data[[ecol]]
  inla(form,family=MAINFAM,data=data,E=EOFF,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE,cpo=FALSE,config=FALSE),control.inla=ci,num.threads="1:1",verbose=FALSE)}
safe<-function(...){
  calls<-list(); strategies<-c("default","gaussian_EB_retry")
  for(k in 1:2){
    warnings_seen<-character(0);err<-NULL;tstart<-as.character(Sys.time())
    LOG("NUMERICAL_ATTEMPT %s start",strategies[k])
    m<-tryCatch(withCallingHandlers(if(k==1)fit(...) else fit(...,robust=TRUE),
      warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));message("INLA_WARNING: ",conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e){err<<-conditionMessage(e);NULL})
    calls[[k]]<-list(strategy=strategies[k],start=tstart,end=as.character(Sys.time()),error=err,warnings=warnings_seen)
    if(!is.null(m)){
      attr(m,"KH_execution_audit")<-list(attempts=calls,selected_strategy=strategies[k],panel_sha256="7d30263a6b6f06633c7a9257a7747711d2da2fbef375a517f1e54a137d722ea4",graph_sha256="8eaecb43a32c724b6a3e8713572fc89e562f4dcc22200552ce5085497dd5edef",years=c(2011,2024),n_obs=3206,threads="1:1",partial_period_fallback=FALSE)
      LOG("NUMERICAL_ATTEMPT %s returned mode.status=%s ok=%s",strategies[k],paste(m$mode$mode.status,collapse=","),paste(m$ok,collapse=","))
      return(m)
    }
    LOG("NUMERICAL_ATTEMPT %s ERROR: %s",strategies[k],err)
  }
  jsonlite::write_json(calls,file.path(Sys.getenv("KH_RUN_LOGDIR"),paste0("failed_attempts_",format(Sys.time(),"%Y%m%d_%H%M%S"),".json")),auto_unbox=TRUE,pretty=TRUE)
  NULL
}
irr<-function(m,v){r<-m$summary.fixed[v,,drop=FALSE];data.frame(cov=v,IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),cred=(r$`0.025quant`>0|r$`0.975quant`<0),stringsAsFactors=FALSE)}
hp<-function(m,pat){h<-m$summary.hyperpar;k<-grep(pat,rownames(h));if(length(k))round(h[k[1],"mean"],3) else NA}
save_part<-function(an,grp,m,vars,data,t0,note="",cmp="",added=character(0)){
  if(is.null(m)){write.csv(data.frame(analysis=an,group=grp,cov=NA,IRR=NA,lo=NA,hi=NA,cred=NA,WAIC=NA,note=paste("FIT FAILED",note)),file.path(RES,sprintf("S43_part_%s.csv",an)),row.names=FALSE);LOG("%s FAILED",an);return(invisible())}
  saveRDS(m,file.path(RES,sprintf("S43_fit_%s.rds",an)))
  note<-paste(note,"[KH repair v1; contact593 cc2 graph; component constraints; no period fallback]")
  r<-do.call(rbind,lapply(vars,function(v)irr(m,v)));r$role<-ifelse(r$cov%in%HYPV,"hypothesis",ifelse(r$cov%in%FORCED,"forced","added"))
  out<-data.frame(analysis=an,group=grp,r,WAIC=round(m$waic$waic,1),DIC=round(m$dic$dic,1),n_obs=sum(!is.na(data$cases)),n_regions=length(unique(data$area[!is.na(data$cases)])),
    years=sprintf("%d-%d",min(data$year[!is.na(data$cases)]),max(data$year[!is.na(data$cases)])),phi=hp(m,"^Phi for area"),comparable_group=cmp,
    minutes=round(as.numeric(difftime(Sys.time(),t0,units="mins")),1),note=note,stringsAsFactors=FALSE)
  write.csv(out,file.path(RES,sprintf("S43_part_%s.csv",an)),row.names=FALSE,fileEncoding="UTF-8")
  LOG("%-22s WAIC %9.1f  cred hyp %2d/%d  walking %.3f (%.3f-%.3f)  · %.1f min",an,out$WAIC[1],sum(r$cred[r$role=="hypothesis"]),length(HYPV),r$IRR[r$cov==FOCAL],r$lo[r$cov==FOCAL],r$hi[r$cov==FOCAL],out$minutes[1])}
done<-function(an){
  if(file.exists(file.path(RES,sprintf("S43_part_%s.csv",an))))
    stop("Existing output needs provenance review; refusing silent skip: ",an)
  FALSE
}
run<-function(an,grp,expr_fun,vars=MAINV,data=DD,note="",cmp=""){if(done(an)){LOG("%s exists, skip",an);return(invisible())};t0<-Sys.time();m<-expr_fun();save_part(an,grp,m,vars,data,t0,note,cmp);invisible(m)}

## ───────── 1. 이웃 정의 ─────────
if(GROUP=="nb"){
  m<-run("main_queen","1_neighbour",function()safe(MAINV),note="주모형 재적합(S3 와 같은 사양)",cmp="neighbour")
  if(!is.null(m)){mu<-m$summary.fitted.values$mean   # 지도(9)용 적합 RR 저장
    fx<-m$summary.fixed$mean;names(fx)<-gsub("`","",rownames(m$summary.fixed));X<-as.matrix(DD[,MAINV]);lpf<-fx["(Intercept)"]+X%*%fx[MAINV]
    write.csv(data.frame(year=DD$year,fbd_region=DD$fbd_region,window=DD$window,cases=DD$cases,E=DD$E,fitted_RR=mu,fixed_lp=as.numeric(lpf)),file.path(RES,"S20_fitted_main.csv"),row.names=FALSE,fileEncoding="UTF-8")}
  for(k in 3:5){an<-sprintf("knn%d",k);if(done(an))next
    nbk<-make.sym.nb(knn2nb(knearneigh(cen,k=k)));gf<-file.path(RES,sprintf("S20_nb_knn%d.graph",k));nb2INLA(gf,nbk);gk<-inla.read.graph(gf)
    run(an,"1_neighbour",function()safe(MAINV,sp=SPf("gk")),note=sprintf("KNN k=%d(중심점, 대칭화), 평균 이웃 %.2f, 부분그래프 %d",k,mean(card(nbk)),n.comp.nb(nbk)$nc),cmp="neighbour")}
}

## ───────── 2. 공간 구조·prior ─────────
if(GROUP=="struct"){
  run("icar_besag","2_structure",function()safe(MAINV,sp="f(area,model='besag',graph=g,scale.model=TRUE,constr=TRUE,hyper=pcp)"),note="ICAR(besag)만",cmp="structure")
  run("iid_only","2_structure",function()safe(MAINV,sp="f(area,model='iid',hyper=pcp)"),note="공간 iid 만",cmp="structure")
  run("no_spatial","2_structure",function()safe(MAINV,sp=""),note="공간효과 없음(RW1+연도 iid 만)",cmp="structure")
  pcA<-pcU(1,0.01,0.5,0.5);pcB<-pcU(2,0.05,0.5,0.5);pcC<-pcU(0.5,0.01,0.5,0.5)
  run("pc_phi_flat","2_structure",function()safe(MAINV,sp=SPf(hyper="pcC")),note="BYM2 PC: prec(0.5,0.01) 유지, phi P(phi<0.5)=0.5(중립)",cmp="structure")
  run("pc_relaxA","2_structure",function()safe(MAINV,sp=SPf(hyper="pcA")),note="BYM2 PC 완화 A: prec P(sd>1)=0.01, phi P(phi<0.5)=0.5",cmp="structure")
  run("pc_relaxB","2_structure",function()safe(MAINV,sp=SPf(hyper="pcB")),note="BYM2 PC 완화 B: prec P(sd>2)=0.05, phi P(phi<0.5)=0.5",cmp="structure")
}

## ───────── 3. Knorr-Held II·III·IV ─────────
if(GROUP=="kh"){
  khfit<-function(type,d) kh_fit_repaired(type,d,MAINV)
  for(ty in c("II","III","IV")){an<-sprintf("KH_type%s",ty);if(done(an))next
    inla.setOption(inla.timeout=if(ty=="IV")5400 else 0)
    t0<-Sys.time();m<-khfit(ty,DD);note<-sprintf("Knorr-Held Type %s(BYM2+RW1+iid 위에 공간×시간 상호작용, 합 0 제약)",ty);d<-DD
    save_part(an,"3_spacetime",m,MAINV,d,t0,note,"spacetime")}
  inla.setOption(inla.timeout=0)
}


## ───────── kh4x: Knorr-Held IV 주모형 후보 추가 점검(2026-10-21 /333 반영) ─────────
if(GROUP=="kh4x"){
  khfit<-function(type,d) kh_fit_repaired(type,d,MAINV)
  ROBV<-{x<-read.csv(file.path(RES,"S7_multiverse_summary.csv"),stringsAsFactors=FALSE);x$cov[x$robust]}
  khv<-function(vars,d) kh_fit_repaired("IV",d,vars)
  inla.setOption(inla.timeout=5400)
  SUB<-Sys.getenv("SUB","all")
  if(SUB%in%c("f","g")) stop("pop80 source/sex definition requires separate reviewed BRCA dependency; legacy MENTAL source blocked")
  if(SUB=="n" && file.exists(file.path(RES,"S46_KH4_spatialRR.csv")))
    stop("Existing S46 needs provenance review")
  if(SUB%in%c("all","a")&&!done("KH4_bw_twoway")){d<-DD;BW<-character(0)
    for(v in ROBV){mi<-ave(d[[v]],d$area,FUN=mean);mt<-ave(d[[v]],d$year,FUN=mean);d[[paste0(v,"__between")]]<-mi;d[[paste0(v,"__within")]]<-d[[v]]-mi-mt+mean(d[[v]]);BW<-c(BW,paste0(v,c("__between","__within")))}
    vv<-c(setdiff(MAINV,ROBV),BW);t0<-Sys.time();save_part("KH4_bw_twoway","3_spacetime",khv(vv,d),vv,d,t0,"KH IV + 강건 집합 이원 중심화 분해")}
  dsl<-grep("dyslip",MAINV,value=TRUE);scr<-grep("screening_uptake",MAINV,value=TRUE)
  if(SUB%in%c("all","b")&&length(dsl)&&!done("KH4_drop_dyslipidemia")){vv<-setdiff(MAINV,dsl);t0<-Sys.time();save_part("KH4_drop_dyslipidemia","3_spacetime",khv(vv,DD),vv,DD,t0,"KH IV, 이상지질혈증(2018 정의 변경) 제외")}
  if(SUB%in%c("all","c")&&length(scr)&&!done("KH4_drop_screening")){vv<-setdiff(MAINV,scr);t0<-Sys.time();save_part("KH4_drop_screening","3_spacetime",khv(vv,DD),vv,DD,t0,"KH IV, 검진 수검률 강제통제 제외")}
  if(SUB%in%c("all","d")&&!done("KH4_S5_reproductive")&&!is.null(SETS[["S5_reproductive"]])){vv<-unique(c(SETS[["S5_reproductive"]],FORCED));t0<-Sys.time();save_part("KH4_S5_reproductive","3_spacetime",khv(vv,DD),vv,DD,t0,"KH IV, 간명 세트 S5")}
  if(SUB%in%c("all","e")&&length(scr)&&!done("T0_drop_screening")){vv<-setdiff(MAINV,scr);t0<-Sys.time();save_part("T0_drop_screening","1_base",safe(vv,data=DD),vv,DD,t0,"기본형, 검진 수검률 강제통제 제외")}
  P80<-file.path(ROOT,"04_분석","MENTAL","results_std3","S23_panel_extra.csv")
  if(SUB%in%c("f","g")){X<-read.csv(P80,fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);X$fbd_region<-enc2utf8(X$fbd_region);d80<-merge(DD,X[,c("year","fbd_region","pop80_share")],by=c("year","fbd_region"),all.x=TRUE,sort=FALSE);d80<-d80[order(d80$id),];stopifnot(!anyNA(d80$pop80_share));vv<-c(MAINV,"pop80_share")
    if(SUB=="f"&&!done("T0_add_pop80")){t0<-Sys.time();save_part("T0_add_pop80","1_base",safe(vv,data=d80),vv,d80,t0,"기본형 + 80세 이상 인구 비율 강제통제(std3형, 전체 인구 기준)")}
    if(SUB=="g"&&!done("KH4_add_pop80")){t0<-Sys.time();save_part("KH4_add_pop80","3_spacetime",khv(vv,d80),vv,d80,t0,"KH IV + 80세 이상 인구 비율 강제통제")}}
  if(SUB%in%c("j","k")){X<-read.csv(file.path(RES,"S45_panel_longlag.csv"),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);X$fbd_region<-enc2utf8(X$fbd_region);dl<-merge(DD,X,by=c("year","fbd_region"),all.x=TRUE,sort=FALSE);dl<-dl[order(dl$id),]
    FV2<-c("total_fertility_rate","mean_age_first_marriage_wife")
    if(SUB=="j"&&!done("T0_lag10_fertility_marriage")){vv<-c(setdiff(MAINV,FV2),paste0(FV2,"__lag10"));t0<-Sys.time();save_part("T0_lag10_fertility_marriage","1_base",safe(vv,data=dl),vv,dl,t0,"합계출산율·초혼연령을 10년 전 값으로(2001-2014)")}
    if(SUB=="k"&&!done("T0_base0004_fertility_marriage")){vv<-c(setdiff(MAINV,FV2),paste0(FV2,"__base0004"));t0<-Sys.time();save_part("T0_base0004_fertility_marriage","1_base",safe(vv,data=dl),vv,dl,t0,"합계출산율·초혼연령을 2000-2004 지역 평균(고정 노출)으로")}}
  if(SUB%in%c("l","m")){X<-read.csv(file.path(RES,"S45_panel_longlag.csv"),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);X$fbd_region<-enc2utf8(X$fbd_region);dl<-merge(DD,X,by=c("year","fbd_region"),all.x=TRUE,sort=FALSE);dl<-dl[order(dl$id),]
    FV2<-c("total_fertility_rate","mean_age_first_marriage_wife");suf<-if(SUB=="l")"__lag10" else "__base0004";an<-paste0("KH4",suf)
    if(!done(an)){vv<-c(setdiff(MAINV,FV2),paste0(FV2,suf));t0<-Sys.time();save_part(an,"3_spacetime",khv(vv,dl),vv,dl,t0,paste("KH IV +",suf))}}
  if(SUB=="n"&&!file.exists(file.path(RES,"S46_KH4_spatialRR.csv"))){m<-khv(MAINV,DD);if(is.null(m))stop("KH4 spatial fit failed; no partial-period fallback");saveRDS(m,file.path(RES,"S46_KH4_full_fit.rds"));sr<-m$summary.random$area
    p1<-sapply(1:n,function(i)1-inla.pmarginal(0,m$marginals.random$area[[i]]))
    out<-data.frame(fbd_region=G$region,RR=round(exp(sr$mean[1:n]),4),P_gt1=round(p1,4),class=ifelse(p1>0.95,"hotspot",ifelse(p1<0.05,"coldspot","neutral")))
    write.csv(out,file.path(RES,"S46_KH4_spatialRR.csv"),row.names=FALSE,fileEncoding="UTF-8");LOG("KH4 spatial: hot %d cold %d",sum(out$class=="hotspot"),sum(out$class=="coldspot"))}
  if(SUB=="h"&&!done("T0_obs_2015_noemp")){d<-DD[DD$year>=2015,];d$time<-d$year-2014;d$yiid<-d$time;vv<-setdiff(MAINV,c("regular_employee_share","female_employment_rate"))
    t0<-Sys.time();save_part("T0_obs_2015_noemp","1_base",safe(vv,data=d),vv,d,t0,"2015-2024(1인가구 전 연도 관측), 대체 비중 큰 고용 2변수 제외")}
  if(SUB=="i"){FV<-intersect(c("total_fertility_rate","mean_age_first_marriage_wife","fertility_rate_35_39","mean_maternal_age"),MAINV)
    for(fv in FV){an<-paste0("T0_single_",fv);if(done(an))next;vv<-c(setdiff(MAINV,setdiff(FV,fv)));t0<-Sys.time();save_part(an,"1_base",safe(vv,data=DD),vv,DD,t0,paste("출산 변수 하나만:",fv))}}
  inla.setOption(inla.timeout=0)
}

## ───────── 4~6. 공간 교란·음성대조·도시성·의료 이용 ─────────
if(GROUP=="conf"){
  # RSR: BYM2 총 공간효과(잠재벡터 앞 n 개)를 시군구 평균 공변량(+절편은 constr)에 직교하도록 추가 제약
  if(!done("RSR")){Xb<-apply(as.matrix(DD[,MAINV]),2,function(x)tapply(x,DD$area,mean))[as.character(1:n),];Xb<-scale(Xb,scale=FALSE)
    RSRC<<-list(A=cbind(t(Xb),matrix(0,ncol(Xb),n)),e=rep(0,ncol(Xb)))
    run("RSR","4_spatial_confounding",function()safe(MAINV,sp=SPf(extra=",extraconstr=RSRC"),env=.GlobalEnv),note=sprintf("RSR: BYM2 총효과 ⟂ 시군구 평균 공변량 %d개(제약)",ncol(Xb)),cmp="main_data")}
  run("negctl_road","5_negcontrol",function()safe(c(MAINV,"road_pavement_rate")),vars=c(MAINV,"road_pavement_rate"),note="도로포장률(음성대조)을 통제로 추가",cmp="main_data")
  DD$urban_q5f<-factor(DD$urban_q5);for(q in 2:5)DD[[sprintf("urbQ%d",q)]]<-as.integer(DD$urban_q5==q)
  V0<-setdiff(MAINV,"urban_population_share")
  run("urban_quintile","5_negcontrol",function()safe(c(V0,sprintf("urbQ%d",2:5))),vars=c(V0,sprintf("urbQ%d",2:5)),note="도시지역 인구 비율 선형항 → 5분위 범주(기준 Q1)",cmp="main_data",data=DD)
  run("urban_rw2","5_negcontrol",function()safe(V0,extra="f(urban_rw2_bin,model='rw2',scale.model=TRUE,constr=TRUE,hyper=pcp)",data=DD),vars=V0,note="도시지역 인구 비율 선형항 → 20분위 rw2",cmp="main_data",data=DD)
  run("urban_q5_plus_road","5_negcontrol",function()safe(c(V0,sprintf("urbQ%d",2:5),"road_pavement_rate")),vars=c(V0,sprintf("urbQ%d",2:5),"road_pavement_rate"),note="도시성 5분위 + 도로포장률",cmp="main_data",data=DD)
  run("access_outside","6_care_use",function()safe(c(MAINV,"outside_care_share")),vars=c(MAINV,"outside_care_share"),note="관외진료 비율(거주지 기준 전체 질환 관외/(관내+관외) 진료실인원, 1 SD) 추가",cmp="main_data")
  run("access_clinics","6_care_use",function()safe(c(MAINV,"clinics_per_100k")),vars=c(MAINV,"clinics_per_100k"),note="의원 수/10만(1 SD) 추가",cmp="main_data")
  run("access_both","6_care_use",function()safe(c(MAINV,"outside_care_share","clinics_per_100k")),vars=c(MAINV,"outside_care_share","clinics_per_100k"),note="관외진료 비율 + 의원 밀도 추가",cmp="main_data")
}

## ───────── 7~8. 기간·지역·치매 offset ─────────
if(GROUP=="period"){
  d<-DD;d$cases[d$year%in%2020:2021]<-NA
  run("excl_2020_2021","7_period_region",function()safe(MAINV,data=d),data=d,note="2020-2021 결측 처리(시간 색인 유지)",cmp="excl_2020_2021")
  d<-DD[DD$year>=2015,];d$time<-d$year-2014;d$yiid<-d$time
  run("start_2015","7_period_region",function()safe(MAINV,data=d),data=d,note="2015-2024만(공변량은 전체 기간 z 그대로)",cmp="start_2015")
  d<-DD;d$cases[d$boundary_change==1]<-NA
  run("excl_boundary_change","7_period_region",function()safe(MAINV,data=d),data=d,note="경계변경 4곳(세종·군위·청주·창원) 결측 처리(그래프 유지)",cmp="excl_bch")
  d<-DD;d$cases[d$boundary_any==1]<-NA
  run("excl_boundary_any","7_period_region",function()safe(MAINV,data=d),data=d,note=sprintf("boundary_flag 있는 시군구 전체 + 세종·군위 %d곳 결측 처리(개명·표기정정·하위구 합산 포함)",length(unique(DD$area[DD$boundary_any==1]))),cmp="excl_bany")
  run("dementia_E_nondem","8_dementia",function()safe(MAINV,ecol="E_sens"),note=sprintf("offset = %s(치매 112 제외 연령 프로파일, S10 재현)",INFO$offset_sensitivity),cmp="E_nondem")
}

## ───────── 9. 관측·적합·잔차 지도 ─────────
if(GROUP=="maps"){suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(patchwork)})
  F<-read.csv(file.path(RES,"S20_fitted_main.csv"),fileEncoding="UTF-8",stringsAsFactors=FALSE);F$fbd_region<-enc2utf8(F$fbd_region);F$mu<-F$fitted_RR*F$E
  per<-c(pooled="2011-2024",setNames(sapply(names(WIN),function(k)sprintf("%d-%d",WIN[[k]][1],WIN[[k]][2])),names(WIN)))
  agg<-do.call(rbind,lapply(names(per),function(p){x<-if(p=="pooled")F else F[F$window==p,]
    # 공변량 예측 RR: exp(고정효과 선형예측)을 기간 내 Σ E·pred = Σ cases 로 재척도(무작위효과 제외)
    x$pf<-exp(x$fixed_lp);x$pf<-x$pf*sum(x$cases)/sum(x$E*x$pf);x$Epf<-x$E*x$pf   # (dplyr 함정: summarise 안에서 E 를 먼저 합치면 뒤의 E*pf 가 합계 E 를 씀 → 행 단위로 미리 계산)
    a<-x%>%group_by(fbd_region)%>%summarise(cases=sum(cases),Esum=sum(E),mu=sum(mu),Epf=sum(Epf),.groups="drop")%>%rename(E=Esum)
    data.frame(period=p,years=per[[p]],fbd_region=a$fbd_region,cases=a$cases,E=round(a$E,1),SIR_obs=round(a$cases/a$E,4),RR_fit=round(a$mu/a$E,4),RR_cov=round(a$Epf/a$E,4),
      resid_fit=round(a$cases/a$mu,4),resid_cov=round(a$cases/a$Epf,4),stringsAsFactors=FALSE)}))
  write.csv(agg,file.path(RES,"S20_obs_fit_resid.csv"),row.names=FALSE,fileEncoding="UTF-8")
  sm<-agg%>%group_by(period)%>%summarise(cor_obs_fit=round(cor(log(SIR_obs),log(RR_fit)),3),cor_obs_cov=round(cor(log(SIR_obs),log(RR_cov)),3),resid_fit_min=min(resid_fit),resid_fit_max=max(resid_fit),resid_cov_p10=round(quantile(resid_cov,.1),3),resid_cov_p90=round(quantile(resid_cov,.9),3),.groups="drop")
  write.csv(sm,file.path(RES,"S20_obs_fit_resid_summary.csv"),row.names=FALSE,fileEncoding="UTF-8");print(as.data.frame(sm))
  RRB<-c(-Inf,0.67,0.8,0.9,1,1.1,1.25,1.5,Inf);RRLAB<-c("<0.67","0.67-0.8","0.8-0.9","0.9-1","1-1.1","1.1-1.25","1.25-1.5",">=1.5")
  PAL8<-setNames(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"),RRLAB)
  RB<-c(-Inf,0.95,0.97,0.99,1,1.01,1.03,1.05,Inf);RLAB<-c("<0.95","0.95-0.97","0.97-0.99","0.99-1","1-1.01","1.01-1.03","1.03-1.05",">=1.05");PALR<-setNames(unname(PAL8),RLAB)
  G$sido<-substr(as.character(G$sigungu_cd),1,2);SIDO<-suppressWarnings(G%>%group_by(sido)%>%summarise(.groups="drop"))
  thm<-theme_void(base_size=10)+theme(strip.text=element_text(face="bold",size=9,hjust=0),legend.title=element_text(face="bold",size=8),legend.text=element_text(size=7),plot.margin=margin(4,4,4,4))
  lv<-function(SP,labs,col){d<-SP[seq_along(labs),];d[[col]]<-factor(labs,levels=labs);geom_sf(data=d,aes(fill=.data[[col]]),alpha=0,color=NA,show.legend=TRUE,inherit.aes=FALSE)}
  mk<-function(pp){SP<-merge(G[,c("region","geometry")],agg[agg$period%in%pp,],by.x="region",by.y="fbd_region")
    L<-data.frame(var=c("SIR_obs","RR_fit","RR_cov","resid_fit"),lab=c("Observed SIR","Fitted RR (full model)","Covariate-predicted RR (fixed effects)","Residual ratio, observed / fitted"),stringsAsFactors=FALSE)
    plots<-list()
    for(p in pp)for(i in 1:4){v<-L$var[i];S2<-SP[SP$period==p,];res<-v=="resid_fit"
      S2$cl<-factor(cut(S2[[v]],if(res)RB else RRB,labels=if(res)RLAB else RRLAB),levels=if(res)RLAB else RRLAB);S2$panel<-sprintf("%s, %s",L$lab[i],per[[p]])
      plots[[length(plots)+1]]<-ggplot(S2)+geom_sf(aes(fill=cl),color="grey70",linewidth=0.05)+lv(S2,if(res)RLAB else RRLAB,"cl")+geom_sf(data=SIDO,fill=NA,color="grey35",linewidth=0.15,inherit.aes=FALSE)+
        facet_wrap(~panel)+scale_fill_manual(values=if(res)PALR else PAL8,breaks=rev(if(res)RLAB else RRLAB),limits=if(res)RLAB else RRLAB,drop=FALSE,name=if(res)"Observed /\nfitted" else "SIR or RR")+
        guides(fill=guide_legend(nrow=2,override.aes=list(alpha=1,color="grey40")))+thm+theme(legend.position="bottom")}
    wrap_plots(plots,ncol=4)}
  ggsave(file.path(FIG,"S20_obs_fit_resid_pooled.png"),mk("pooled"),width=15,height=5.4,dpi=300,bg="white")
  ggsave(file.path(FIG,"S20_obs_fit_resid_windows.png"),mk(names(WIN)),width=15,height=15,dpi=300,bg="white")
  LOG("maps saved")
}

## ───────── combine ─────────
if(GROUP=="combine"){
  fs<-list.files(RES,pattern="^S43_part_.*\\.csv$",full.names=TRUE);P<-do.call(rbind,lapply(fs,function(f)read.csv(f,fileEncoding="UTF-8",stringsAsFactors=FALSE)))
  ad<-function(f,an,grp,note,cmp){x<-read.csv(file.path(RES,f),stringsAsFactors=FALSE);data.frame(analysis=an,group=grp,cov=x$cov,IRR=x$IRR,lo=x$lo,hi=x$hi,cred=x$cred,
    role=ifelse(x$cov%in%HYPV,"hypothesis",ifelse(x$cov%in%FORCED,"forced","added")),WAIC=x$WAIC,DIC=NA,n_obs=INFO$n_rows,n_regions=INFO$n_sgg,years="2011-2024",phi=NA,comparable_group=cmp,minutes=NA,note=note,stringsAsFactors=FALSE)}
  P<-rbind(P,ad("S9_typeI_IRR.csv","KH_typeI_S9","3_spacetime","Knorr-Held Type I(기존 S9, 81)","spacetime"),ad("S10_offset_sensitivity_IRR.csv","dementia_E_nondem_S10","8_dementia",sprintf("기존 S10(81) offset %s",INFO$offset_sensitivity),"E_nondem"))
  # 비교 가능 집단(같은 자료·같은 offset): 주모형 자료 = main_data + neighbour + structure + spacetime
  P$cmp2<-ifelse(P$comparable_group%in%c("neighbour","structure","spacetime","main_data"),"main_data",P$comparable_group)
  w<-unique(P[,c("analysis","comparable_group","cmp2","WAIC")]);mainW<-w$WAIC[w$analysis=="main_queen"]
  P$dWAIC_vs_main<-ifelse(P$cmp2=="main_data",round(P$WAIC-mainW,1),NA)
  P$min_WAIC_neighbour<-P$comparable_group=="neighbour"&P$WAIC==min(w$WAIC[w$comparable_group=="neighbour"],na.rm=TRUE)
  sw<-w[w$comparable_group=="structure"|w$analysis=="main_queen",];P$min_WAIC_structure<-P$analysis==sw$analysis[which.min(sw$WAIC)]
  M<-read.csv(file.path(RES,"S3_main_IRR.csv"),stringsAsFactors=FALSE);mm<-setNames(log(M$IRR),M$cov);mc<-setNames(M$cred,M$cov)
  P$same_sign_as_main<-ifelse(P$cov%in%M$cov,sign(log(P$IRR))==sign(mm[P$cov])|P$IRR==1,NA);P$main_cred<-mc[P$cov]
  ROB<-read.csv(file.path(RES,"S7_multiverse_summary.csv"),stringsAsFactors=FALSE);ROBV<-ROB$cov[ROB$robust];P$robust8<-P$cov%in%ROBV
  P<-P[order(P$group,P$analysis,match(P$role,c("hypothesis","forced","added")),match(P$cov,c(HYPV,FORCED))),]
  P$cmp2<-NULL;write.csv(P,file.path(RES,"S43_sensitivity_IRR.csv"),row.names=FALSE,fileEncoding="UTF-8")
  # 분석별 요약
  sm<-do.call(rbind,lapply(split(P,P$analysis),function(x){h<-x[x$role=="hypothesis",];r<-h[h$robust8,];wk<-h[h$cov==FOCAL,]
    data.frame(analysis=x$analysis[1],group=x$group[1],WAIC=x$WAIC[1],dWAIC_vs_main=x$dWAIC_vs_main[1],comparable_group=x$comparable_group[1],n_obs=x$n_obs[1],years=x$years[1],phi=x$phi[1],
      n_cred_hyp=sum(h$cred,na.rm=TRUE),robust8_cred_same_sign=sum(r$cred&r$same_sign_as_main,na.rm=TRUE),robust8_lost=paste(r$cov[!(r$cred&r$same_sign_as_main)],collapse=";"),
      walking=if(nrow(wk))sprintf("%.3f (%.3f-%.3f)%s",wk$IRR,wk$lo,wk$hi,if(isTRUE(wk$cred))"*" else "") else NA,
      sign_flips_vs_main=paste(h$cov[!is.na(h$same_sign_as_main)&!h$same_sign_as_main&h$cred],collapse=";"),
      min_WAIC_neighbour=any(x$min_WAIC_neighbour),min_WAIC_structure=any(x$min_WAIC_structure),note=x$note[1],stringsAsFactors=FALSE)}))
  sm<-sm[order(sm$group,sm$analysis),];write.csv(sm,file.path(RES,"S43_sensitivity_fit.csv"),row.names=FALSE,fileEncoding="UTF-8");print(sm[,c("analysis","WAIC","dWAIC_vs_main","n_cred_hyp","robust8_cred_same_sign","walking","robust8_lost")])
}
LOG("done")
