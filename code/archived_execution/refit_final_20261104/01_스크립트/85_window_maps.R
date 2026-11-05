# 85 — 창별 공간 RR·초과확률(핫/콜드) 지도 + 창별 연도별 Moran (JKMS 01 §3 T17·T28 창별판, 저자 지시 2026-10-19)
#  사용: (01_스크립트/ 에서) Rscript 85_window_maps.R <TAG> <VER>     VER = v1(=results/, figs/) | std | v2 ...(=results_<VER>/, figs_<VER>/)
#  전제: 81 이 끝나 HM_best_<TAG><SUF>.txt·S1_likelihood.csv·S4_window_IRR.csv 가 있음(주 세트·우도를 그대로 씀). 81·80 은 읽기만 함.
#  모형: 81 S4 와 동일 — cases ~ offset(log E) + forced + 주 세트 + BYM2(PC, queen+섬 최근접) + RW1(연도) + iid(연도), 창마다 재적합, INLA 1스레드.
#  산출: results<SUF>/S4b_window_spatialRR.csv(창×시군구 RR·CrI·P(RR>1)·class·창 관측 SIR) · S4b_window_fit.csv(창별 핫/콜드 수·phi·WAIC)
#        · S4b_window_moran.csv(창별 연도별 + 창 합산 log SIR Moran I, 999 순열) · figs<SUF>/<TAG>_F4b_window_spatialRR.png · <TAG>_F4c_window_exceedance.png
suppressPackageStartupMessages({library(INLA);library(sf);library(spdep);library(jsonlite);library(ggplot2);library(dplyr)})
args<-commandArgs(trailingOnly=TRUE);TAG<-args[1];VER<-if(length(args)>=2)args[2] else "v1";stopifnot(!is.na(TAG))
SUF<-if(VER %in% c("","v1"))"" else paste0("_",VER)
ROOT<-normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE))),".."))
SHP<-"${HOME}/Library/CloudStorage/[REDACTED_EMAIL]/내 드라이브/S.K/기타/R_Studio/FBD_DATA_ZIP/final.shp"
OUTD<-file.path(ROOT,"04_분석",TAG);RES<-file.path(OUTD,paste0("results",SUF));FIG<-file.path(OUTD,paste0("figs",SUF));dir.create(FIG,showWarnings=FALSE)
inla.setOption(num.threads="1:1")   # 1스레드 고정(체육 플레이북 §4)
LOG<-function(...)cat(sprintf("[85 %s %s%s] ",format(Sys.time(),"%H:%M:%S"),TAG,SUF),sprintf(...),"\n",sep="")
S<-fromJSON(file.path(OUTD,sprintf("config_%s%s.json",TAG,SUF)),simplifyVector=TRUE)
INFO<-fromJSON(file.path(OUTD,sprintf("panel_INFO_%s%s.json",TAG,SUF)))
DD<-read.csv(file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);DD$fbd_region<-enc2utf8(DD$fbd_region)
HB<-readLines(file.path(OUTD,sprintf("HM_best_%s%s.txt",TAG,SUF)));kv<-function(k)sub(paste0("^",k,"="),"",grep(paste0("^",k,"="),HB,value=TRUE)[1])
MAINFAM<-kv("main_family");BEST<-kv("best_set");FORCED<-S$forced;SETS<-lapply(S$sets_after_vif,as.character)
MAINV<-unique(c(SETS[[BEST]],FORCED));stopifnot(all(MAINV%in%names(DD)))
WIN<-INFO$windows;LOG("family=%s set=%s nvar=%d windows=%s",MAINFAM,BEST,length(MAINV),paste(names(WIN),collapse=","))
## 그래프: 81 과 같은 규칙(queen + 고립 폴리곤 최근접 1개 대칭 연결), 81 이 저장한 nb_<TAG>.graph 와 대조
G<-st_read(SHP,quiet=TRUE);G$region<-enc2utf8(G$region);n<-nrow(G)
nb<-poly2nb(G,queen=TRUE);isl<-which(card(nb)==0)
if(length(isl)){cen<-st_coordinates(st_centroid(st_geometry(G)));kn<-knearneigh(cen,k=1)$nn
  for(i in isl){j<-kn[i,1];nb[[i]]<-sort(unique(c(nb[[i]][nb[[i]]>0],j)));nb[[j]]<-sort(unique(c(nb[[j]][nb[[j]]>0],i)))}}
GF<-file.path(RES,sprintf("nb_%s.graph",TAG));g<-inla.read.graph(GF);stopifnot(g$n==n)
for(i in seq_len(n))nb[[i]]<-as.integer(sort(g$nbs[[i]]))
attr(nb,"sym")<-TRUE
stopifnot(all(g$nnbs==card(nb)))
stopifnot(all(vapply(seq_len(n),function(i) identical(sort(as.integer(g$nbs[[i]])),sort(as.integer(nb[[i]]))),logical(1))))
OBSR<-G$region%in%unique(DD$fbd_region)
moran_obs<-function(x,regions){k<-!is.na(x)&is.finite(x);if(sum(k)<10)return(list(statistic=NA,p.value=NA));kk<-G$region%in%regions[k];nbk<-subset(nb,kk);lwk<-nb2listw(nbk,style="W",zero.policy=TRUE);moran.mc(x[k][match(G$region[kk],regions[k])],lwk,999,zero.policy=TRUE)}
DD$area<-match(DD$fbd_region,G$region)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
form<-as.formula(paste("cases ~ 1 +",paste(sprintf("`%s`",MAINV),collapse="+"),"+ f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc) + f(time,model='rw1',scale.model=TRUE,hyper=pcp) + f(yiid,model='iid',hyper=pcp)"))
fitw<-function(d,robust=FALSE){ci<-if(robust)list(control.vb=list(emergency=50),int.strategy="eb",strategy="gaussian") else list(control.vb=list(emergency=50))
  inla(form,family=MAINFAM,data=d,E=d$E,control.fixed=CF,control.compute=list(waic=TRUE,dic=FALSE,cpo=FALSE,config=FALSE),control.inla=ci,verbose=FALSE)}
F1<-file.path(RES,"S4b_window_spatialRR.csv");F2<-file.path(RES,"S4b_window_fit.csv");F3<-file.path(RES,"S4b_window_moran.csv")
RRL<-NULL;FIT<-NULL;MOR<-NULL;set.seed(20261019)
for(wn in names(WIN)){yy<-WIN[[wn]];d<-DD[DD$year>=yy[1]&DD$year<=yy[2],];d$time<-d$year-min(d$year)+1;d$yiid<-d$time
  t0<-Sys.time();m<-tryCatch(fitw(d),error=function(e){LOG("retry robust: %s",conditionMessage(e));fitw(d,TRUE)})
  ur<-m$summary.random$area[1:n,];pep<-sapply(m$marginals.random$area[1:n],function(mg)1-inla.pmarginal(0,mg))
  ag<-aggregate(cbind(cases,E)~fbd_region,d,sum)
  rr<-data.frame(window=wn,years=sprintf("%d-%d",yy[1],yy[2]),area=1:n,fbd_region=G$region,sigungu_cd=G$sigungu_cd,RR=round(exp(ur$mean),4),lo=round(exp(ur$`0.025quant`),4),hi=round(exp(ur$`0.975quant`),4),P_gt1=round(pep,4))
  rr$class<-ifelse(rr$P_gt1>0.95,"hotspot",ifelse(rr$P_gt1<0.05,"coldspot","neutral"));rr$class[!OBSR]<-"unobserved"
  rr$SIR_obs<-round((ag$cases/ag$E)[match(rr$fbd_region,ag$fbd_region)],4);RRL<-rbind(RRL,rr)
  hp<-m$summary.hyperpar;phi<-if("Phi for area"%in%rownames(hp))hp["Phi for area",] else NULL
  FIT<-rbind(FIT,data.frame(window=wn,years=sprintf("%d-%d",yy[1],yy[2]),family=MAINFAM,set=BEST,n_obs=nrow(d),cases=sum(d$cases),WAIC=round(m$waic$waic,1),
    phi=round(if(is.null(phi))NA else phi$mean,3),phi_lo=round(if(is.null(phi))NA else phi$`0.025quant`,3),phi_hi=round(if(is.null(phi))NA else phi$`0.975quant`,3),
    n_hotspot=sum(rr$class=="hotspot"),n_coldspot=sum(rr$class=="coldspot"),RR_min=min(rr$RR[OBSR]),RR_max=max(rr$RR[OBSR]),RR_p10=round(quantile(rr$RR[OBSR],0.1),3),RR_p90=round(quantile(rr$RR[OBSR],0.9),3),
    top5_RR=paste(head(rr$fbd_region[order(-rr$RR)],5),collapse=";"),bottom5_RR=paste(head(rr$fbd_region[order(rr$RR)],5),collapse=";")))
  for(y in yy[1]:yy[2]){dy<-d[d$year==y,];x<-log((dy$cases+0.5)/dy$E)[match(G$region[OBSR],dy$fbd_region)];mm<-moran_obs(x,G$region[OBSR])
    MOR<-rbind(MOR,data.frame(window=wn,year=as.character(y),quantity="observed_logSIR_year",moran_I=round(unname(mm$statistic),4),p_perm=round(mm$p.value,4)))}
  x<-log(ag$cases/ag$E)[match(G$region[OBSR],ag$fbd_region)];mm<-moran_obs(x,G$region[OBSR])
  MOR<-rbind(MOR,data.frame(window=wn,year=sprintf("%d-%d",yy[1],yy[2]),quantity="observed_logSIR_window",moran_I=round(unname(mm$statistic),4),p_perm=round(mm$p.value,4)))
  x<-log(rr$RR[OBSR]);mm<-moran_obs(x,G$region[OBSR])
  MOR<-rbind(MOR,data.frame(window=wn,year=sprintf("%d-%d",yy[1],yy[2]),quantity="bym2_logRR_window",moran_I=round(unname(mm$statistic),4),p_perm=round(mm$p.value,4)))
  LOG("%s: hot %d cold %d RR %.3f-%.3f phi %.3f · %.1f min",wn,sum(rr$class=="hotspot"),sum(rr$class=="coldspot"),min(rr$RR),max(rr$RR),FIT$phi[nrow(FIT)],as.numeric(difftime(Sys.time(),t0,units="mins")))
  rm(m)}
# 창 간 일치: 핫스팟 유지·신규·소멸
w<-names(WIN);H<-sapply(w,function(k)RRL$class[RRL$window==k]=="hotspot");C<-sapply(w,function(k)RRL$class[RRL$window==k]=="coldspot")
FIT$hot_kept_from_prev<-c(NA,sapply(2:length(w),function(k)sum(H[,k]&H[,k-1])));FIT$hot_new_vs_prev<-c(NA,sapply(2:length(w),function(k)sum(H[,k]&!H[,k-1])));FIT$hot_lost_vs_prev<-c(NA,sapply(2:length(w),function(k)sum(!H[,k]&H[,k-1])))
FIT$hot_all_windows<-sum(rowSums(H)==length(w));FIT$cold_all_windows<-sum(rowSums(C)==length(w))
FIT$spearman_logRR_vs_prev<-c(NA,sapply(2:length(w),function(k)round(cor(log(RRL$RR[RRL$window==w[k]]),log(RRL$RR[RRL$window==w[k-1]]),method="spearman"),3)))
write.csv(RRL,F1,row.names=FALSE,fileEncoding="UTF-8");write.csv(FIT,F2,row.names=FALSE,fileEncoding="UTF-8");write.csv(MOR,F3,row.names=FALSE,fileEncoding="UTF-8")
## 지도(3패널 a·b·c): 공통 로그대칭 8계급 + 핫(P>0.95) 검은 외곽선 / 콜드(P<0.05) 남색 외곽선
RRB<-c(-Inf,0.67,0.8,0.9,1,1.1,1.25,1.5,Inf);RRLAB<-c("<0.67","0.67-0.8","0.8-0.9","0.9-1","1-1.1","1.1-1.25","1.25-1.5",">=1.5")
PAL8<-c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B")
PEPB<-c(-Inf,0.05,0.5,0.8,0.95,Inf);PEPL<-c("<0.05","0.05-0.5","0.5-0.8","0.8-0.95",">0.95");PAL5<-c("#2166AC","#D1E5F0","#FEE0B6","#FDAE61","#B2182B")
G$sido<-substr(as.character(G$sigungu_cd),1,2);SIDO<-suppressWarnings(G%>%group_by(sido)%>%summarise(.groups="drop"))
pl<-setNames(sprintf("(%s) %s",letters[seq_along(w)],sapply(w,function(k)sprintf("%d-%d",WIN[[k]][1],WIN[[k]][2]))),w)
SP<-merge(G[,c("region","geometry")],RRL,by.x="region",by.y="fbd_region");SP$panel<-factor(pl[SP$window],levels=pl)
SP$cls<-factor(cut(SP$RR,RRB,labels=RRLAB),levels=RRLAB);SP$pcls<-factor(cut(SP$P_gt1,PEPB,labels=PEPL),levels=PEPL)
thm<-theme_void(base_size=11)+theme(strip.text=element_text(face="bold",size=11,hjust=0),legend.title=element_text(face="bold",size=9),legend.text=element_text(size=8),plot.margin=margin(6,6,6,6))
ol<-list(geom_sf(data=SP[SP$class=="hotspot",],aes(color="Hotspot, P(RR > 1) > 0.95"),fill=NA,linewidth=0.35,inherit.aes=FALSE),
         geom_sf(data=SP[SP$class=="coldspot",],aes(color="Coldspot, P(RR > 1) < 0.05"),fill=NA,linewidth=0.35,inherit.aes=FALSE),
         scale_color_manual(name=NULL,values=c("Hotspot, P(RR > 1) > 0.95"="black","Coldspot, P(RR > 1) < 0.05"="#3182BD")))
lvl<-function(lv,col){d<-SP[seq_along(lv),];d[[col]]<-factor(lv,levels=lv);geom_sf(data=d,aes(fill=.data[[col]]),alpha=0,color=NA,show.legend=TRUE,inherit.aes=FALSE)}   # 빈 계급도 범례에 표시
g1<-ggplot(SP)+geom_sf(aes(fill=cls),color="grey70",linewidth=0.06)+lvl(RRLAB,"cls")+geom_sf(data=SIDO,fill=NA,color="grey35",linewidth=0.2,inherit.aes=FALSE)+ol+facet_wrap(~panel,nrow=1)+
  scale_fill_manual(values=setNames(PAL8,RRLAB),breaks=rev(RRLAB),limits=RRLAB,drop=FALSE,name="Spatial RR\n(BYM2, window)")+guides(fill=guide_legend(override.aes=list(alpha=1,color="grey40")),color=guide_legend(override.aes=list(fill=NA)))+thm
ggsave(file.path(FIG,sprintf("%s_F4b_window_spatialRR.png",TAG)),g1,width=13,height=4.9,dpi=300,bg="white")
g2<-ggplot(SP)+geom_sf(aes(fill=pcls),color="grey70",linewidth=0.06)+lvl(PEPL,"pcls")+geom_sf(data=SIDO,fill=NA,color="grey35",linewidth=0.2,inherit.aes=FALSE)+facet_wrap(~panel,nrow=1)+
  scale_fill_manual(values=setNames(PAL5,PEPL),breaks=rev(PEPL),limits=PEPL,drop=FALSE,name="P(RR > 1)")+guides(fill=guide_legend(override.aes=list(alpha=1,color="grey40")))+thm
ggsave(file.path(FIG,sprintf("%s_F4c_window_exceedance.png",TAG)),g2,width=13,height=4.9,dpi=300,bg="white")
cap<-list(F4b=list(file=sprintf("%s_F4b_window_spatialRR.png",TAG),caption=sprintf("District spatial relative risk by window (BYM2 total spatial effect, refitted in each window with the main model: %s likelihood, set %s, expected counts from %s). Black outline = hotspot (posterior P(RR > 1) > 0.95); blue outline = coldspot (P < 0.05). Common classes across panels.",MAINFAM,BEST,if(SUF=="_std")"indirect standardization using female age distributions (E_std)" else "crude national rate of the year")),
          F4c=list(file=sprintf("%s_F4c_window_exceedance.png",TAG),caption="Posterior exceedance probability P(RR > 1) by district and window, same models as F4b."))
writeLines(enc2utf8(as.character(toJSON(cap,auto_unbox=TRUE,pretty=TRUE))),file.path(FIG,sprintf("%s_F4bc_captions.json",TAG)),useBytes=TRUE)
LOG("DONE → %s, %s, %s",basename(F1),basename(F2),basename(F3))
