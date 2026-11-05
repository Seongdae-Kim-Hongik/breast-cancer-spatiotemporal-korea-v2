suppressPackageStartupMessages({library(INLA);library(Matrix);library(jsonlite);library(sf);library(spdep)})
D<-normalizePath('BRCA_재투고준비_20261104/kh_repair/final_cc2');N<-file.path(D,'numerical_review');dep<-fromJSON(file.path(D,'dependency_manifest.json'));source(file.path(D,'kh_constraints.R'))
B<-dirname(dep$panel_source);R<-file.path(B,'results_std');g<-inla.read.graph(file.path(D,'chosen_nb_BRCA.graph'));regions<-fromJSON(file.path(D,'region_order.json'));n<-g$n
kh_assert_dependencies(g,file.path(D,'chosen_nb_BRCA.graph'),regions,dep$panel_source,D)
S<-fromJSON(file.path(B,'config_BRCA_std.json'));HB<-readLines(file.path(B,'HM_best_BRCA_std.txt'));kv<-function(k)sub(paste0('^',k,'='),'',grep(paste0('^',k,'='),HB,value=TRUE)[1])
MAINFAM<-kv('main_family');BEST<-kv('best_set');stopifnot(MAINFAM=='nbinomial',BEST=='S1_narrative_full')
FORCED<-S$forced;HYPV<-S$sets_after_vif[[BEST]];MAINV<-unique(c(HYPV,FORCED));stopifnot(length(MAINV)==24)
DD<-read.csv(dep$panel_source,check.names=FALSE);DD<-DD[order(DD$year,match(DD$fbd_region,regions)),];DD$area<-match(DD$fbd_region,regions);DD$time<-DD$year-2010;DD$yiid<-DD$time;DD$st<-(DD$time-1)*n+DD$area;DD$EOFF<-DD$E
kh_check_panel(DD,n,MAINV)
pcU<-function(u=.5,a=.01,pu=.5,pa=2/3)list(prec=list(prior='pc.prec',param=c(u,a)),phi=list(prior='pc',param=c(pu,pa)))
pc<-pcU();pcp<-list(prec=list(prior='pc.prec',param=c(.5,.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=.001)
TI<-"f(time,model='rw1',scale.model=TRUE,hyper=pcp)";IT<-"f(yiid,model='iid',hyper=pcp)"
SPf<-function(gr='g',hyper='pc')sprintf("f(area,model='bym2',graph=%s,scale.model=TRUE,constr=TRUE,hyper=%s)",gr,hyper)
formula_for<-function(sp=SPf(),extra=character(0))as.formula(paste('cases ~',paste(c('1',sprintf('`%s`',MAINV),extra,sp[nzchar(sp)],TI,IT),collapse='+')),env=.GlobalEnv)
inla.setOption(num.threads='1:1',inla.timeout=900)
LOG<-function(...)cat(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),sprintf(...),'\n')
fit_once<-function(id,outdir,sp=SPf(),extra=character(0),ci=list(control.vb=list(emergency=50)),cm=NULL,graphnote='canonical_contact593',A=NULL,graph_file=file.path(D,'chosen_nb_BRCA.graph')){
 target<-file.path(outdir,paste0(id,'.rds'));stopifnot(!file.exists(target));started<-as.character(Sys.time());warn<-character(0);err<-NULL
 LOG('START %s',id)
 m<-tryCatch(withCallingHandlers(inla(formula_for(sp,extra),family=MAINFAM,data=DD,E=EOFF,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE,cpo=FALSE,config=FALSE),control.inla=ci,control.mode=cm,num.threads='1:1',verbose=FALSE),warning=function(w){warn<<-c(warn,conditionMessage(w));message(conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){err<<-conditionMessage(e);NULL})
 audit<-list(id=id,start=started,end=as.character(Sys.time()),error=err,R_warnings=warn,panel_sha256=dep$authorized_panel_sha256,graph_sha256=dep$graph_sha256,graph_policy=graphnote,model_graph_sha256=digest::digest(file=graph_file,algo='sha256'),years=c(2011,2024),n_obs=3206,control_inla=ci,control_mode=if(is.null(cm))NULL else list(theta=cm$theta,restart=cm$restart,fixed=cm$fixed,latent_initial_supplied=!is.null(cm$x)),num_threads='1:1')
 if(!is.null(m)){
  attr(m,'KH_execution_audit')<-audit
  same<-all(vapply(c('year','fbd_region','cases','E',MAINV),function(v)isTRUE(all.equal(m$.args$data[[v]],DD[[v]],check.attributes=FALSE,tolerance=0)),logical(1)))
  audit$data_exact<-same;audit$ok<-m$ok;audit$mode_status<-m$mode$mode.status;audit$misc_mode_status<-m$misc$mode.status;audit$warnings<-m$misc$warnings
  audit$fixed_finite<-all(is.finite(as.matrix(m$summary.fixed)));audit$hyper_finite<-all(is.finite(as.matrix(m$summary.hyperpar)));audit$pointwise_finite<-length(m$waic$local.waic)==3206&&all(is.finite(m$waic$local.waic));audit$cov_eigen_min<-min(m$misc$cov.intern.eigenvalues);audit$WAIC<-m$waic$waic;audit$DIC<-m$dic$dic;audit$log_posterior_mode<-m$mode$log.posterior.mode
  if(!is.null(A))audit$constraint_error<-max(abs(as.numeric(A%*%m$summary.random$st$mean))/sqrt(rowSums(A^2)))
  audit$diagnostics_pass<-same&&isTRUE(m$ok)&&m$mode$mode.status==0&&m$misc$mode.status==0&&audit$fixed_finite&&audit$hyper_finite&&audit$pointwise_finite&&is.finite(audit$cov_eigen_min)&&audit$cov_eigen_min>0&&(is.null(A)||audit$constraint_error<1e-6)
  audit$clean<-audit$diagnostics_pass&&length(m$misc$warnings)==0&&length(warn)==0
  attr(m,'numerical_review')<-audit;saveRDS(m,target)
  f<-m$summary.fixed;write.csv(data.frame(cov=gsub('`','',rownames(f)),f,IRR=exp(f$mean),lo=exp(f$`0.025quant`),hi=exp(f$`0.975quant`),check.names=FALSE),file.path(outdir,paste0(id,'_fixed.csv')),row.names=FALSE)
  write.csv(data.frame(parameter=rownames(m$summary.hyperpar),m$summary.hyperpar,check.names=FALSE),file.path(outdir,paste0(id,'_hyper.csv')),row.names=FALSE)
  writeLines(m$logfile,file.path(outdir,paste0(id,'_internal.log')))
  LOG('END %s WAIC=%.9f status=%s warnings=%s clean=%s',id,m$waic$waic,m$mode$mode.status,paste(m$misc$warnings,collapse=';'),audit$clean)
 }else LOG('ERROR %s %s',id,err)
 write_json(audit,file.path(outdir,paste0(id,'_diagnostics.json')),auto_unbox=TRUE,pretty=TRUE,digits=16,null='null',na='null')
 invisible(m)
}
