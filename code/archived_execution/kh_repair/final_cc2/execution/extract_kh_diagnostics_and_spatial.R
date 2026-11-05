suppressPackageStartupMessages({library(INLA);library(jsonlite)})
D<-normalizePath('BRCA_재투고준비_20261104/kh_repair/final_cc2');RUN<-file.path(D,'execution');meta<-fromJSON(file.path(D,'dependency_manifest.json'))
B<-dirname(meta$panel_source);R<-file.path(B,'results_std');n<-229L;regions<-fromJSON(file.path(D,'region_order.json'))
stopifnot(digest::digest(file=meta$panel_source,algo='sha256')==meta$authorized_panel_sha256,
 digest::digest(file=file.path(R,'nb_BRCA.graph'),algo='sha256')==meta$graph_sha256)
expected<-read.csv(meta$panel_source,check.names=FALSE)
expected<-expected[order(expected$year,match(expected$fbd_region,regions)),]
cfg<-fromJSON(file.path(B,'config_BRCA_std.json'))
comparecols<-c('year','fbd_region','cases','E',unique(c(cfg$sets_after_vif[[cfg$rules$main_set]],cfg$forced)))
diags<-list()
for(ty in c('II','III','IV')){
 f<-file.path(R,paste0('S43_fit_KH_type',ty,'.rds'))
 if(!file.exists(f)){diags[[ty]]<-list(status='NO_SAVED_FIT');next}
 m<-readRDS(f);audit<-attr(m,'KH_execution_audit')
 stopifnot(identical(audit$panel_sha256,meta$authorized_panel_sha256),identical(audit$graph_sha256,meta$graph_sha256),audit$n_obs==3206,identical(as.numeric(audit$years),c(2011,2024)))
 dataequal<-all(vapply(comparecols,function(v)isTRUE(all.equal(m$.args$data[[v]],expected[[v]],check.attributes=FALSE,tolerance=0)),logical(1)))
 stopifnot(dataequal,identical(m$.args$data$EOFF,m$.args$data$E))
 pointwisefinite<-length(m$waic$local.waic)==3206&&all(is.finite(m$waic$local.waic))
 covpositive<-length(m$misc$cov.intern.eigenvalues)>0&&all(is.finite(m$misc$cov.intern.eigenvalues))&&min(m$misc$cov.intern.eigenvalues)>0
 fixedfinite<-all(is.finite(as.matrix(m$summary.fixed)));hyperfinite<-all(is.finite(as.matrix(m$summary.hyperpar)))
 statusok<-isTRUE(m$ok)&&isTRUE(m$mode$mode.status==0)&&isTRUE(m$misc$mode.status==0)
 valid<-statusok&&fixedfinite&&hyperfinite&&pointwisefinite&&covpositive&&is.finite(m$waic$waic)&&is.finite(m$dic$dic)&&nrow(m$summary.fitted.values)==3206
 means<-m$summary.random$st$mean;zm<-readRDS(file.path(D,paste0('KH_',ty,'_matrices.rds')))
 stopifnot(length(means)==3206,identical(as.integer(m$summary.random$st$ID),1:3206))
 cerr<-max(abs(zm$A%*%means));cerrnorm<-max(abs(as.numeric(zm$A%*%means))/sqrt(rowSums(zm$A^2)));valid<-valid&&cerrnorm<1e-6
 haswarn<-length(m$misc$warnings)>0
 diags[[ty]]<-list(status=if(valid&&haswarn)'COMPLETED_WITH_OPTIMIZATION_WARNING' else if(valid)'COMPLETED_DIAGNOSTICS_PASS' else 'REQUIRES_NUMERICAL_REVIEW',diagnostics_pass=valid,provisional=haswarn||!valid,ok=m$ok,mode_status=m$mode$mode.status,misc_mode_status=m$misc$mode.status,WAIC=m$waic$waic,DIC=m$dic$dic,exact_data_columns_equal_final_panel=dataequal,pointwise_waic_finite=pointwisefinite,hyperparameter_covariance_positive=covpositive,summary_fixed_finite=fixedfinite,summary_hyper_finite=hyperfinite,posterior_mean_constraint_max_abs=cerr,posterior_mean_constraint_max_rowL2_normalized=cerrnorm,constraint_diagnostic_tolerance=1e-6,fit_warnings=m$misc$warnings,stored_control=list(strategy=m$.args$control.inla$strategy,int_strategy=m$.args$control.inla$int.strategy,optimiser=m$.args$control.inla$optimiser,num_threads=m$.args$num.threads,diagonal=m$.args$control.inla$diagonal),cpu=m$cpu.used,execution=audit,rds_sha256=digest::digest(file=f,algo='sha256'),cpo='not computed in original88k specification; no CPO claim',convergence_scope='INLA optimization status and finite summaries/constraint checks; not MCMC convergence or global-optimum proof')
}
write_json(diags,file.path(RUN,'kh_diagnostics.json'),auto_unbox=TRUE,pretty=TRUE,digits=16,null='null',na='null')
if(!isTRUE(diags$IV$diagnostics_pass))stop('KH IV unavailable/not cleared; no spatial extraction')
f<-file.path(R,'S46_KH4_spatialRR.csv');if(file.exists(f))stop('Refuse overwrite existing S46')
mapone<-function(m){sr<-m$summary.random$area;stopifnot(nrow(sr)>=n,identical(as.integer(sr$ID[1:n]),seq_len(n)))
 p<-vapply(seq_len(n),function(i)1-inla.pmarginal(0,m$marginals.random$area[[i]]),numeric(1))
 data.frame(area=1:n,fbd_region=regions,RR=exp(sr$mean[1:n]),lo=exp(sr$`0.025quant`[1:n]),hi=exp(sr$`0.975quant`[1:n]),P_gt1=p,class=ifelse(p>.95,'hotspot',ifelse(p<.05,'coldspot','neutral')))}
kh<-readRDS(file.path(R,'S43_fit_KH_typeIV.rds'));out<-mapone(kh)
s6<-read.csv(file.path(R,'S6_spatial_RR.csv'),colClasses=c(sigungu_cd='character'))
stopifnot(identical(as.character(s6$fbd_region),regions),identical(as.integer(s6$area),seq_len(n)))
out$sigungu_cd<-s6$sigungu_cd
out$provisional<-isTRUE(diags$IV$provisional)
main<-mapone(readRDS(file.path(R,'refit_main.rds')));stopifnot(identical(main$class,s6$class))
write.csv(main,file.path(RUN,'main_spatial_full_precision.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(out,f,row.names=FALSE,fileEncoding='UTF-8')
write_json(list(provisional=diags$IV$provisional,numerical_warning=diags$IV$fit_warnings,source='S43_fit_KH_typeIV.rds reused; NO duplicate fit',graph_sha256=meta$graph_sha256,panel_sha256=meta$authorized_panel_sha256,RR='exp(posterior mean of log total BYM2 area effect), excludes temporal and interaction components',years='2011-2024 full-panel fit',hotspot='P(area effect >0)>0.95',coldspot='P(area effect >0)<0.05',spatial_file_sha256=digest::digest(file=f,algo='sha256')),file.path(RUN,'S46_extraction_provenance.json'),auto_unbox=TRUE,pretty=TRUE)
cat('S46 extracted from saved KH IV; no duplicate model fit\n')
