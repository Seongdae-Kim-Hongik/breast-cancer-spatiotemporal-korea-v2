suppressPackageStartupMessages({library(INLA);library(jsonlite)})
D<-normalizePath('BRCA_재투고준비_20261104/kh_repair/final_cc2');N<-file.path(D,'numerical_review');crit<-fromJSON(file.path(N,'acceptance_criteria.json'))
paths<-c(initial=file.path(N,'initial_snapshot/S43_fit_KH_typeIV.rds'),A=file.path(N,'retry/retry_A.rds'),B=file.path(N,'retry/retry_B.rds'))
if(!all(file.exists(paths))){write_json(list(adopt=FALSE,reason='One or more bounded retries failed/no saved fit',initial_excluded=TRUE),file.path(N,'KH_IV_decision.json'),auto_unbox=TRUE,pretty=TRUE);quit(status=0)}
mods<-lapply(paths,readRDS)
compare<-function(a,b){
 x<-a$summary.fixed;y<-b$summary.fixed;stopifnot(identical(rownames(x),rownames(y)))
 sd<-(x$sd+y$sd)/2
 rrA<-exp(a$summary.random$area$mean[1:229]);rrB<-exp(b$summary.random$area$mean[1:229])
 pa<-vapply(1:229,function(i)1-inla.pmarginal(0,a$marginals.random$area[[i]]),numeric(1));pb<-vapply(1:229,function(i)1-inla.pmarginal(0,b$marginals.random$area[[i]]),numeric(1))
 cl<-function(p)ifelse(p>.95,'hot',ifelse(p<.05,'cold','neutral'))
 f<-data.frame(cov=rownames(x),mean_a=x$mean,mean_b=y$mean,mean_delta=y$mean-x$mean,delta_over_sd=abs(y$mean-x$mean)/sd,lo_a=x$`0.025quant`,lo_b=y$`0.025quant`,hi_a=x$`0.975quant`,hi_b=y$`0.975quant`)
 list(metrics=list(delta_WAIC=b$waic$waic-a$waic$waic,abs_delta_log_posterior_mode=abs(b$mode$log.posterior.mode-a$mode$log.posterior.mode),max_fixed_delta_over_sd=max(f$delta_over_sd),max_CrI_endpoint_delta_over_sd=max(abs(cbind(f$lo_b-f$lo_a,f$hi_b-f$hi_a))/sd),credible_classifications_match=identical(x$`0.025quant`>0|x$`0.975quant`<0,y$`0.025quant`>0|y$`0.975quant`<0),max_spatial_RR_relative_difference=max(abs(rrB/rrA-1)),spatial_class_agreement=mean(cl(pa)==cl(pb))),fixed=f)
}
AB<-compare(mods$A,mods$B);IA<-compare(mods$initial,mods$A);IB<-compare(mods$initial,mods$B)
write.csv(AB$fixed,file.path(N,'retry_AB_fixed_comparison.csv'),row.names=FALSE);write.csv(IA$fixed,file.path(N,'initial_vs_A_fixed_comparison.csv'),row.names=FALSE);write.csv(IB$fixed,file.path(N,'initial_vs_B_fixed_comparison.csv'),row.names=FALSE)
hyper<-do.call(rbind,lapply(names(mods),function(k)data.frame(fit=k,parameter=rownames(mods[[k]]$summary.hyperpar),mods[[k]]$summary.hyperpar,check.names=FALSE)))
write.csv(hyper,file.path(N,'hyperparameter_comparison.csv'),row.names=FALSE)
auditA<-attr(mods$A,'numerical_review');auditB<-attr(mods$B,'numerical_review');v<-AB$metrics
checks<-list(both_clean=isTRUE(auditA$clean)&&isTRUE(auditB$clean),WAIC=abs(v$delta_WAIC)<=crit$max_abs_delta_WAIC_between_retries,log_posterior=v$abs_delta_log_posterior_mode<=crit$max_abs_log_posterior_mode_difference,fixed=v$max_fixed_delta_over_sd<=crit$max_fixed_mean_delta_over_mean_posterior_SD,CrI=v$max_CrI_endpoint_delta_over_sd<=crit$max_CrI_endpoint_delta_over_mean_posterior_SD,credible=isTRUE(v$credible_classifications_match),spatial_RR=v$max_spatial_RR_relative_difference<=crit$max_spatial_RR_relative_difference,spatial_class=v$spatial_class_agreement>=crit$min_spatial_class_agreement)
# Explicitly access the actual metric name; no empty-vector acceptance.
checks$log_posterior<-v$abs_delta_log_posterior_mode<=crit$max_abs_log_posterior_mode_difference
stopifnot(all(lengths(checks)==1));adopt<-all(unlist(checks))
decision<-list(adopt=adopt,checks=checks,criteria=crit,comparison_AB=v,comparison_initial_A=IA$metrics,comparison_initial_B=IB$metrics,selected_fit=if(adopt)paths[['A']] else NULL,selection_rule='A warm-start selected a priori only if both clean and agreement criteria pass; NOT minimum WAIC',initial_fit_excluded=TRUE,final_inference=if(adopt)'Use validated retry_A with preserved original/retry_B provenance' else 'OMIT KH IV and its derived S46 from final inference',diagnostics_A=auditA,diagnostics_B=auditB)
write_json(decision,file.path(N,'KH_IV_decision.json'),auto_unbox=TRUE,pretty=TRUE,digits=16,null='null')
print(checks);print(v);cat('ADOPT=',adopt,'\n')
