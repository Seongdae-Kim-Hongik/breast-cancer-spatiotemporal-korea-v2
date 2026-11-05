from pathlib import Path
W=Path(__file__).resolve().parent.parent;R=W/'refit_final_20261104/01_스크립트'
p=R/'85_window_maps.R';s=p.read_text();old='g<-inla.read.graph(GF);stopifnot(g$n==n,all(g$nnbs==card(nb)))';new='g<-inla.read.graph(GF);stopifnot(g$n==n)\nfor(i in seq_len(n))nb[[i]]<-as.integer(sort(g$nbs[[i]]))\nattr(nb,"sym")<-TRUE\nstopifnot(all(g$nnbs==card(nb)))';assert old in s;s=s.replace(old,new);p.write_text(s)
p=R/'88g_sens_suite.R';s=p.read_text();old='msign<-setNames(sign(log(M3$IRR)),M3$cov)';new='MF<-read.csv(file.path(RES,"refit_main_fixed_unrounded.csv"));msign<-setNames(sign(MF$beta),MF$cov)';assert old in s;s=s.replace(old,new)
s=s.replace('IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),P_gt1=round(pg,3)','IRR=exp(r$mean),lo=exp(r$`0.025quant`),hi=exp(r$`0.975quant`),P_gt1=pg,beta=r$mean')
s=s.replace('  m<-tryCatch(inla(form,','  strategy<-"default"\n  m<-tryCatch(inla(form,',1)
s=s.replace('  if(is.null(m))m<-tryCatch(inla(form,','  if(is.null(m)){strategy<-"eb_gaussian_retry";m<-tryCatch(inla(form,',1)
s=s.replace('message("robust fail: ",conditionMessage(e));NULL})\n  m}', 'message("robust fail: ",conditionMessage(e));NULL})}\n  if(!is.null(m))attr(m,"numerical_strategy")<-strategy\n  m}',1)
s=s.replace('if(is.null(m)){LOG("%s FAILED",an);return(invisible())}', 'if(is.null(m)){writeLines(toJSON(list(analysis=an,success=FALSE,reason="Both default and EB/Gaussian attempts failed; see group log"),auto_unbox=TRUE),file.path(RES,paste0("S41_status_",an,".json")));LOG("%s FAILED",an);return(invisible())}')
s=s.replace('  r<-do.call(rbind,lapply(intersect(vars,rownames(m$summary.fixed))', '  writeLines(toJSON(list(analysis=an,success=TRUE,strategy=attr(m,"numerical_strategy"),n_obs=nrow(data)),auto_unbox=TRUE),file.path(RES,paste0("S41_status_",an,".json")))\n  r<-do.call(rbind,lapply(intersect(vars,rownames(m$summary.fixed))')
p.write_text(s)
print('Patched85 canonical graph and88g precision/retry audit only. Spatial+ proposal still awaiting deployment.')
