from pathlib import Path
import shutil,json
W=Path(__file__).resolve().parent.parent
root=W/'refit_final_20261104';stage=W/'population_final_audit/final_stage'
assert not root.exists(),'Never overwrite a prior final run'
root.mkdir()
for name in ['02_패널','04_분석']:shutil.copytree(stage/name,root/name)
(root/'01_스크립트').mkdir()
s=(W/'refit_corrected/01_스크립트/81_fit_inla.R').read_text()
needle='stopifnot(all(card(nb)>0));nb2INLA(GF,nb);'
replacement='''if(nzchar(Sys.getenv("BRCA_GRAPH"))){
 supplied<-inla.read.graph(Sys.getenv("BRCA_GRAPH"));stopifnot(supplied$n==nrow(G))
 for(i in seq_along(nb))nb[[i]]<-as.integer(sort(supplied$nbs[[i]]))
 attr(nb,"sym")<-TRUE
}
stopifnot(all(card(nb)>0));nb2INLA(GF,nb);'''
assert needle in s;s=s.replace(needle,replacement)
marker='## ── S7 다중우주:'
assert marker in s
hook='''if(Sys.getenv("ST_ONLY_CORE")=="1"){
 main_saved<-readRDS(file.path(RES,"refit_main.rds"));rr<-main_saved$summary.fixed
 write.csv(data.frame(cov=rownames(rr),beta=rr$mean,sd=rr$sd,lo=rr$`0.025quant`,hi=rr$`0.975quant`),file.path(RES,"refit_main_fixed_unrounded.csv"),row.names=FALSE)
 LOG("CORE ONLY DONE — multiverse deferred to isolated workers")
 quit(save="no",status=0)
}

'''
s=s.replace(marker,hook+marker)
(root/'01_스크립트/81_fit_inla.R').write_text(s)
extras=W/'refit_corrected/extra_prepared/01_스크립트'
for p in extras.iterdir():
 if p.is_file() and p.name!='81_fit_inla.R':shutil.copy2(p,root/'01_스크립트'/p.name)
# Core input source needed only for sensitivity imputation masks, not old fits.
shutil.copy2(W/'refit_corrected/extra_prepared/02_패널/MASTER_시군구_연도_wide.parquet',root/'02_패널/MASTER_시군구_연도_wide.parquet')
(root/'NOT_YET_AUTHORIZED_TO_RUN.txt').write_text('Data directly verified; graph selection and complete integration gate pending. No fit outputs imported.\n')
print(root)
