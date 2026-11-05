from pathlib import Path
import json,hashlib,shutil
import pandas as pd
import numpy as np
O=Path(__file__).resolve().parent;W=O.parent;F=W/'refit_final_20261104';B=O/'04_분석/BRCA'
for d in [B,O/'panels',O/'graphs',O/'scripts',O/'results']:d.mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
ready=json.loads((F/'DATA_READY.json').read_text())
for name,h in ready['sha256'].items():
 if name.startswith('04_분석/BRCA/panel') or name=='input_graph.graph':assert sha(F/name)==h,name
for name in ['panel_BRCA_std.csv','panel_META_BRCA_std.csv','config_BRCA_std.json','panel_INFO_BRCA_std.json']:
 shutil.copyfile(F/'04_분석/BRCA'/name,B/name)
shutil.copyfile(F/'input_graph.graph',O/'graphs/contact.graph')
original=W.parent.parent/'체육_시공간_데이터셋/04_분석/BRCA/results_std/nb_BRCA.graph'
shutil.copyfile(original,O/'graphs/original.graph')
p=pd.read_csv(B/'panel_BRCA_std.csv');meta=pd.read_csv(B/'panel_META_BRCA_std.csv').set_index('variable')
vs=['total_fertility_rate','mean_age_first_marriage_wife'];cs=[meta.loc[v,'source_column'] for v in vs]
master=F/'02_패널/MASTER_시군구_연도_wide.parquet'
m=pd.read_parquet(master,columns=['year','fbd_region']+cs)
assert not m.duplicated(['year','fbd_region']).any()
report={};source_rows=[];jobs=[]
for name in ['lag10','base0004_complete5']:
 x=p[['year','fbd_region']].copy()
 for v,c in zip(vs,cs):
  tab=m.pivot(index='fbd_region',columns='year',values=c)
  if name=='lag10':x[v]=[tab.at[r,y-10] if r in tab.index and y-10 in tab.columns else np.nan for r,y in zip(p.fbd_region,p.year)]
  else:x[v]=p.fbd_region.map(tab.reindex(columns=range(2000,2005)).dropna().mean(axis=1))
  for r,y in zip(p.fbd_region,p.year):
   for sy in ([y-10] if name=='lag10' else range(2000,2005)):
    source_rows.append(dict(cohort=name,year=y,fbd_region=r,variable=v,source_column=c,source_year=sy,value=tab.at[r,sy] if r in tab.index and sy in tab.columns else np.nan))
 keep=x[vs].notna().all(axis=1);hist=p.loc[keep].copy();curr=p.loc[keep].copy();scales={}
 for v in vs:
  val=x.loc[keep,v];sd=val.std(ddof=0);assert sd>0
  scales[v]={'mean':float(val.mean()),'sd_ddof0':float(sd)}
  hist[v]=(val-val.mean())/sd
  # historical raw columns must describe the actually fitted historical exposure.
  hist[v+'_raw']=val
 for label,d in [('historical',hist),('matched_current',curr)]:
  job=f'{name}_{label}';dest=O/'panels'/f'{job}.csv';d.to_csv(dest,index=False)
  jobs.append(dict(name=job,panel=str(dest),graph=str(O/'graphs/contact.graph'),comparison_group=name))
 pd.testing.assert_frame_equal(hist[['year','fbd_region','cases','E']],curr[['year','fbd_region','cases','E']])
 x.assign(included=keep).to_csv(O/'panels'/f'{name}_exposure_eligibility.csv',index=False)
 report[name]={'n_obs':len(hist),'n_districts':hist.fbd_region.nunique(),'excluded_obs':int((~keep).sum()),'historical_scaling':scales,'current_covariates_scale':'final full-panel SD retained','E':'Final full-panel E retained without renormalization','imputation':'None; exact t-10 or all5 years2000-04 observed','historical_raw_columns':'Updated to matching historical exposure; parent preparation leaves these columns current, only z enters model'}
for name,frame,graph in [('original_graph_main',p,'original'),('post2014_contact',p[p.year>=2014],'contact')]:
 dest=O/'panels'/f'{name}.csv';frame.to_csv(dest,index=False)
 jobs.append(dict(name=name,panel=str(dest),graph=str(O/f'graphs/{graph}.graph'),comparison_group='full_panel_graph' if name=='original_graph_main' else 'post2014_only'))
pd.DataFrame(source_rows).to_csv(O/'historical_source_rows.csv',index=False)
src=(W/'audit_fit.R').read_text()
start=src.index('ROOT<-');end=src.index('\nSHP<-',start)
src=src[:start]+'ROOT<-normalizePath(Sys.getenv("AUDIT_ROOT"),mustWork=TRUE)'+src[end:]
start=src.index('GF<-file.path');end=src.index('\n## ── 데이터 색인',start)
src=src[:start]+'''GF<-Sys.getenv("AUDIT_GRAPH");stopifnot(nzchar(GF),file.exists(GF))
g<-inla.read.graph(GF);stopifnot(g$n==nrow(G))
nb<-g$nbs;class(nb)<-"nb";attr(nb,"region.id")<-as.character(seq_len(g$n));attr(nb,"sym")<-TRUE
stopifnot(all(lengths(nb)>0))
for(i in seq_len(g$n))stopifnot(all(vapply(nb[[i]],function(j)i%in%nb[[j]],logical(1))))
file.copy(GF,file.path(RES,"actual_input_graph.graph"),overwrite=FALSE)
write.csv(data.frame(n_regions=g$n,n_links=sum(lengths(nb)),n_subgraphs=n.comp.nb(nb)$nc,graph_source=GF),file.path(RES,"actual_graph_metadata.csv"),row.names=FALSE)
write.csv(data.frame(area=seq_len(g$n),fbd_region=G$region),file.path(RES,"graph_node_order.csv"),row.names=FALSE)
''' +src[end:]
src=src.replace('if(nzchar(Sys.getenv("AUDIT_GRAPH")))g<-inla.read.graph(Sys.getenv("AUDIT_GRAPH"))','')
src=src.replace('m<-fitf(MAINV,"nbinomial",extra=form_extra)','''write.csv(DD,file.path(RES,"indexed_model_input.csv"),row.names=FALSE)
stopifnot(all(complete.cases(DD[,c(MAINV,"cases","E")])),all(DD$E>0))
X<-model.matrix(as.formula(paste("~",paste(bt(MAINV),collapse="+"))),DD)
stopifnot(qr(X)$rank==ncol(X))
m<-fitf(MAINV,"nbinomial",extra=form_extra)''')
src+='''
beta<-data.frame(variable=rownames(m$summary.fixed),m$summary.fixed,check.names=FALSE)
write.csv(beta,file.path(RES,"beta_unrounded.csv"),row.names=FALSE)
u<-data.frame(variable=rownames(m$summary.fixed),IRR=exp(m$summary.fixed$mean),lo=exp(m$summary.fixed$`0.025quant`),hi=exp(m$summary.fixed$`0.975quant`))
write.csv(u,file.path(RES,"IRR_unrounded.csv"),row.names=FALSE)
write_json(list(WAIC=m$waic$waic,DIC=m$dic$dic,n=nrow(DD),year_min=min(DD$year),year_max=max(DD$year),family="nbinomial",graph=GF,main_variables=MAINV,inla_version=as.character(packageVersion("INLA")),threads="1:1",fallback_used=FALSE,cpo_failed=sum(m$cpo$failure>0,na.rm=TRUE)),file.path(RES,"fit_unrounded.json"),auto_unbox=TRUE,pretty=TRUE,digits=NA)
stopifnot(length(m$waic$local.waic)==nrow(DD),all(is.finite(m$waic$local.waic)))
writeLines(capture.output(sessionInfo()),file.path(RES,"sessionInfo.txt"))
writeLines("complete",file.path(RES,"COMPLETE"))
'''
(O/'scripts/fit.R').write_text(src)
manifest={'final_root':str(F),'definitions_source':str(W/'corrected_extra_plan/prepare_strict_historical.py'),'parent_audit_source_sha256':sha(W/'audit_fit.R'),'master_source':str(master),'master_sha256':sha(master),'definitions_sha256':sha(W/'corrected_extra_plan/prepare_strict_historical.py'),'cohorts':report,'jobs':jobs,'input_sha256':{str(f.relative_to(O)):sha(f) for folder in ['04_분석','graphs','panels','scripts'] for f in (O/folder).rglob('*') if f.is_file()},'time_index':'min(year) of each exact panel =1; gaps retained by calendar year;post2014 starts2014','fitting':'NB, same main set/prior,default INLA,no automatic EB retry,one fit at a time','no_WAIC_ranking_across_different_rows':True}
(O/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
print(json.dumps(report,ensure_ascii=False,indent=2))
