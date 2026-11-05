from pathlib import Path
import json,hashlib,shutil,numpy as np,pandas as pd
from scipy.stats import spearmanr,pearsonr
N=Path(__file__).resolve().parent;D=N.parent;O=N/'adopted';dep=json.loads((D/'dependency_manifest.json').read_text());dec=json.loads((N/'KH_IV_decision.json').read_text());assert dec['adopt']
B=Path(dep['panel_source']).parent;R=B/'results_std';ROOT=B.parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(Path(dep['panel_source']))==dep['authorized_panel_sha256'] and sha(R/'nb_BRCA.graph')==dep['graph_sha256']
# Guard concurrent changes; initial result files are archived and never replaced.
initial=json.loads((N/'input_manifest.json').read_text())['initial_sha256']
for name,h in initial.items():assert sha(R/name)==h,name
kh=pd.read_csv(O/'S46_KH4_spatialRR.csv');main=pd.read_csv(D/'execution/main_spatial_full_precision.csv');assert len(kh)==len(main)==229
p=kh.merge(main,on='fbd_region',suffixes=('_kh','_main'),validate='one_to_one');assert len(p)==229
hk=set(kh.loc[kh['class']=='hotspot','fbd_region']);hc=set(kh.loc[kh['class']=='coldspot','fbd_region']);hm=set(main.loc[main['class']=='hotspot','fbd_region'])
col='canc__암_유방__연령표준화발생률';source=ROOT/'02_패널/MASTER_시군구_연도_wide.parquet';m=pd.read_parquet(source,columns=['year','fbd_region',col]);q=m[m.year.isin([2013,2018])].pivot(index='fbd_region',columns='year',values=col).reindex(kh.fbd_region)
a=q.dropna(subset=[2013,2018]).copy();assert len(a)==227;a['asir_mean_both_periods']=(a[2013]+a[2018])/2
aligned=kh.merge(a.rename(columns={2013:'asir_2009_2013',2018:'asir_2014_2018'}).reset_index(),on='fbd_region',validate='one_to_one');assert len(aligned)==227
missing=q[q.isna().any(axis=1)];assert set(missing.index)=={'경상남도창원시','세종시세종시'}
row=dict(n_hot=len(hk),n_cold=len(hc),hot_capital=sum(x.startswith(('서울','경기','인천')) for x in hk),jaccard_hot=len(hk&hm)/len(hk|hm),rho_RR_main=float(spearmanr(p.RR_kh,p.RR_main).statistic),rho_incidence=float(spearmanr(aligned.RR,aligned.asir_mean_both_periods).statistic),n_RR_main=229,n_registry=227,registry_periods='2009-2013 and 2014-2018',registry_average='unweighted arithmetic mean; BOTH periods observed',registry_exclusions='경상남도창원시;세종시세종시',pearson_log_incidence=float(pearsonr(np.log(aligned.RR),np.log(aligned.asir_mean_both_periods)).statistic),RR_definition='exp posterior mean log BYM2 area effect; full2011-2024; validated retry_A',provisional=False,numerical_status='VALIDATED_TWO_START_REPRODUCIBILITY',numerical_warning='none in either bounded retry',selected_fit='numerical_review/adopted/S43_fit_KH_typeIV_validated.rds')
pd.DataFrame([row]).to_csv(O/'S46_KH4_spatial_summary.csv',index=False);aligned.to_csv(O/'S46_registry_alignment_227.csv',index=False);missing.to_csv(O/'S46_registry_missing_periods.csv')
# Authorized S46 refresh only; preserve both initial S43 files byte-for-byte.
for name in ['S46_KH4_spatialRR.csv','S46_KH4_spatial_summary.csv']:shutil.copy2(O/name,R/name)
for name in ['S43_fit_KH_typeIV.rds','S43_part_KH_typeIV.csv']:assert sha(R/name)==initial[name]
selection={'KH_IV_inference_status':'ADOPT_VALIDATED_A','selected_fit':str((O/'S43_fit_KH_typeIV_validated.rds').resolve()),'selected_effects':str((O/'S43_part_KH_typeIV_validated.csv').resolve()),'selected_fit_sha256':sha(O/'S43_fit_KH_typeIV_validated.rds'),'excluded_initial_fit':str(R/'S43_fit_KH_typeIV.rds'),'excluded_initial_effects':str(R/'S43_part_KH_typeIV.csv'),'originals_preserved':True,'S46_refreshed':True,'S46_summary':row,'S46_hashes':{name:sha(R/name) for name in ['S46_KH4_spatialRR.csv','S46_KH4_spatial_summary.csv']},'panel_sha256':dep['authorized_panel_sha256'],'graph_sha256':dep['graph_sha256'],'registry_source_sha256':sha(source),'acceptance_decision':str((N/'KH_IV_decision.json').resolve()),'warning':'Initial S43 remains provisional and must not be selected by filename glob. Use the explicit validated selection. Registry227 includes Jeju, excludes Changwon/Sejong; it is not the227-node graph component.'}
(N/'inference_selection.json').write_text(json.dumps(selection,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
print(json.dumps(row,ensure_ascii=False))
