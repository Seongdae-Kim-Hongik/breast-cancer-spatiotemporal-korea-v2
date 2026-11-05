from pathlib import Path
import json,hashlib,numpy as np,pandas as pd
from scipy.stats import spearmanr,pearsonr
RUN=Path(__file__).resolve().parent;D=RUN.parent;dep=json.loads((D/'dependency_manifest.json').read_text());B=Path(dep['panel_source']).parent;R=B/'results_std';ROOT=B.parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(Path(dep['panel_source']))==dep['authorized_panel_sha256']
k=pd.read_csv(R/'S46_KH4_spatialRR.csv');main=pd.read_csv(RUN/'main_spatial_full_precision.csv')
assert len(k)==len(main)==229 and k.fbd_region.is_unique and main.fbd_region.is_unique
p=k.merge(main,on='fbd_region',suffixes=('_kh','_main'),validate='one_to_one');assert len(p)==229
hk=set(k.loc[k['class']=='hotspot','fbd_region']);hc=set(k.loc[k['class']=='coldspot','fbd_region']);hm=set(main.loc[main['class']=='hotspot','fbd_region']);union=hk|hm
col='canc__암_유방__연령표준화발생률';source=ROOT/'02_패널/MASTER_시군구_연도_wide.parquet'
m=pd.read_parquet(source,columns=['year','fbd_region',col]);q=m[m.year.isin([2013,2018])].pivot(index='fbd_region',columns='year',values=col)
q=q.reindex(k.fbd_region);missing=q[q[[2013,2018]].isna().any(axis=1)].copy();a=q.dropna(subset=[2013,2018]).copy();assert len(a)==227
assert set(missing.index)=={'경상남도창원시','세종시세종시'}
a['asir_mean_both_periods']=(a[2013]+a[2018])/2
aligned=k.merge(a.rename(columns={2013:'asir_2009_2013',2018:'asir_2014_2018'}).reset_index(),on='fbd_region',validate='one_to_one');assert len(aligned)==227 and aligned.RR.gt(0).all() and aligned.asir_mean_both_periods.gt(0).all()
rho_main=float(spearmanr(p.RR_kh,p.RR_main).statistic);rho_inc=float(spearmanr(aligned.RR,aligned.asir_mean_both_periods).statistic)
diagnostics=json.loads((RUN/'kh_diagnostics.json').read_text())['IV']
row=dict(provisional=diagnostics['provisional'],numerical_status=diagnostics['status'],numerical_warning=diagnostics['fit_warnings'],n_hot=len(hk),n_cold=len(hc),hot_capital=sum(x.startswith(('서울','경기','인천')) for x in hk),jaccard_hot=len(hk&hm)/len(union) if union else np.nan,rho_RR_main=rho_main,rho_incidence=rho_inc,n_RR_main=229,n_registry=227,registry_periods='2009-2013 and 2014-2018',registry_average='unweighted arithmetic mean; both periods observed required',registry_exclusions='경상남도창원시;세종시세종시',pearson_log_incidence=float(pearsonr(np.log(aligned.RR),np.log(aligned.asir_mean_both_periods)).statistic),RR_definition='exp posterior mean log BYM2 area effect from full2011-2024 KH IV; not total spatiotemporal risk',hot_capital_definition='count of hotspots whose region starts 서울/경기/인천',same_data_graph=True)
output=R/'S46_KH4_spatial_summary.csv';assert not output.exists();pd.DataFrame([row]).to_csv(output,index=False)
aligned.to_csv(RUN/'S46_registry_alignment_227.csv',index=False);missing.to_csv(RUN/'S46_registry_missing_periods.csv')
meta=dict(summary=row,panel_sha256=dep['authorized_panel_sha256'],graph_sha256=dep['graph_sha256'],registry_source_sha256=sha(source),fit_sha256=sha(R/'S43_fit_KH_typeIV.rds'),spatial_rr_sha256=sha(R/'S46_KH4_spatialRR.csv'),summary_sha256=sha(output),warning='Registry n227 excludes Changwon/Sejong and INCLUDES Jeju; it is NOT the227-node graph component. Periods differ; this is ecological spatial concordance, not incidence prediction or external causal validation.')
(RUN/'S46_summary_provenance.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2,allow_nan=False)+'\n');print(json.dumps(row,ensure_ascii=False))
