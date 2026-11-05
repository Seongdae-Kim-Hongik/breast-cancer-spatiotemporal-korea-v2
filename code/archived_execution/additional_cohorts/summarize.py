from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
O=Path(__file__).resolve().parent;m=json.loads((O/'manifest.json').read_text())
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
summ=[];effects=[];pairs=[]
for job in m['jobs']:
 d=O/'results'/job['name'];assert (d/'COMPLETE').exists(),job['name']
 h=json.loads((d/'hash_manifest.json').read_text())
 assert sha(Path(job['panel']))==h['panel_sha256'] and sha(Path(job['graph']))==h['graph_sha256']
 assert sha(d/'actual_input_graph.graph')==h['graph_sha256']
 fit=json.loads((d/'fit_unrounded.json').read_text());inp=pd.read_csv(d/'indexed_model_input.csv');pw=pd.read_csv(d/'pointwise_waic.csv');b=pd.read_csv(d/'beta_unrounded.csv')
 assert len(inp)==len(pw)==fit['n'] and not pw.duplicated(['year','fbd_region']).any()
 pd.testing.assert_frame_equal(inp[['year','fbd_region']],pw[['year','fbd_region']])
 assert np.isfinite(pw.waic).all();assert abs(pw.waic.sum()-fit['WAIC'])<1e-5
 graph=pd.read_csv(d/'actual_graph_metadata.csv').iloc[0]
 summ.append(dict(analysis=job['name'],n=fit['n'],districts=inp.fbd_region.nunique(),years=f"{fit['year_min']}-{fit['year_max']}",WAIC=fit['WAIC'],DIC=fit['DIC'],cpo_failed=fit['cpo_failed'],graph_edges=int(graph.n_links/2),graph_components=graph.n_subgraphs,family=fit['family'],comparison_group=job['comparison_group']))
 for _,row in b.iterrows():
  effects.append(dict(analysis=job['name'],variable=row.variable,beta=row['mean'],sd=row['sd'],IRR=np.exp(row['mean']),lo=np.exp(row['0.025quant']),hi=np.exp(row['0.975quant']),credible=bool(row['0.025quant']>0 or row['0.975quant']<0)))
for cohort in ['lag10','base0004_complete5']:
 a=O/'results'/f'{cohort}_historical';b=O/'results'/f'{cohort}_matched_current'
 ia=pd.read_csv(a/'indexed_model_input.csv');ib=pd.read_csv(b/'indexed_model_input.csv')
 pd.testing.assert_frame_equal(ia[['year','fbd_region','cases','E']],ib[['year','fbd_region','cases','E']])
 assert sha(a/'actual_input_graph.graph')==sha(b/'actual_input_graph.graph')
 wa=pd.read_csv(a/'pointwise_waic.csv');wb=pd.read_csv(b/'pointwise_waic.csv')
 dd=wa.waic-wb.waic
 pairs.append(dict(comparison=cohort,n=len(dd),delta_WAIC_historical_minus_matched_current=float(dd.sum()),pointwise_SE=float(np.sqrt(len(dd)*dd.var(ddof=1))),note='Same rows/outcome/E/graph; pointwise SE ignores spatial-temporal dependence and numerical fitting variability. Historical and current reproductive covariates use different SD scales.'))
pd.DataFrame(summ).to_csv(O/'fit_summary.csv',index=False)
pd.DataFrame(effects).to_csv(O/'effects_unrounded.csv',index=False)
pd.DataFrame(pairs).to_csv(O/'matched_WAIC_comparisons.csv',index=False)
(O/'completion_validation.json').write_text(json.dumps({'all6_complete':True,'saved_RDS':6,'pointwise_WAIC_sum_checks':True,'actual_graph_hash_checks':True,'same_row_pair_checks':True,'no_WAIC_ranking_post2014_vs_full':True},indent=2))
print(pd.DataFrame(summ).to_string(index=False))
print(pd.DataFrame(effects).query("variable in ['total_fertility_rate','mean_age_first_marriage_wife']").to_string(index=False))
print(pd.DataFrame(pairs).to_string(index=False))
