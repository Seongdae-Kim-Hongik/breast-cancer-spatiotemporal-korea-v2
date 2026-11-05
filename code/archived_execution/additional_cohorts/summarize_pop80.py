from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
O=Path(__file__).resolve().parent;W=O.parent
M=json.loads((O/'pop80_manifest.json').read_text())
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
B=W/'refit_final_20261104/04_분석/BRCA';R=B/'results_std'
assert sha(B/'panel_BRCA_std.csv')==sha(O/'04_분석/BRCA/panel_BRCA_std.csv')
assert sha(W/'refit_final_20261104/input_graph.graph')==sha(O/'graphs/contact.graph')
baseline=pd.read_csv(R/'refit_main_pointwise_waic.csv')
fit0=pd.read_csv(R/'S3_main_fit.csv');assert fit0.family.iloc[0]=='nbinomial'
main=pd.read_csv(R/'S3_main_IRR.csv')['cov'].tolist()
rows=[];effects=[]
for job in M['jobs']:
 d=O/'results'/job['name'];assert (d/'COMPLETE').exists()
 h=json.loads((d/'hash_manifest.json').read_text())
 for f,v in h['output_sha256'].items():assert sha(d/f)==v
 fit=json.loads((d/'fit_unrounded.json').read_text())
 assert len(fit['main_variables'])==25 and set(fit['main_variables'])==set(main)|{'pop80_share'}
 actual=(d/'actual_MAINV.txt').read_text().splitlines();assert actual==fit['main_variables']
 inp=pd.read_csv(d/'indexed_model_input.csv');p=pd.read_csv(B/'panel_BRCA_std.csv')
 pd.testing.assert_frame_equal(inp[p.columns],p,check_exact=False,check_dtype=False,rtol=1e-12,atol=1e-12)
 assert len(inp)==3206 and sha(d/'actual_input_graph.graph')==sha(O/'graphs/contact.graph')
 pw=pd.read_csv(d/'pointwise_waic.csv');cmp=pw.merge(baseline,on=['year','fbd_region'],validate='one_to_one',suffixes=('_extra','_main'))
 assert len(cmp)==3206 and abs(pw.waic.sum()-fit['WAIC'])<1e-5
 delta=cmp.waic_extra-cmp.waic_main
 rows.append(dict(analysis=job['name'],role=M['definitions'][job['name']]['role'],n=fit['n'],WAIC=fit['WAIC'],baseline_WAIC_pointwise_sum=baseline.waic.sum(),delta_WAIC_vs_main=delta.sum(),pointwise_SE=np.sqrt(len(delta)*delta.var(ddof=1)),cpo_failed=fit['cpo_failed'],main_variables=len(main),extra_variables=1))
 b=pd.read_csv(d/'beta_unrounded.csv')
 for _,r in b.iterrows():effects.append(dict(analysis=job['name'],variable=r.variable,beta=r['mean'],sd=r['sd'],IRR=np.exp(r['mean']),lo=np.exp(r['0.025quant']),hi=np.exp(r['0.975quant']),credible=bool(r['0.025quant']>0 or r['0.975quant']<0)))
pd.DataFrame(rows).to_csv(O/'pop80_fit_summary.csv',index=False)
pd.DataFrame(effects).to_csv(O/'pop80_effects_unrounded.csv',index=False)
(O/'pop80_validation.json').write_text(json.dumps({'both_complete':True,'actual_MAINV_verified':25,'all24_main_variables_preserved':True,'same3206_rows_outcomes_E_main_covariates':True,'actual_graph_hash_verified':True,'full_output_hashes_verified':True,'baseline_sources':{str(f):sha(f) for f in [B/'panel_BRCA_std.csv',R/'refit_main_pointwise_waic.csv',R/'S3_main_fit.csv',R/'S3_main_IRR.csv']},'comparison_limit':'Descriptive pointwise WAIC differences only; nominal SE does not account for spatial-temporal dependence or numerical fitting variability'},indent=2))
print(pd.DataFrame(rows).to_string(index=False))
print(pd.DataFrame(effects).query("variable in ['pop80_share','total_fertility_rate','mean_age_first_marriage_wife']").to_string(index=False))
