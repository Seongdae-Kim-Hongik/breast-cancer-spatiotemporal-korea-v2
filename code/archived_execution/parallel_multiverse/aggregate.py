from pathlib import Path
import json,sys,hashlib
import numpy as np,pandas as pd
root=Path(sys.argv[1]).resolve();b=root/'04_분석/BRCA';r=b/'results_std';workers=root/'multiverse_workers'
specs=pd.read_csv(r/'S7_multiverse_specs.csv');meta=pd.read_csv(b/'panel_META_BRCA_std.csv').set_index('variable');cfg=json.loads((b/'config_BRCA_std.json').read_text())
main=pd.read_csv(r/'refit_main_fixed_unrounded.csv').set_index('cov');hyp=cfg['sets_after_vif'][cfg['rules']['main_set']]
status=[];frames=[]
for s in specs.itertuples():
 found=list(workers.glob(f'worker_*/spec_{s.spec:03d}_status.json'));assert len(found)==1,(s.spec,'Missing/duplicate status')
 st=json.loads(found[0].read_text());status.append(st)
 if not st['success']:continue
 p=found[0].with_name(f'spec_{s.spec:03d}.csv');x=pd.read_csv(p)
 assert set(x['cov'])==set(s.vars.split('|')) and len(x)==len(s.vars.split('|'))
 assert x.spec.eq(s.spec).all() and np.isfinite(x[['IRR','lo','hi','beta']]).all().all()
 frames.append(x)
assert frames,'No successful fits'
mv=pd.concat(frames,ignore_index=True);summary=[];pool=sorted(set('|'.join(specs.vars).split('|')))
for v in pool:
 x=mv[mv['cov']==v];planned=sum(v in s.split('|') for s in specs.vars)
 if x.empty: summary.append(dict(cov=v,n_specs=0,n_specs_planned_for_cov=planned,robust=False,provisional=True));continue
 ref=np.sign(main.loc[v,'beta']) if v in main.index else np.sign(np.median(x.beta))
 cf=float(x.cred.mean());sa=float((np.sign(x.beta)==ref).mean())
 summary.append(dict(cov=v,label_en=meta.loc[v,'label_en'],label_ko=meta.loc[v,'label_ko'],in_main_set=v in hyp,n_specs=len(x),n_specs_planned_for_cov=planned,cred_frac=cf,sign_agree=sa,IRR_median=x.IRR.median(),IRR_min=x.IRR.min(),IRR_max=x.IRR.max(),robust=cf>=.75 and sa>=.95,provisional=len(x)<50,n_specs_total=len(specs),n_specs_successful=sum(z['success'] for z in status)))
for f in ['S7_multiverse_results.csv','S7_multiverse_summary.csv']:assert not (r/f).exists(),'Refuse overwrite'
def write_r_csv(frame,path):
 frame=frame.copy()
 for c in frame.select_dtypes(include=['bool']).columns:frame[c]=frame[c].map({True:'TRUE',False:'FALSE'})
 frame.to_csv(path,index=False)
write_r_csv(mv,r/'S7_multiverse_results.csv');write_r_csv(pd.DataFrame(summary),r/'S7_multiverse_summary.csv')
pd.DataFrame([{'spec':z['spec'],'status':'success' if z['success'] else 'failed','reason':str(z.get('errors',[])) if not z['success'] else '', 'approximation':z['integration']} for z in status]).to_csv(r/'S7_execution_status.csv',index=False)
report={'planned':len(specs),'success':sum(z['success'] for z in status),'failed_spec_ids':[z['spec'] for z in status if not z['success']],'retry_ids':[z['spec'] for z in status if z['integration']!='default'],'precision':'Unrounded beta determines sign; credible status computed before rounding','selection':'Original candidate pool, forced controls and 200 specification IDs retained','seeds':'20261104 + specification id','workers':4,'threads_per_worker':1,'panel_sha256':hashlib.sha256((b/'panel_BRCA_std.csv').read_bytes()).hexdigest(),'statuses':status}
(r/'S7_execution_audit.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='statuses'}))
