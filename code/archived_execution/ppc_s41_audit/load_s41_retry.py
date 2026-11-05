"""Read-only adoption of independently verified isolated S41 numerical retries."""
from pathlib import Path
import json,hashlib
import pandas as pd
import numpy as np
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def load_review(root,analysis,original_part):
 root=Path(root);review=root.parent/'s41_numerical_review';folder=review/analysis;decision=folder/'decision.json'
 if not decision.exists():return None
 j=json.loads(decision.read_text())
 if j.get('status')!='adopt_A':return None
 snap=json.loads((review/'producer_provenance_snapshot.json').read_text())
 for name,h in snap['hashes'].items():assert sha(root/name)==h,'Retry source changed: '+name
 assert sha(folder/'original'/original_part.name)==sha(original_part),'Retry does not address current result'
 A=folder/'A';B=folder/'B';diags=[];fixed=[]
 for d in [A,B]:
  x=json.loads((d/'diagnostics.json').read_text());diags.append(x)
  assert x['clean'] and x['exact_data'] and x['ok'] and x['mode_status']==0 and x['misc_mode_status']==0
  assert x['fixed_finite'] and x['hyper_finite'] and x['pointwise_finite'] and x['cov_eigen_min']>0
  assert not x.get('internal_warnings') and not x.get('R_warnings') and not x.get('error')
  fixed.append(pd.read_csv(d/'fixed.csv').set_index('cov'))
 x,y=diags;f,g=fixed
 for table in fixed:assert np.isfinite(table.to_numpy(dtype=float)).all() and (table.sd>0).all()
 assert x['input_sha256']==y['input_sha256'] and list(f.index)==list(g.index)
 scale=(f.sd+g.sd)/2
 assert abs(x['WAIC']-y['WAIC'])<=1 and abs(x['log_posterior_mode']-y['log_posterior_mode'])<=.1
 assert (abs(f['mean']-g['mean'])/scale).max()<=.05
 for col in ['0.025quant','0.975quant']:assert (abs(f[col]-g[col])/scale).max()<=.1
 classify=lambda d:np.where(d['0.025quant']>0,1,np.where(d['0.975quant']<0,-1,0))
 assert np.array_equal(classify(f),classify(g))
 part=A/original_part.name;status=A/f'S41_status_{analysis}.json'
 assert Path(j['part']).resolve()==part.resolve()
 st=json.loads(status.read_text());st.update(waic=x['WAIC'],waic_finite=True,verified_isolated_retry=True,strategy='strict_A')
 files=[A/'fit.rds',B/'fit.rds',A/'input.rds',B/'input.rds',decision,A/'diagnostics.json',B/'diagnostics.json',A/'fixed.csv',B/'fixed.csv',part,status,review/'producer_provenance_snapshot.json']
 return {'part':part,'status':status,'diagnostics':st,'provenance':{'decision':str(decision),'selected':'A','independent_AB_checks_passed':True,'files':{str(p):sha(p) for p in files}}}
