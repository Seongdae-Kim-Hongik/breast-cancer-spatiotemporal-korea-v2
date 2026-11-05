from pathlib import Path
import json,hashlib,numpy as np,pandas as pd
from spatial_plus_components import read_graph,spatial_basis,residualize,mm
D=Path(__file__).resolve().parent;m=json.loads((D/'dependency_manifest.json').read_text())
pfile=Path(m['panel_source']);B=pfile.parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert m['panel_authorized'] and sha(pfile)==m['authorized_panel_sha256']
for name,h in m['authorized_input_hashes'].items():assert sha(B/name)==h
assert sha(D/'chosen_nb_BRCA.graph')==m['graph_sha256']
config=json.loads((B/'config_BRCA_std.json').read_text());main=config['rules']['main_set'];assert main=='S1_narrative_full'
cols=list(dict.fromkeys(config['sets_after_vif'][main]+config['forced']))
p=pd.read_csv(pfile);order=json.loads((D/'region_order.json').read_text());idx=p.fbd_region.map(dict(zip(order,range(len(order)))))
assert idx.notna().all() and len(p)==3206 and sorted(p.year.unique())==list(range(2011,2025)) and not p[['year','fbd_region']].duplicated().any()
assert len(cols)==24 and np.isfinite(p[['cases','E']].to_numpy()).all()
idx=idx.to_numpy(dtype=int);adj=read_graph(D/'chosen_nb_BRCA.graph');basis,info=spatial_basis(adj,23)
assert info['components']==2 and info['component_sizes']==[227,2] and info['rank']==25
z,r,sd=residualize(p[cols].to_numpy(),idx,basis);obs=basis[idx]
err=float(np.max(np.abs(mm(obs.T,r))));zerr=float(np.max(np.abs(mm(obs.T,z))))
assert err<1e-9 and zerr<1e-8 and np.max(np.abs(z.std(0)-1))<1e-12
rng=np.random.default_rng(20261104);rotation,_=np.linalg.qr(rng.normal(size=(25,25)))
zr,_,_=residualize(p[cols].to_numpy(),idx,mm(basis,rotation));rot=float(np.max(np.abs(zr-z)));assert rot<1e-10
out=p[['year','fbd_region']].copy();out[cols]=z;out.to_csv(D/'S41_panel_splus_cc2_K23.csv',index=False)
# Retain key-specific covariate mapping and no reordering surprises.
back=pd.read_csv(D/'S41_panel_splus_cc2_K23.csv');assert back[['year','fbd_region']].equals(p[['year','fbd_region']]);assert np.allclose(back[cols],z,rtol=1e-14,atol=1e-14)
report=dict(passed=True,fits_run=0,analysis_id='D1_spatial_plus_cc2_K23',input_file='S41_panel_splus_cc2_K23.csv',panel_sha256=sha(pfile),graph_sha256=m['graph_sha256'],main_set=main,covariates=cols,rows=len(p),basis=info,residual_sd_ddof0=dict(zip(cols,sd.tolist())),max_raw_projection_inner_product=err,max_standardized_projection_inner_product=zerr,rotation_invariance_error=rot,output_sha256=sha(D/'S41_panel_splus_cc2_K23.csv'),definition='Project all two component indicators and 23 strictly positive low-frequency Laplacian modes; pooled over 2011-2024 full panel, all24 main covariates; residual population-SD z; outcome/E unchanged and not used in projection.',integration_gate='Check future HM_best best_set matches S1_narrative_full and actual MAINV matches this covariate list; do not silently merge incomplete columns.')
(D/'S41_splus_cc2_K23_definition.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else x,allow_nan=False)+'\n')
print('Spatial+ final cc2 PASS: 3206 rows,24 covariates,rank25, max orthogonality error',err)
