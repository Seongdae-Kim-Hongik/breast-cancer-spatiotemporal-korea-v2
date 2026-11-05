from pathlib import Path
import json,hashlib
import pandas as pd
import numpy as np
O=Path(__file__).resolve().parent;W=O.parent
F=W/'population_final_audit/final_stage/04_분석/BRCA'
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
p=pd.read_csv(O/'04_분석/BRCA/panel_BRCA_std.csv')
config=json.loads((O/'04_분석/BRCA/config_BRCA_std.json').read_text())
main=list(dict.fromkeys(config['sets_after_vif'][config['rules']['main_set']]+config['forced']))
assert len(main)==24 and 'pop80_share' not in main
jobs=[];definitions={};inputs={}
for sex in ['female','total']:
 source=F/f'S23_BRCA_pop80_{sex}.csv';x=pd.read_csv(source)
 assert len(x)==3206 and not x.duplicated(['year','fbd_region']).any()
 raw=x.pop80_share_raw;mu=float(raw.mean());sd=float(raw.std(ddof=1))
 np.testing.assert_allclose(x.pop80_share,(raw-mu)/sd,rtol=1e-12,atol=1e-12)
 q=p.merge(x,on=['year','fbd_region'],validate='one_to_one',how='left',sort=False)
 assert q.pop80_share.notna().all() and len(q)==3206
 pd.testing.assert_frame_equal(q[p.columns],p)
 design=np.column_stack([np.ones(len(q)),q[main+['pop80_share']]])
 assert np.linalg.matrix_rank(design)==26
 name=f'pop80_{sex}_contact';panel=O/'panels'/f'{name}.csv';q.to_csv(panel,index=False)
 jobs.append({'name':name,'panel':str(panel),'graph':str(O/'graphs/contact.graph'),'comparison_group':'full_panel_contact_3206'})
 definitions[name]={'role':'primary residual-age sensitivity' if sex=='female' else 'secondary reproduction sensitivity',
 'raw_definition':f'100 * {sex}-sex population aged80+ / {sex}-sex population all ages' if sex=='female' else '100 * both-sex population aged80+ / both-sex population all ages',
 'age_bands':['80-84','85-89','90-94','95-99','100+'],'scale':'pooled3206 observation sample SD; ddof1','mean_pct':mu,'sd_pct_ddof1':sd,'source':str(source),'source_sha256':sha(source),'main_variables':main+['pop80_share'],'n':len(q),'E':'unchanged final full-panel E'}
 inputs[str(panel.relative_to(O))]=sha(panel)
src=(O/'scripts/fit.R').read_text()
needle='MAINV<-unique(c(SETS[[RU$main_set]],FORCED))'
assert src.count(needle)==1
src=src.replace(needle,needle+'''\nstopifnot(length(MAINV)==24L,!"pop80_share"%in%MAINV)
MAINV<-c(MAINV,"pop80_share")
stopifnot(length(MAINV)==25L,"pop80_share"%in%names(DD))
writeLines(MAINV,file.path(RES,"actual_MAINV.txt"))''')
(O/'scripts/fit_pop80.R').write_text(src)
for f in ['scripts/fit_pop80.R','graphs/contact.graph','04_분석/BRCA/config_BRCA_std.json','04_분석/BRCA/panel_META_BRCA_std.csv','04_분석/BRCA/panel_INFO_BRCA_std.json','04_분석/BRCA/panel_BRCA_std.csv']:inputs[f]=sha(O/f)
manifest={'jobs':jobs,'definitions':definitions,'input_sha256':inputs,'never_both_pop80_variables_in_same_fit':True,'preserved_main_variable_count':24,'fitted_variable_count':25,'parent_core_mutated':False}
(O/'pop80_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
run=(O/'run_sequential.py').read_text().replace("O/'manifest.json'","O/'pop80_manifest.json'").replace("O/'scripts/fit.R'","O/'scripts/fit_pop80.R'").replace('ALL SIX COMPLETE','BOTH POP80 FITS COMPLETE')
(O/'run_pop80_sequential.py').write_text(run)
print('Prepared female-primary and total-secondary;25 actual variables, identical3206 rows/E, sample-SD verified')
