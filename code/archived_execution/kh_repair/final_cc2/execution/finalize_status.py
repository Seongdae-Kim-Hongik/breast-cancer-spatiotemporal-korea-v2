from pathlib import Path
import json,hashlib,pandas as pd
RUN=Path(__file__).resolve().parent;dep=json.loads((RUN.parent/'dependency_manifest.json').read_text());R=Path(dep['panel_source']).parent/'results_std'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
p=RUN/'kh_diagnostics.json';d=json.loads(p.read_text())
for ty,v in d.items():
 assert v['status']=='COMPLETED_DIAGNOSTICS_PASS'
 warns=v['fit_warnings'];haswarn=bool(warns)
 v['diagnostics_pass']=True;v['provisional']=haswarn
 if haswarn:v['status']='COMPLETED_WITH_OPTIMIZATION_WARNING'
 out=R/f'S43_part_KH_type{ty}.csv';tab=pd.read_csv(out);tab['numerical_status']=v['status'];tab['provisional']=v['provisional'];tab['numerical_strategy']=v['execution']['selected_strategy'];tab['numerical_warning']=warns if isinstance(warns,str) else '; '.join(warns)
 tab.to_csv(out,index=False)
p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
iv=d['IV'];assert iv['provisional']
p=R/'S46_KH4_spatialRR.csv';tab=pd.read_csv(p,dtype={'sigungu_cd':str});tab['provisional']=True;tab.to_csv(p,index=False)
p=R/'S46_KH4_spatial_summary.csv';tab=pd.read_csv(p);tab['provisional']=True;tab['numerical_status']=iv['status'];tab['numerical_warning']=iv['fit_warnings'];tab.to_csv(p,index=False)
p=RUN/'S46_extraction_provenance.json';v=json.loads(p.read_text());v.update(provisional=True,numerical_warning=iv['fit_warnings'],spatial_file_sha256=sha(R/'S46_KH4_spatialRR.csv'));p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n')
p=RUN/'S46_summary_provenance.json';v=json.loads(p.read_text());v['summary'].update(provisional=True,numerical_status=iv['status'],numerical_warning=iv['fit_warnings']);v.update(spatial_rr_sha256=sha(R/'S46_KH4_spatialRR.csv'),summary_sha256=sha(R/'S46_KH4_spatial_summary.csv'));p.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
pre=json.loads((RUN/'preflight.json').read_text())
for path,h in pre['input_hashes'].items():assert sha(Path(path))==h,path
outputs=[R/f'S43_{prefix}_KH_type{ty}.{ext}' for ty in ['II','III','IV'] for prefix,ext in [('part','csv'),('fit','rds')]]+[R/'S46_KH4_spatialRR.csv',R/'S46_KH4_spatial_summary.csv']
manifest={'completed':True,'fits_run':3,'max_simultaneous_fits':1,'num_threads':'1:1','same_final_panel_and_graph_verified':True,'no_partial_period_fit':True,'no_external_EB_Gaussian_retry':True,'S46_reused_KH_IV':True,'KH_IV_provisional_optimization_warning':iv['fit_warnings'],'no_longlag_pop80_or_all_combo_fits':True,'output_sha256':{str(f):sha(f) for f in outputs}}
(RUN/'RUN_COMPLETE.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
print('Finalized8 designated result/RDS files; KH IV and S46 marked provisional. All preflight input hashes unchanged.')
