from pathlib import Path
import os,subprocess,json,time,concurrent.futures,hashlib
W=Path(__file__).resolve().parent.parent;root=W/'refit_final_20261104';res=root/'04_분석/BRCA/results_std'
assert (root/'typeI_and_completion.log').exists() and 'ALL DONE BRCA_std' in (root/'typeI_and_completion.log').read_text(),'Core/multiverse/typeI must finish'
a=json.loads((res/'S7_execution_audit.json').read_text());assert a['success']==a['planned']==200,'Investigate missing subsets before sensitivities'
assert not list(res.glob('S41_part_*.csv')),'Never silently reuse/overwrite sensitivity outputs'
provenance_files=[root/'01_스크립트/88g_sens_suite.R',root/'input_graph.graph',root/'04_분석/BRCA/panel_BRCA_std.csv',res/'S7_multiverse_results.csv',res/'S7_multiverse_summary.csv',res/'S7_execution_status.csv',res/'refit_main_fixed_unrounded.csv',res/'S41_panel_splus_cc2_K23.csv',res/'S41_splus_cc2_K23_definition.json']
provenance={'created':time.strftime('%Y-%m-%d %H:%M:%S%z'),'hashes':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in provenance_files},'max_parallel':3,'threads_per_fit':1,'calendar_gaps':'year-min(year)+1; explicit n and values preserve missing calendar positions','status_diagnostics':'mode.status/misc.mode.status/ok/internal warnings/fixed and hyper finite/WAIC'}
(root/'sensitivity_execution_provenance.json').write_text(json.dumps(provenance,indent=2))

e=os.environ.copy();e.update(TAG='BRCA',ST_VER='std',ST_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',TMPDIR='/tmp',FOCAL='mean_age_first_marriage_wife',SMOOTH='sex_ratio,screening_uptake',NCO='',DROPS=json.dumps({'dyslipidemia':['screening_normalB_dyslipidemia_share'],'imputed_hh_emp':['regular_employee_share','single_person_household_share','female_employment_rate']}))
def run(g):
 env=e.copy();env['GROUP']=str(g);t=time.time()
 with (root/f'sensitivity_group{g}.log').open('w') as log:p=subprocess.run(['/usr/local/bin/Rscript','01_스크립트/88g_sens_suite.R'],cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT)
 result={'group':g,'exit_code':p.returncode,'seconds':time.time()-t};print(result,flush=True);return result
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:status=list(pool.map(run,[1,2,3,4,5,6]))
(root/'sensitivity_group_execution.json').write_text(json.dumps(status,indent=2));assert all(x['exit_code']==0 for x in status),status
# Do not combine blindly: independently inspect per-analysis status inventory first.
print('S41 groups terminated. Per-analysis completion audit and combination still required.',flush=True)
