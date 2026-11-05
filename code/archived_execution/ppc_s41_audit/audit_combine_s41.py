"""Read-only parent result audit; combine only into own isolated snapshot."""
from pathlib import Path
import json,hashlib,shutil,subprocess,sys,re
import pandas as pd
import numpy as np
from verify_execution_provenance import validate as validate_execution_provenance
from validate_conditional import validate_status as validate_conditional_status
from load_s41_retry import load_review
O=Path(__file__).resolve().parent;F=O.parent/'refit_final_20261104';R=F/'04_분석/BRCA/results_std'
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
issues=[];status_file=F/'sensitivity_group_execution.json'
if not status_file.exists():
 (O/'S41_conclusion_status.json').write_text(json.dumps({'state':'waiting_for_parent_groups','combined':False},indent=2));sys.exit(0)
execution_provenance=validate_execution_provenance(F)
(O/'S41_execution_provenance_validation.json').write_text(json.dumps(execution_provenance,ensure_ascii=False,indent=2))
if not execution_provenance['valid']:
 report={'state':'blocked_execution_provenance','combined':False,'issues':execution_provenance['issues'],'execution_provenance':execution_provenance}
 (O/'S41_conclusion_status.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
 print('S41 pre-fit provenance validation failed');sys.exit(2)
exits=json.loads(status_file.read_text())
if len(exits)!=6 or any(x['exit_code']!=0 for x in exits):issues.append('Not all six parent groups exited successfully')
s7=pd.read_csv(R/'S7_multiverse_summary.csv');rob=s7.loc[s7.robust,'cov'].tolist()
if len(rob)!=4:issues.append('Expected final robust count4')
if sha(F/'input_graph.graph')!=sha(R/'nb_BRCA.graph'):issues.append('Actual S41 graph differs from canonical graph')
panel=pd.read_csv(F/'04_분석/BRCA/panel_BRCA_std.csv');sido=pd.read_csv(R/'S41_panel_sido.csv',dtype={'sido':str})
expected=['C1_poisson_lognormal','C2_rw2_year','C3_ar1_year','C4_sido_year_iid','B3_smooth_forced','D1_spatial_plus_cc2_K23','E1_lag1','E2_lag2','F3_extreme_R7','G1_kchs_ma3']
expected += ['F1_drop_sido_'+s for s in sorted(sido.sido.unique())]
expected += ['F2_drop_year_'+str(y) for y in sorted(panel.year.unique())]+['H4_odd_years','H4_even_years']
expected += ['G2_noise_r%02d'%i for i in range(1,21)]+['E0_nolag_from2012','E0_nolag_from2013','drop_dyslipidemia','drop_imputed_hh_emp']
if rob:expected+=['bw_twoway']
actual={p.stem.removeprefix('S41_part_') for p in R.glob('S41_part_*.csv')}
missing=sorted(set(expected)-actual);unexpected=sorted(actual-set(expected))
if missing:issues.append('Missing analyses: '+','.join(missing))
if unexpected:issues.append('Unexpected analyses: '+','.join(unexpected))
rows=[];sources=[];diagnostics_by_analysis={};adoptions={};effective_files={}
config=json.loads((F/'04_분석/BRCA/config_BRCA_std.json').read_text())
mainvars=list(dict.fromkeys(config['sets_after_vif'][config['rules']['main_set']]+config['forced']))
def expected_rows(an):
 if an.startswith('F1_drop_sido_'):
  x=panel[['year','fbd_region']].merge(sido,on=['year','fbd_region'],validate='one_to_one');return int((x.sido!=an.removeprefix('F1_drop_sido_')).sum())
 if an.startswith('F2_drop_year_'):return int((panel.year!=int(an.rsplit('_',1)[1])).sum())
 if an in ['H4_odd_years','H4_even_years']:return int((panel.year%2==(1 if an=='H4_odd_years' else 0)).sum())
 if an.startswith('E0_nolag_from'):return int((panel.year>=int(an[-4:])).sum())
 names={'D1_spatial_plus_cc2_K23':'splus_cc2_K23','E1_lag1':'lag1','E2_lag2':'lag2','F3_extreme_R7':'extreme','G1_kchs_ma3':'ma3'}
 if an in names:
  x=pd.read_csv(R/f'S41_panel_{names[an]}.csv');cols=[c for c in x.columns[2:] if c in mainvars]
  d=panel[['year','fbd_region']].merge(x[['year','fbd_region']+cols],on=['year','fbd_region'],how='left',validate='one_to_one');return int(d[cols].notna().all(axis=1).sum())
 return len(panel)
for an in expected:
 f=R/f'S41_part_{an}.csv';st=R/f'S41_status_{an}.json'
 if not f.exists():continue
 if not st.exists():issues.append('No status: '+an);continue
 j=json.loads(st.read_text())
 if j.get('internal_warnings') or j.get('success') is not True or j.get('strategy')!='default' or j.get('mode_status') not in (0,None) or j.get('misc_mode_status') not in (0,None):
  try:reviewed=load_review(F,an,f)
  except Exception as error:issues.append('Retry adoption verification failed: '+an+': '+str(error));reviewed=None
  if reviewed:
   f=reviewed['part'];st=reviewed['status'];j=reviewed['diagnostics'];adoptions[an]=reviewed['provenance']
 d=pd.read_csv(f);effective_files[an]=str(f)
 diagnostics_by_analysis[an]=j
 required=['ok','mode_status','misc_mode_status','internal_warnings','fixed_finite','hyper_finite','waic_finite','waic']
 for key in required:
  if key not in j:issues.append('Missing diagnostic '+key+': '+an)
 for key in ['fixed_finite','hyper_finite','waic_finite']:
  if j.get(key) is not True:issues.append('Failed numerical diagnostic '+key+': '+an)
 if j.get('ok') not in (None,True):issues.append('INLA ok not true: '+an)
 for key in ['mode_status','misc_mode_status']:
  value=j.get(key)
  if value is not None and any(x!=0 for x in (value if isinstance(value,list) else [value])):issues.append('Nonzero '+key+': '+an)
 warnings=j.get('internal_warnings')
 if warnings:issues.append('Internal warnings require scientific review: '+an)
 if j.get('strategy')!='default' and not j.get('verified_isolated_retry'):issues.append('Nondefault numerical strategy requires review: '+an)
 if not isinstance(j.get('waic'),(int,float)) or not np.isfinite(j['waic']):issues.append('Invalid diagnostic WAIC: '+an)
 wanted=set(mainvars)
 if an=='B3_smooth_forced':wanted-=set(['sex_ratio','screening_uptake'])
 if an=='drop_dyslipidemia':wanted.discard('screening_normalB_dyslipidemia_share')
 if an=='drop_imputed_hh_emp':wanted-=set(['regular_employee_share','single_person_household_share','female_employment_rate'])
 if an=='bw_twoway':wanted=(wanted-set(rob))|{v+suffix for v in rob for suffix in ['__between','__within']}
 if set(d['cov'])!=wanted:issues.append('Coefficient inventory mismatch: '+an)
 expected_n=expected_rows(an)
 if j.get('n_obs')!=expected_n or set(d.n_obs)!={expected_n}:issues.append('Observation inventory mismatch: '+an)
 if j.get('success') is not True:issues.append('Fit failed: '+an)
 if set(d.analysis)!={an} or d['cov'].duplicated().any():issues.append('Invalid coefficient keys: '+an)
 if set(d.n_robust)!={len(rob)}:issues.append('Robust-count metadata mismatch: '+an)
 if not {'beta','IRR','lo','hi','cred'}.issubset(d):issues.append('Missing precise coefficient fields: '+an)
 else:
  if not np.isfinite(d[['beta','IRR','lo','hi']]).all().all():issues.append('Nonfinite coefficients: '+an)
  if not np.allclose(d.IRR,np.exp(d.beta),rtol=1e-12,atol=1e-12):issues.append('Rounded/inconsistent IRR: '+an)
 rows.append({'analysis':an,'success':j.get('success'),'strategy':j.get('strategy'),'n_obs':j.get('n_obs'),'expected_n_obs':expected_n,'coefficient_rows':len(d),'sha256':sha(f)})
 sources += [f,st]
 hot=f.parent/f'S41_hot_{an}.csv'
 if hot.exists():sources.append(hot)
if not (R/'S41_oos2024.csv').exists():issues.append('Missing conditional2024 check')
code=(F/'01_스크립트/88g_sens_suite.R').read_text()
conditional_status=R/'S41_status_conditional_2024.json'
conditional_output=R/'S41_oos2024.csv'
conditional_adoption=None
cp=O/'conditional_review/adoption_map.json'
if cp.exists():
 ca=json.loads(cp.read_text())
 if ca.get('status')=='adopt_A':
  try:
   assert sha(conditional_status)==ca['original_status_sha256'] and sha(conditional_output)==ca['original_output_sha256']
   for name,h in ca['hashes'].items():assert sha(Path(name))==h
   assert ca['metrics']['same_input'] and ca['metrics']['same_formula'] and ca['metrics']['same_credible_classes']
   assert abs(ca['metrics']['delta_WAIC'])<=1 and ca['metrics']['max_beta_delta_sd']<=.05 and ca['metrics']['max_endpoint_delta_sd']<=.1
   for d in [ca['A'],ca['B']]:assert not d.get('internal_warnings') and d['ok'] and d['mode_status']==0 and d['misc_mode_status']==0
   conditional_status=Path(ca['status_file']);conditional_output=Path(ca['output']);conditional_adoption={'map':str(cp),'sha256':sha(cp),'precision':ca['precision']}
  except Exception as e:issues.append('Conditional retry adoption invalid: '+str(e))
conditional_audit={'status_file':str(conditional_status),'expected_n_observed':2977,'expected_n_heldout':229}
if not conditional_status.exists():issues.append('Missing separate conditional2024 status')
else:
 j=json.loads(conditional_status.read_text());j['verified_isolated_retry']=conditional_adoption is not None;conditional_audit['diagnostics']=j;conditional_audit['adoption']=conditional_adoption
 issues.extend(validate_conditional_status(j));sources.append(conditional_status)
 conditional_audit['status_sha256']=sha(conditional_status)
 if int((panel.year==2024).sum())!=229 or int((panel.year!=2024).sum())!=2977:issues.append('Conditional2024 panel year counts mismatch')
 if conditional_output.exists():
  oos=pd.read_csv(conditional_output)
  cols=['n','coverage95','median_abs_pct_error','spearman_obs_pred']
  if len(oos)!=1 or not set(cols).issubset(oos):issues.append('Invalid conditional2024 output schema')
  elif not np.isfinite(oos[cols]).all().all() or oos.n.iloc[0]!=229 or not 0<=oos.coverage95.iloc[0]<=1 or not -1<=oos.spearman_obs_pred.iloc[0]<=1 or oos.median_abs_pct_error.iloc[0]<0:issues.append('Invalid conditional2024 metrics')
conditional_audit['r_seed_202611042']=bool(re.search(r'set.seed\(202611042L?\)',code))
conditional_audit['explicit_inla_posterior_seed']=bool(re.search(r'inla.posterior.sample\(200,m,[^;\n]*seed\s*=\s*202611042L?',code))
conditional_audit['interpretation']='Conditional held-out counts with 2024 E/covariates retained, not prospective forecasting'
conditional_audit['effective_precision']='unrounded_selected_retry' if conditional_adoption else 'original_display_rounded'
conditional_audit['selected_retry_seed']=202611042 if conditional_adoption else None
if not conditional_audit['r_seed_202611042'] and not conditional_audit['explicit_inla_posterior_seed']:issues.append('Conditional2024 seed not recorded in hashed source')
status_expected=set(expected)|{'conditional_2024'}
status_actual={p.stem.removeprefix('S41_status_') for p in R.glob('S41_status_*.json')}
if status_actual!=status_expected:issues.append('Status inventory mismatch: missing='+str(sorted(status_expected-status_actual))+' unexpected='+str(sorted(status_actual-status_expected)))

waic_rounded='WAIC=round(m$waic$waic' in code
# jsonlite toJSON defaults to four decimal places: raw expression alone is insufficient.
status_precision=bool(re.search(r'toJSON\(diag,[^\n]*digits\s*=\s*NA',code))
waic_precision_label='full_precision_status' if status_precision else 'status_JSON_4_decimal_places'
# Authorized precision: retain status four decimals; never describe this as unrounded.
if not rob:issues.append('Empty robust set requires explicit zero-row schema adaptation before original89g')
report={'state':'blocked_validation' if issues else 'validated_ready_to_combine','combined':False,'expected_count':len(expected),'expected_status_count':len(status_expected),'conditional_2024':conditional_audit,'actual_count':len(actual),'expected':sorted(expected),'missing':missing,'unexpected':unexpected,'issues':issues,'robust_variables':rob,'NCO':'intentionally excluded by parent NCO empty; not validated negative-control outcomes','full_precision_WAIC':status_precision,'WAIC_precision':waic_precision_label,'diagnostics':diagnostics_by_analysis,'execution_provenance':execution_provenance,'isolated_retry_adoption_map':adoptions,'effective_part_files':effective_files,'per_analysis':rows,'parent_writer_sha256':sha(F/'01_스크립트/88g_sens_suite.R')}
(O/'S41_inventory_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
pd.DataFrame(rows).to_csv(O/'S41_analysis_status.csv',index=False)
if issues:
 (O/'S41_conclusion_status.json').write_text(json.dumps(report,ensure_ascii=False,indent=2));print('S41 blocked:',len(issues),'issues');sys.exit(2)
# Never write into parent results, and never copy old-data results.
C=O/'S41_combine_snapshot';C.mkdir(exist_ok=True)
RAW=C/'raw_sources';RAW.mkdir(exist_ok=True)
sources += [R/n for n in ['S7_multiverse_summary.csv','refit_main_fixed_unrounded.csv']]+[conditional_output]
for f in sources:
 raw=RAW/f.name;dst=C/f.name
 if raw.exists():assert sha(raw)==sha(f),'Raw snapshot source changed: '+f.name
 else:shutil.copyfile(f,raw)
 shutil.copyfile(raw,dst)
# Replace display-one-decimal WAIC with status value; retain explicit precision label.
if waic_rounded:
 for an in expected:
  target=C/f'S41_part_{an}.csv';part=pd.read_csv(target)
  part['WAIC']=diagnostics_by_analysis[an]['waic'];part['WAIC_precision']='retry_diagnostics_16_digits' if an in adoptions else waic_precision_label;part.to_csv(target,index=False)
main=pd.read_csv(C/'refit_main_fixed_unrounded.csv')
main=main.assign(IRR=np.exp(main.beta),lo=np.exp(main.lo),hi=np.exp(main.hi))
main['cred']=(main.lo>1)|(main.hi<1)
main.to_csv(C/'S3_main_IRR.csv',index=False)
reference=main.set_index('cov').beta
for an in expected:
 target=C/f'S41_part_{an}.csv';part=pd.read_csv(target)
 kept=part[part['cov'].isin(rob)].copy()
 count=int((kept.cred & (np.sign(kept.beta)==np.sign(kept['cov'].map(reference)))).sum())
 part['robust_kept']=count;part['n_robust']=len(rob);part.to_csv(target,index=False)
# Revalidate immediately before combination to catch changes during the audit.
second=validate_execution_provenance(F)
if not second['valid'] or second['provenance_sha256']!=execution_provenance['provenance_sha256']:
 report.update(state='blocked_provenance_changed_during_audit',combined=False,second_validation=second)
 (O/'S41_conclusion_status.json').write_text(json.dumps(report,ensure_ascii=False,indent=2));sys.exit(2)
shutil.copyfile(F/'sensitivity_execution_provenance.json',C/'sensitivity_execution_provenance.json')
source=F/'01_스크립트/89g_sens_combine.py';s=source.read_text()
start=s.index('ROOT =');end=s.index('\nP =',start)
s=s[:start]+'R = os.environ["COMBINE_INPUT"]'+s[end:]
s=s.replace('round(math.exp(qb), 3)','math.exp(qb)').replace('round(math.exp(qb - 1.96 * math.sqrt(T)), 3)','math.exp(qb - 1.96 * math.sqrt(T))').replace('round(math.exp(qb + 1.96 * math.sqrt(T)), 3)','math.exp(qb + 1.96 * math.sqrt(T))')
s=s.replace('def r2(x): return float(Decimal(str(x)).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))','def r2(x): return float(x)')
# Fixed sign should use unrounded stored beta, not rounded IRR.
s=s.replace('sgn = np.sign(np.log(M3.IRR))','sgn = np.sign(M3.beta)')
s=s.replace('np.sign(np.log(a.IRR))','np.sign(a.beta)').replace('np.sign(np.log(b.IRR))','np.sign(b.beta)')
s='\n'.join(line for line in s.splitlines() if 'nc.to_csv(' not in line)+'\n'
s=s.replace('## 비교 결과\\n{nc.to_markdown(index=False)}\\n\\n','')
s=s.replace('q = np.log(b.IRR)', 'q = b.beta')
script=O/'89g_validated_combine.py';script.write_text(s)
import os
env=os.environ.copy();env['COMBINE_INPUT']=str(C)
with (O/'combine89g.log').open('w') as log:rc=subprocess.run([sys.executable,str(script),'BRCA','std'],env=env,stdout=log,stderr=subprocess.STDOUT).returncode
report.update(state='combined' if rc==0 else 'combine_failed',combined=rc==0,combine_exit_code=rc,output_directory=str(C),source_hashes={str(f):sha(f) for f in sources},limitations=['Rubin-style noise perturbation pooling uses CrI-based log-normal SE approximation; not proper MI','E-values ecological perSD associations, not proof of causal robustness','conditional2024 includes outcome-calibrated E and full-panel preprocessing; adopted conditional output unrounded when documented, otherwise original metrics display-rounded','No NCO fits were run; NCO export is intentionally omitted, not evidence'])
report['output_hashes']={p.name:sha(p) for p in C.iterdir() if p.is_file()}
report['precision_notes']={'coefficients':'source numeric precision retained; signs use beta','WAIC_original':'status JSON four decimal places','WAIC_adopted_retry':'diagnostics sixteen digits','hotspot_exports':'producer RR/probability rounded four decimals; threshold counts were computed before rounding','hot_jaccard':'producer three decimals','conditional2024':'selected unrounded retry summary'}
report['combiner_sha256']=sha(script)
report['combiner_original_sha256']=sha(source)
(O/'S41_conclusion_status.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
sys.exit(rc)
