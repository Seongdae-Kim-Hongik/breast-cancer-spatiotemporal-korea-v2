"""Validate immutable pre-fit evidence against current final-root inputs."""
from pathlib import Path
import hashlib,json,re
REQUIRED={
 '01_스크립트/88g_sens_suite.R','input_graph.graph','04_분석/BRCA/panel_BRCA_std.csv',
 *('04_분석/BRCA/results_std/'+n for n in ['S7_multiverse_results.csv','S7_multiverse_summary.csv','S7_execution_status.csv','refit_main_fixed_unrounded.csv','S41_panel_splus_cc2_K23.csv','S41_splus_cc2_K23_definition.json'])}
def validate(root):
 root=Path(root).resolve();p=root/'sensitivity_execution_provenance.json'
 result={'valid':False,'issues':[],'files':[]}
 if not p.is_file():result['issues'].append('Missing pre-fit sensitivity_execution_provenance.json');return result
 raw=p.read_bytes();result['provenance_sha256']=hashlib.sha256(raw).hexdigest()
 try:record=json.loads(raw)
 except (ValueError,UnicodeError) as e:result['issues'].append('Invalid provenance JSON: '+str(e));return result
 hashes=record.get('hashes',{});result['created']=record.get('created')
 if not isinstance(hashes,dict):result['issues'].append('hashes must be mapping');return result
 for name in sorted(REQUIRED-set(hashes)):result['issues'].append('Missing required hash: '+name)
 if record.get('max_parallel')!=3 or record.get('threads_per_fit')!=1:result['issues'].append('Unexpected execution parallelism metadata')
 for name,expected in sorted(hashes.items()):
  target=(root/name).resolve()
  if Path(name).is_absolute() or root not in target.parents:result['issues'].append('Non-finalroot provenance path: '+name);continue
  actual=hashlib.sha256(target.read_bytes()).hexdigest() if target.is_file() else None
  match=isinstance(expected,str) and bool(re.fullmatch('[0-9a-f]{64}',expected)) and actual==expected
  result['files'].append({'path':name,'expected_sha256':expected,'actual_sha256':actual,'match':match})
  if not match:result['issues'].append('Missing or changed input: '+name)
 result['valid']=not result['issues'];return result
