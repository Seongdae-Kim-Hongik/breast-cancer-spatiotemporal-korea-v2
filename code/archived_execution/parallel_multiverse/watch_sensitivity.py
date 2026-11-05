from pathlib import Path
import subprocess,time,json
w=Path(__file__).resolve().parent.parent
r=w/'refit_final_20261104'
deadline=time.monotonic()+8*3600
while time.monotonic()<deadline:
 p=r/'typeI_and_completion.log'
 if p.exists() and 'ALL DONE BRCA_std' in p.read_text():
  result=subprocess.run(['python3',str(w/'parallel_multiverse/run_sensitivity_groups.py')])
  raise SystemExit(result.returncode)
 exits=r/'worker_exit_codes.json'
 if exits.exists() and any(x['exit_code']!=0 for x in json.loads(exits.read_text())): raise RuntimeError('multiverse worker failed')
 time.sleep(15)
raise TimeoutError('Core completion deadline exceeded; no sensitivities started')
