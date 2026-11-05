from pathlib import Path
import json,os,subprocess,hashlib,datetime
O=Path(__file__).resolve().parent
m=json.loads((O/'pop80_manifest.json').read_text())
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
env=os.environ.copy();env.update(AUDIT_ROOT=str(O),ST_VER='std',ST_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',TMPDIR='/tmp',AUDIT_MODE='main')
for job in m['jobs']:
 name=job['name'];out=O/'results'/name
 if (out/'COMPLETE').exists():
  print(name,'already complete; preserve',flush=True);continue
 assert not out.exists(),f'Partial output exists: review before restarting {name}'
 for f,h in m['input_sha256'].items():assert sha(O/f)==h,f
 out.mkdir();env.update(AUDIT_OUT=str(out),AUDIT_PANEL=job['panel'],AUDIT_GRAPH=job['graph'])
 print('START',name,datetime.datetime.now().isoformat(),flush=True)
 with (out/'fit.log').open('w') as log:
  result=subprocess.run(['/usr/local/bin/Rscript',str(O/'scripts/fit_pop80.R'),'BRCA'],env=env,stdout=log,stderr=subprocess.STDOUT)
 if result.returncode or not (out/'COMPLETE').exists():raise SystemExit(f'Fit failed: {name}; inspect private fit.log; sequential queue stopped')
 artifact={'job':job,'input_manifest_sha256':sha(O/'pop80_manifest.json'),'panel_sha256':sha(Path(job['panel'])),'graph_sha256':sha(Path(job['graph'])),'output_sha256':{f.name:sha(f) for f in out.iterdir() if f.is_file()}}
 (out/'hash_manifest.json').write_text(json.dumps(artifact,indent=2))
 print('DONE',name,datetime.datetime.now().isoformat(),flush=True)
print('BOTH POP80 FITS COMPLETE',flush=True)
