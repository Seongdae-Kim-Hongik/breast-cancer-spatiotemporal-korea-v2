from pathlib import Path
import json,os,subprocess,time,hashlib,sys
W=Path(__file__).resolve().parent.parent
root=Path(sys.argv[1]).resolve();ready=root/'DATA_READY.json';assert ready.exists(),'Inputs/graph must be explicitly validated first'
d=json.loads(ready.read_text());b=root/'04_분석/BRCA';r=b/'results_std'
for rel,h in d['sha256'].items():assert hashlib.sha256((root/rel).read_bytes()).hexdigest()==h,rel
manifest=[]
def run(cmd,label,env):
 t=time.time()
 with (root/f'{label}.log').open('w') as log:x=subprocess.run(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT)
 manifest.append({'stage':label,'exit_code':x.returncode,'seconds':time.time()-t});(root/'execution_manifest.json').write_text(json.dumps(manifest,indent=2));print(manifest[-1],flush=True)
 assert x.returncode==0,label
base=os.environ.copy();base.update(ST_VER='std',ST_THREADS='1',OMP_NUM_THREADS='1',TMPDIR='/tmp',BRCA_GRAPH=str(root/'input_graph.graph'))
e=base.copy();e['ST_ONLY_CORE']='1';run(['/usr/local/bin/Rscript','01_스크립트/81_fit_inla.R','BRCA'],'core',e)
procs=[]
for i in range(1,5):
 out=root/'multiverse_workers'/f'worker_{i}';out.mkdir(parents=True,exist_ok=False)
 e=base.copy();e.update(BRCA_INPUT_ROOT=str(root),AUDIT_OUT=str(out),MV_WORKER=str(i),MV_WORKERS='4',MV_SPECS=str(r/'S7_multiverse_specs.csv'))
 log=(out/'worker.log').open('w');p=subprocess.Popen(['/usr/local/bin/Rscript',str(W/'parallel_multiverse/worker.R'),'BRCA'],cwd=root,env=e,stdout=log,stderr=subprocess.STDOUT);procs.append((i,p,log))
status=[]
for i,p,log in procs:status.append({'worker':i,'exit_code':p.wait()});log.close()
(root/'worker_exit_codes.json').write_text(json.dumps(status,indent=2));assert all(x['exit_code']==0 for x in status),status
run(['python3',str(W/'parallel_multiverse/aggregate.py'),str(root)],'aggregate',base)
run(['/usr/local/bin/Rscript','01_스크립트/81_fit_inla.R','BRCA'],'typeI_and_completion',base)
print('CORE + MULTIVERSE + TYPE I COMPLETE; additional sensitivities and final package still required.',flush=True)
