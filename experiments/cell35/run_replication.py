from pathlib import Path
from itertools import combinations
from concurrent.futures import ThreadPoolExecutor
import numpy as np, subprocess, json, time
from analysis35 import DT, short_paths, checkpoint, support, jsonable
import argparse
p=argparse.ArgumentParser();p.add_argument('--run',default='DCU_Cell35_run');args=p.parse_args()
root=Path(args.run).resolve();out=root/'replication';out.mkdir(exist_ok=True)
protocol='Fixed earliest length-3 delta-2 pair in each early source, carried unchanged to late epoch. 32 fresh continuations per source epoch; complete-burst local32 and global8192 readings. Added after first pilot world inspection; every original result retained.'
(out/'SCOPE_BEFORE_RUN.txt').write_text(protocol+'\nSeeds 20360920+1000*world_index+replicate; paired numbers across early/late, state-dependent draws diverge. No phase fit, pair reselection, changed sampler or changed costs.\n')
tasks=[]
for wi,(H,seed) in enumerate((H,s) for H in (0,6) for s in range(20350920,20350924)):
 prefix=root/f'H{H}_s{seed}';nodes=np.loadtxt(prefix.with_suffix('.nodes.tsv'),dtype=np.int64)
 ev=np.memmap(prefix.with_suffix('.events.bin'),dtype=DT,mode='r');taus=np.array(ev['tau'],dtype=np.int64)
 idx=int(np.searchsorted(taus,128));cp=checkpoint(ev,nodes,idx,H);pp=short_paths(nodes,cp['n'],cp['step'])[3]
 pair=None
 for i,j in combinations(range(len(pp)),2):
  if len(set(support(pp[i])) ^ set(support(pp[j])))==2:
   pair=(i,j);break
 if pair is None:raise RuntimeError('No fixed early pair; do not substitute a favorable later object.')
 pa,pb=(pp[k] for k in pair);a,b=support(pa),support(pb)
 for epoch in (128,4194304):
  idx=int(np.searchsorted(taus,epoch));cp=checkpoint(ev,nodes,idx,H)
  fn=out/f'H{H}_s{seed}_age{epoch}.input.tsv';dest=fn.with_suffix('.output.tsv')
  with fn.open('w') as f:
   f.write(' '.join(map(str,(H,cp['F'],cp['P'],cp['step'],cp['tau'],cp['n'])))+'\n')
   for row in nodes[:cp['n']]:
    vals=list(map(int,row[1:8]))
    if vals[-1]>cp['step']:vals[-2:]=[0,0]
    f.write(' '.join(map(str,vals))+'\n')
   f.write(' '.join(map(str,(len(a),*a,len(b),*b)))+'\n')
  tasks.append(dict(H=H,source_seed=seed,epoch=epoch,checkpoint=cp,paths=(pa,pb),input_file=str(fn),output_file=str(dest),fork_seed_start=20360920+1000*wi))

def run(task):
 cmd=[str(root/'native35'),'--fork',task['input_file'],task['output_file'],str(task['fork_seed_start']),'32','200000000']
 subprocess.run(cmd,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 rows=[]
 for line in Path(task['output_file']).read_text().splitlines():
  parts=line.split();rr,ss,mode,complete,steps,tau,na,nb,un,n,F=parts
  row=dict(rep=int(rr),seed=int(ss),mode=mode,complete=bool(int(complete)),elapsed_steps=int(steps),elapsed_tau=int(tau),N_A=int(na),N_B=int(nb),union_services=int(un),n=int(n),F=int(F))
  d=row['N_A']-row['N_B'];row['difference']=d;row['correct_word']={str(Q):float((1+2*np.cos(2*np.pi*(d%Q)/Q))**2/9) for Q in (27,81)};rows.append(row)
 task['readouts']=rows;print('FORK DONE',task['H'],task['source_seed'],task['epoch'],flush=True);return task
start=time.perf_counter()
with ThreadPoolExecutor(max_workers=4) as ex:results=list(ex.map(run,tasks))
(out/'replication_results.json').write_text(json.dumps(jsonable(dict(protocol=protocol,tasks=results,wall_seconds=time.perf_counter()-start)),indent=2))
