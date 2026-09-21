"""Read-only replay of the compact Cell35 results; no new native trajectory."""
from pathlib import Path
from collections import Counter
import gzip,json,math
import numpy as np
ROOT=Path(__file__).resolve().parent
counts=Counter()
files=sorted((ROOT/'results').glob('H*.analysis.json.gz'))
assert len(files)==8
for file in files:
 with gzip.open(file,'rt') as f:world=json.load(f)
 nodes=np.loadtxt(file.with_name(file.name.replace('.analysis.json.gz','.nodes.tsv')),dtype=np.int64)
 for target,ep in world['epochs'].items():
  if not ep['reached']:continue
  cp=ep['checkpoint'];assert int(target)<=cp['tau']<=int(target)+2
  for co in ep['cohorts'].values():
   paths=co['paths']
   for root,tokens in paths:
    here=root
    for z,e in tokens:
     assert z<cp['n'] and int(nodes[z,e+1])==here and 0<int(nodes[z,7])<=cp['step'];here=z
    counts['complete_native_paths']+=1
   for pair in co['pair_readouts']:
    a,b=({z for z,e in paths[i][1]} for i in (pair['i'],pair['j']))
    assert len(a)==len(b)==pair['length'] and len(a^b)==pair['delta']
    for mode in ('global_windows','local_windows'):
     for threshold,rd in pair[mode].items():
      na,nb=rd['N_A'],rd['N_B'];d=na-nb
      assert d==rd['difference'] and na+nb==rd['private_activity']+2*rd['shared_activity']
      assert rd['union_activity']==rd['private_activity']+rd['shared_activity']
      if rd['complete']:
       exposure=rd['elapsed_tau'] if mode=='global_windows' else rd['union_activity']
       assert int(threshold)<=exposure<=int(threshold)+2
      for Q in (27,81):
       omega=2*math.pi*(d%Q)/Q
       p=(3+4*math.cos(omega)+2*math.cos(2*omega))/9
       assert abs(p-rd['correct_word'][str(Q)])<3e-14
      if a==b:assert d==0
      counts['phase_readout_records']+=1
   for cc in co['collective_availability'].values():
    for ex in cc.get('examples',{}).values():
     ss=[tuple(z for z,e in p[1]) for p in ex['paths']];assert len(set(ss))==4
     balance=Counter();entry=Counter()
     for sign,path in zip((1,1,-1,-1),ex['paths']):
      for z,e in path[1]:balance[z]+=sign;entry[z,e]+=sign
     assert all(v==0 for v in balance.values())
     assert ex['entry_token_balance']==all(v==0 for v in entry.values())
     assert ex['sum_difference']==0
     counts['four_distinct_path_certificates']+=1
rep=json.loads((ROOT/'replication/replication_results.json').read_text())
assert len(rep['tasks'])==16
for task in rep['tasks']:
 assert len(task['readouts'])==64
 for rd in task['readouts']:
  d=rd['N_A']-rd['N_B'];assert d==rd['difference']
  if rd['complete']:
   value=rd['elapsed_tau'] if rd['mode']=='global8192' else rd['union_services']
   goal=8192 if rd['mode']=='global8192' else 32
   assert goal<=value<=goal+2
  for Q in (27,81):
   p=(1+2*math.cos(2*math.pi*(d%Q)/Q))**2/9
   assert abs(p-rd['correct_word'][str(Q)])<3e-14
  counts['fresh_confirmation_readouts']+=1
print('PASS: compact exact-input/readout replay')
print(json.dumps(dict(counts),indent=2))
print('Native generation is NOT rerun; no matter or physical epoch is inferred.')
