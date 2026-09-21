from pathlib import Path
from itertools import combinations,product
from fractions import Fraction as Q
from collections import Counter
from math import comb
import ast,gzip,json,time
import numpy as np
from analysis35 import DT,support,short_paths,aggregate
ROOT=Path(__file__).resolve().parent
import os
WORK=Path(os.environ.get('DCU35_RAW_RUN',str(ROOT/'DCU_Cell35_run')))

def original_definitions():
 env={};p=ROOT/'reference/native_constructor.py';tree=ast.parse(p.read_text())
 exec(compile(ast.Module([x for x in tree.body if isinstance(x,(ast.Import,ast.ImportFrom,ast.ClassDef,ast.FunctionDef))],type_ignores=[]),str(p),'exec'),env)
 p=ROOT/'reference/DCU_Mass_Cell_25.py';tree=ast.parse(p.read_text());f=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='_run_dcu_mass_25');step=next(x for x in f.body if isinstance(x,ast.FunctionDef) and x.name=='step')
 env['gamma']=3;exec(compile(ast.Module([step],type_ignores=[]),str(p),'exec'),env)
 return env

class Draw:
 def __init__(self,labels):self.labels=labels
 def sample(self,pool,q):
  assert len(self.labels)==q and len(set(self.labels))==q and all(0<=z<len(pool) for z in self.labels)
  return self.labels

def quiet(F,P,g,H):
 if H==0:
  assert F>=3*g;return F-3*g,P+6*g
 while g:
  if P>=24 and F>=9:
   k=min(g,F//9);F-=9*k;g-=k
   if not g:break
  assert F>=3
  fm=F-3;pm=P+6;quota=2*((max(1,pm//6)+1)//2)
  v=min(quota,H,fm,pm) if fm>=3 and pm>=6 else 0
  F,P=fm-v,pm-v;g-=1
 return F,P

def validate_history(prefix,counters):
 H=int(json.loads(prefix.with_suffix('.meta.json').read_text())['H'])
 env=original_definitions();st=env['DCUStructure']();F=P=tau=0
 for line in prefix.with_suffix('.audit.tsv').read_text().splitlines():
  t,n,f,p,*r=map(int,line.split());d=[z for z in r[:3] if z>=0];exp=r[3:]
  assert (len(st),F,P)==(n,f,p)
  F,P,event=env['step'](st,F,P,Draw(d),H);tau+=len(event['served'])
  assert [len(st),F,P,tau,event['cost'],event['relief_removed']]==exp
  counters['original_complete_steps']+=1
 nodes=np.loadtxt(prefix.with_suffix('.nodes.tsv'),dtype=np.int64,ndmin=2)
 assert tuple(st.parents)==tuple(None if i<2 else tuple(map(int,nodes[i,1:3])) for i in range(len(st)))
 assert set(st.recorded)=={i for i in range(len(st)) if 0<int(nodes[i,7])<=t}
 ev=np.memmap(prefix.with_suffix('.events.bin'),dtype=DT,mode='r');steps=np.array(ev['step'],dtype=np.int64)
 # Every missing record is an all-forced service. Audit the complete ledger, in blocks.
 for start in range(0,len(ev)-1,200000):
  end=min(start+200000,len(ev)-1);a=ev[start:end];b=ev[start+1:end+1]
  s=np.array(a['s'],dtype=np.int64);q=np.full(len(a),3,dtype=np.int64)
  if start==0:q[0]=2
  f=q-s;fm=np.array(a['F'],dtype=np.int64)-f;pm=np.array(a['P'],dtype=np.int64)+2*f
  quota=2*((np.maximum(1,pm//6)+1)//2)
  v=np.where((fm>=3)&(pm>=6),np.minimum.reduce([quota,np.full(len(a),H),fm,pm]),0)
  af=fm-v+np.array(a['cost'],dtype=np.int64);ap=pm-v
  gap=np.array(b['step']-a['step']-1,dtype=np.int64)
  if H==0:
   xf=af-3*gap;xp=ap+6*gap
  else:
   saturated=(ap>=24)&(af>=9*gap)
   xf=af-9*gap;xp=ap.copy()
   for k in np.flatnonzero(~saturated):xf[k],xp[k]=quiet(int(af[k]),int(ap[k]),int(gap[k]),H)
  assert np.all(xf==b['F']) and np.all(xp==b['P'])
  assert np.all(b['tau']-a['tau']==b['s'])
  counters['all_positive_service_ledger_links']+=len(a)
 # Independent full ancestry bitsets, all node chain weights and recorded times.
 ancestors=[1,2];chains=[1,1]
 for i in range(2,len(nodes)):
  a,b=map(int,nodes[i,1:3]);ancestors.append(ancestors[a]|ancestors[b]|(1<<i));chains.append(chains[a]+chains[b])
  assert chains[i]==int(nodes[i,3])
  counters['all_native_node_chains']+=1
 groups=defaultdict_list(nodes)
 ids=sorted(groups);chosen=sorted(set(ids[:64]+[ids[(len(ids)-1)*j//511] for j in range(512)]))
 for step in chosen:
  born=groups[step];idx=int(np.searchsorted(steps,step));e=ev[idx];assert int(e['step'])==step
  W=sorted(int(z) for z in e['W'][:e['s']]);n=int(e['n']);before_pairs={(int(a),int(b)) for a,b in nodes[2:n,1:3]}
  expected=[ab for ab in combinations(W,2) if ab not in before_pairs]
  assert expected==[tuple(map(int,nodes[i,1:3])) for i in born]
  repeat=0;allmask=0
  for u,v in expected:
   mask=ancestors[u]|ancestors[v];allmask|=mask
   while mask:
    b=mask&-mask;z=b.bit_length()-1;mask-=b;repeat+=2*(2*chains[z]-2)
  premium=0
  while allmask:
   b=allmask&-allmask;z=b.bit_length()-1;allmask-=b
   if int(nodes[z,7])>=step:assert int(nodes[z,7])==step;premium+=9*(2*chains[z]-2)
  assert repeat+premium==int(e['cost'])
  counters['independent_late_complete_bursts']+=1

def defaultdict_list(nodes):
 d={}
 for i in range(2,len(nodes)):d.setdefault(int(nodes[i,5]),[]).append(i)
 return d

def exact_local_law(counters):
 # Conditional union subset: no iid replacement or Poisson approximation.
 for L in (3,4,8):
  for c in range(L+1):
   A=set(range(L));B=set(range(c))|set(range(L,2*L-c));U=sorted(A|B);u=len(U);delta=len(A^B)
   for j in range(1,min(3,u)+1):
    diff=[sum((z in A)-(z in B) for z in subset) for subset in combinations(U,j)]
    mean=sum(map(Q,diff))/len(diff);var=sum(Q(x*x) for x in diff)/len(diff)
    expected=Q(j*delta*(u-j),u*(u-1))
    assert mean==0 and var==expected
    counters['exact_union_subset_laws']+=1

def quantum_tests(counters):
 # Independent dense Fourier circuits; no readout formula reused.
 omega=np.zeros(9,complex);omega[[0,4,8]]=1/np.sqrt(3);r=np.arange(3)
 for modulus in (27,81):
  for a,b in ((0,0),(7,2),(19,45),(128,127),(81,27)):
   u=np.kron(np.exp(2j*np.pi*r*a/modulus),np.exp(-2j*np.pi*r*b/modulus))
   for k in range(3):
    state=omega*np.kron(np.exp(2j*np.pi*r*k/3),np.ones(3))*u
    probs=np.zeros(3)
    for x,y in product(range(3),repeat=2):
     vx=np.exp(-2j*np.pi*r*x/3)/np.sqrt(3);vy=np.exp(2j*np.pi*r*y/3)/np.sqrt(3)
     probs[(y-x)%3]+=abs(np.vdot(np.kron(vx,vy),state))**2
    predicted=np.array([(1+2*np.cos(2*np.pi*((a-b)/modulus+(k-word)/3)))**2/9 for word in range(3)])
    assert np.max(abs(probs-predicted))<3e-14
    counters['independent_dense_code_probabilities']+=3

def validate_analysis(counters):
 files=sorted((ROOT/'results').glob('H*.analysis.json.gz'))
 assert len(files)==8
 for p in files:
  with gzip.open(p,'rt') as f:r=json.load(f)
  nodes=np.loadtxt(p.with_name(p.name.replace('.analysis.json.gz','.nodes.tsv')),dtype=np.int64,ndmin=2)
  prefix=WORK/p.name.replace('.analysis.json.gz','')
  events=np.memmap(prefix.with_suffix('.events.bin'),dtype=DT,mode='r')
  fullsteps=np.array(events['step'],dtype=np.int64)
  for target,e in r['epochs'].items():
   if not e['reached']:continue
   cp=e['checkpoint'];assert int(target)<=cp['tau']<=int(target)+2
   for label,co in e['cohorts'].items():
    for root,tokens in co['paths']:
     # Validate backwards, independently of the forward enumerator.
     cur=tokens[-1][0]
     for z,entry in tokens[::-1]:assert cur==z;cur=int(nodes[z,entry+1]);assert 0<nodes[z,7]<=cp['step']
     assert cur==root;counters['selected_complete_paths']+=1
    for pair in co['pair_readouts']:
     pa,pb=(co['paths'][i] for i in (pair['i'],pair['j']));a,b=set(support(pa)),set(support(pb))
     assert len(a)==len(b)==pair['length'] and len(a^b)==pair['delta']
     for mode in ('global_windows','local_windows'):
      for target,read in pair[mode].items():
       if read['complete']:
        actual=read['elapsed_tau'] if mode=='global_windows' else read['union_activity']
        assert int(target)<=actual<=int(target)+2
       d=read['difference'];assert d==read['N_A']-read['N_B']
       assert read['union_activity']==read['private_activity']+read['shared_activity']
       assert read['N_A']+read['N_B']==read['private_activity']+2*read['shared_activity']
       for Q0 in (27,81):
        expected=1/3+4/9*np.cos(2*np.pi*(d%Q0)/Q0)+2/9*np.cos(4*np.pi*(d%Q0)/Q0)
        assert abs(expected-read['correct_word'][str(Q0)])<2e-14
       if pair['delta']==0:assert d==0
       counters['readout_records_checked']+=1
    if co['pair_readouts']:
     pair=co['pair_readouts'][0]
     pa,pb=(co['paths'][i] for i in (pair['i'],pair['j']))
     a,b=set(support(pa)),set(support(pb))
     start=int(np.searchsorted(fullsteps,cp['step'],side='right'))
     for mode,window in (('global_windows','8192'),('local_windows','32'),('local_windows','128')):
      rd=pair[mode][window];end=int(np.searchsorted(fullsteps,cp['step']+rd['elapsed_steps'],side='right'))
      W=events['W'][start:end]
      na=int(np.isin(W,list(a)).sum());nb=int(np.isin(W,list(b)).sum());un=int(np.isin(W,list(a|b)).sum())
      assert (na,nb,un)==(rd['N_A'],rd['N_B'],rd['union_activity'])
      counters['independent_raw_count_windows']+=1
    for cl in co['collective_availability'].values():
     for ex in cl.get('examples',{}).values():
      pp=ex['paths'];supports=[set(support(p)) for p in pp]
      assert len({tuple(sorted(s)) for s in supports})==4
      balance=Counter()
      for sign,s in zip((1,1,-1,-1),supports):
       for z in s:balance[z]+=sign
      assert not any(balance.values()) and ex['sum_difference']==0
      tokens=Counter()
      for sign,(_,ts) in zip((1,1,-1,-1),pp):
       for t in ts:tokens[tuple(t)]+=sign
      assert ex['entry_token_balance']==(not any(tokens.values()))
      if ex['entry_token_balance']:assert ex['entry_sensitive_code_phases']==[0,0,0]
      counters['four_distinct_path_certificates']+=1

if __name__=='__main__':
 start=time.perf_counter();c=Counter()
 for H in (0,6):
  for s in range(20350920,20350924):
   validate_history(WORK/f'H{H}_s{s}',c);print('validated native',H,s,flush=True)
 exact_local_law(c);quantum_tests(c);validate_analysis(c)
 report={'counts':dict(c),'seconds':time.perf_counter()-start,'all_checks_passed':True}
 (ROOT/'DCU_Mass_Cell_35_TESTS.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
