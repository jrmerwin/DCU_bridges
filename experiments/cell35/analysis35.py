from pathlib import Path
from itertools import combinations
from collections import defaultdict, Counter
from random import Random
import hashlib, json, gzip, math, time
import numpy as np

EPOCHS=(128,1024,8192,65536,524288,1048576,4194304)
GLOBAL_WINDOWS=(128,512,2048,8192)
LOCAL_WINDOWS=(8,32,128)
DT=np.dtype([('step','<u8'),('tau','<u8'),('n','<u8'),('F','<u8'),('P','<u8'),('cost','<u8'),('s','<u4'),('W','<u4',(3,))])
assert DT.itemsize==64

def jsonable(v):
 if isinstance(v,np.ndarray):return v.tolist()
 if isinstance(v,np.generic):return v.item()
 if isinstance(v,dict):return {str(k):jsonable(x) for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [jsonable(x) for x in v]
 return v

def selector(*parts):
 return Random(int.from_bytes(hashlib.sha256('|'.join(map(str,parts)).encode()).digest()[:8],'big'))

def short_paths(nodes,n,step,maxlength=4):
 cache={0:[(0,())],1:[(1,())]}
 out={k:[] for k in (3,4)}
 for z in range(2,n):
  cache[z]=[(root,tok+((z,e),)) for e,p in enumerate(nodes[z,1:3]) for root,tok in cache[int(p)] if len(tok)<maxlength]
  if 0<int(nodes[z,7])<=step:
   for path in cache[z]:
    if len(path[1]) in out:out[len(path[1])].append(path)
 return {k:sorted(v) for k,v in out.items()}

def recent_paths(nodes,n,step,target,seed,H):
 candidates=[z for z in range(2,n) if 0<int(nodes[z,7])<=step and target//8<int(nodes[z,6])]
 rng=selector('recent',seed,H,target)
 chosen=sorted(rng.sample(candidates,min(64,len(candidates))))
 raw=[]
 for z in chosen:
  for _ in range(2):
   k=rng.randrange(int(nodes[z,3]));cur=z;tokens=[]
   while cur>=2:
    a,b=map(int,nodes[cur,1:3]);ca=int(nodes[a,3]);e=0 if k<ca else 1
    if e:k-=ca
    tokens.append((cur,e));cur=(a,b)[e]
   raw.append((cur,tuple(reversed(tokens))))
 return sorted(set(raw)),dict(eligible_recent_endpoints=len(candidates),selected_endpoints=chosen,raw_path_draws=len(raw),distinct_paths=len(set(raw)))

def support(path):return tuple(z for z,e in path[1])

def census_and_pairs(paths,cap,seed,H,target,label):
 # Same lexicographic population and RNG selection as materializing every pair,
 # but keep a byte-valued overlap table instead of millions of Python tuples.
 bylength=defaultdict(list)
 for i,p in enumerate(paths):bylength[len(p[1])].append(i)
 selected=[];summaries={}
 for length,ids in sorted(bylength.items()):
  if len(ids)<2:continue
  occurrence=defaultdict(list)
  for i,idx in enumerate(ids):
   for z in support(paths[idx]):occurrence[z].append(i)
  occurrence={z:np.asarray(v,dtype=np.int64) for z,v in occurrence.items()}
  rows=[];counts=[]
  for i,idx in enumerate(ids[:-1]):
   shared=np.zeros(len(ids)-i-1,dtype=np.int16)
   for z in support(paths[idx]):
    occ=occurrence[z];occ=occ[int(np.searchsorted(occ,i+1)):]
    shared[occ-i-1]+=1
   unique=(length-shared).astype(np.uint8 if length<256 else np.uint16)
   rows.append(unique);counts.append(np.bincount(unique,minlength=length+1))
  cumulative=np.cumsum(np.asarray(counts,dtype=np.int64),axis=0)
  for m,total in enumerate(cumulative[-1]):
   total=int(total)
   if not total:continue
   delta=2*m;rng=selector('pairs',seed,H,target,label,length,delta)
   positions=sorted(rng.sample(range(total),min(cap,total)));chosen=[]
   for position in positions:
    i=int(np.searchsorted(cumulative[:,m],position,side='right'))
    prior=int(cumulative[i-1,m]) if i else 0
    j=i+1+int(np.flatnonzero(rows[i]==m)[position-prior])
    chosen.append((ids[i],ids[j]))
   key=f'L{length}_delta{delta}'
   summaries[key]=dict(length=length,delta=delta,shared=length-delta//2,eligible_pair_count=total,selected_pair_count=len(chosen))
   for i,j in sorted(chosen):selected.append(dict(stratum=key,i=i,j=j,length=length,delta=delta))
 return selected,summaries

def collective_census(paths,max_supports=2048):
 # Fixed lexicographically first path per support. No choice based on quantum outputs.
 bylen=defaultdict(dict)
 for p in paths:bylen[len(p[1])].setdefault(support(p),p)
 out={}
 for length,unique in sorted(bylen.items()):
  ps=[unique[s] for s in sorted(unique)];n=len(ps)
  row=dict(path_length=length,distinct_supports=n,complete=n<=max_supports)
  if n>max_supports:
   row['reason']='Collective census resource cap, not evidence of absence';out[str(length)]=row;continue
  reggroups=defaultdict(list);tokgroups=defaultdict(list)
  for i,j in combinations(range(n),2):
   a,b=set(support(ps[i])),set(support(ps[j]));key=(tuple(sorted(a|b)),tuple(sorted(a&b)))
   reggroups[key].append((i,j))
   a,b=set(ps[i][1]),set(ps[j][1]);key=(tuple(sorted(a|b)),tuple(sorted(a&b)))
   tokgroups[key].append((i,j))
  row['register_balance_relations']=sum(math.comb(len(v),2) for v in reggroups.values())
  row['entry_token_balance_relations']=sum(math.comb(len(v),2) for v in tokgroups.values())
  row['examples']={}
  for name,groups in (('register',reggroups),('entry_token',tokgroups)):
   group=next((v for k,v in sorted(groups.items()) if len(v)>1),None)
   if group:
    ids=group[0]+group[1]
    assert len(set(ids))==4
    pp=[ps[i] for i in ids]
    bal=Counter(support(pp[0]));bal.update(support(pp[1]));bal.subtract(support(pp[2]));bal.subtract(support(pp[3]))
    assert all(v==0 for v in bal.values())
    tokenbal=Counter(pp[0][1]);tokenbal.update(pp[1][1]);tokenbal.subtract(pp[2][1]);tokenbal.subtract(pp[3][1])
    row['examples'][name]=dict(paths=pp,register_balance=True,entry_token_balance=all(v==0 for v in tokenbal.values()))
  out[str(length)]=row
 return out

def service_index(events,n):
 ids=[np.asarray(events['W'][:,0],dtype=np.uint32)];ts=[np.arange(len(events),dtype=np.uint32)]
 for k in (1,2):
  mask=events['s']>k;ids.append(np.asarray(events['W'][mask,k],dtype=np.uint32));ts.append(np.flatnonzero(mask).astype(np.uint32))
 ids=np.concatenate(ids);ts=np.concatenate(ts)
 assert len(ids)==int(events[-1]['tau'])
 order=np.argsort(ids,kind='stable');counts=np.bincount(ids,minlength=n);offset=np.r_[0,np.cumsum(counts)]
 times=ts[order];del ids,ts,order
 return [np.sort(times[offset[z]:offset[z+1]]) for z in range(n)]

def checkpoint(events,nodes,idx,H):
 e=events[idx];step,tau=int(e['step']),int(e['tau']);n=int(np.searchsorted(nodes[:,5],step,side='right'))
 q=min(3,int(e['n']+e['F']));f=q-int(e['s']);fm=int(e['F'])-f;pm=int(e['P'])+2*f
 quota=2*((max(1,pm//6)+1)//2);v=min(quota,H,fm,pm) if fm>=3 and pm>=6 else 0
 F=fm-v+int(e['cost']);P=pm-v
 return dict(step=step,tau=tau,n=n,F=F,P=P,queue_ratio=F/n,recorded_objects=int(np.sum((nodes[:n,7]>0)&(nodes[:n,7]<=step))))

def pair_readouts(row,paths,idx,events,services,taus):
 a,b=map(set,(support(paths[row['i']]),support(paths[row['j']])))
 union=sorted(a|b);u=len(union);delta=len(a^b)
 phase_start=[int(np.searchsorted(services[z],idx,side='right')) for z in union]
 seqs=[services[z][k:] for z,k in zip(union,phase_start)]
 merged=np.concatenate(seqs);uniq,j=np.unique(merged,return_counts=True)
 total=np.cumsum(j,dtype=np.int64)
 bracket=np.cumsum(j*(u-j)/(u-1),dtype=float) if u>1 else np.zeros(len(j))
 stau,sstep=int(events[idx]['tau']),int(events[idx]['step'])
 def record(endidx,mode,requested):
  if endidx is None:
   endidx=len(events)-1;complete=False
  else:complete=True
  per=[int(np.searchsorted(services[z],endidx,side='right'))-k for z,k in zip(union,phase_start)]
  nd=dict(zip(union,per));na=sum(nd[z] for z in a);nb=sum(nd[z] for z in b);d=na-nb;end=events[endidx]
  k=int(np.searchsorted(uniq,endidx,side='right'))-1
  activity=int(total[k]) if k>=0 else 0;breal=float(bracket[k]) if k>=0 else 0
  assert activity==sum(per)
  return dict(complete=complete,window_type=mode,requested=requested,N_A=na,N_B=nb,difference=d,
   union_activity=activity,private_activity=sum(nd[z] for z in a^b),shared_activity=sum(nd[z] for z in a&b),
   multi_hit_bursts=int(np.sum(j[:k+1]>1)),union_hit_bursts=k+1,
   predictable_variance=delta/u*breal,squared_difference=d*d,
   elapsed_tau=int(end['tau'])-stau,elapsed_steps=int(end['step'])-sstep,
   correct_word={str(Q):float((1+2*math.cos(2*math.pi*(d%Q)/Q))**2/9) for Q in (27,81)})
 out=dict(row,union_size=u,unique_supports=delta>0,global_windows={},local_windows={})
 for w in GLOBAL_WINDOWS:
  endidx=int(np.searchsorted(taus,stau+w,side='left'));endidx=None if endidx>=len(events) else endidx
  out['global_windows'][str(w)]=record(endidx,'global_tau',w)
 for m in LOCAL_WINDOWS:
  k=int(np.searchsorted(total,m,side='left'));endidx=int(uniq[k]) if k<len(uniq) else None
  out['local_windows'][str(m)]=record(endidx,'union_service',m)
 return out

def analyze_history(prefix,outdir):
 started=time.perf_counter();meta=json.loads(prefix.with_suffix('.meta.json').read_text());H=meta['H'];seed=meta['seed']
 events=np.memmap(prefix.with_suffix('.events.bin'),dtype=DT,mode='r');nodes=np.loadtxt(prefix.with_suffix('.nodes.tsv'),dtype=np.int64,ndmin=2)
 assert len(nodes)==meta['n'] and int(events[-1]['tau'])==meta['tau'];assert np.all(np.diff(events['tau'].astype(np.int64))>0)
 assert all(tuple(nodes[i,1:3])==(-1,-1) for i in (0,1))
 assert np.all(nodes[2:,1]<nodes[2:,2]) and np.all(nodes[2:,2]<np.arange(2,len(nodes)))
 services=service_index(events,len(nodes));ages={};taus=np.array(events['tau'],dtype=np.int64)
 for target in EPOCHS:
  idx=int(np.searchsorted(taus,target,side='left'))
  if idx>=len(events):ages[str(target)]=dict(reached=False);continue
  cp=checkpoint(events,nodes,idx,H);cp.update(target_tau=target,overshoot=cp['tau']-target)
  assert 0<=cp['overshoot']<=2
  short=short_paths(nodes,cp['n'],cp['step']);recent,rs=recent_paths(nodes,cp['n'],cp['step'],target,seed,H)
  cohorts={}
  for name,paths,cap in [('short_L3',short[3],16),('short_L4',short[4],16),('recent_full',recent,4)]:
   for root,tokens in paths:
    p=root
    for z,e in tokens:
     assert z<cp['n'] and int(nodes[z,7])<=cp['step'] and int(nodes[z,7])>0 and int(nodes[z,e+1])==p;p=z
   sel,counts=census_and_pairs(paths,cap,seed,H,target,name)
   rr=[pair_readouts(r,paths,idx,events,services,taus) for r in sel]
   coll=collective_census(paths)
   for v in coll.values():
    for ex in v.get('examples',{}).values():
     endidx=int(np.searchsorted(taus,cp['tau']+8192,side='left'));endidx=min(endidx,len(events)-1)
     regs=sorted(set(z for pp in ex['paths'] for z,e in pp[1]));nd={z:int(np.searchsorted(services[z],endidx,side='right')-np.searchsorted(services[z],idx,side='right')) for z in regs}
     cnt=[sum(nd[z] for z,e in pp[1]) for pp in ex['paths']]
     ex['counts_after_global_8192']=cnt;ex['sum_difference']=cnt[0]+cnt[1]-cnt[2]-cnt[3];assert ex['sum_difference']==0
     entry_phase=[]
     for r in range(3):
      val=0
      for sign,pp in zip((1,1,-1,-1),ex['paths']):
       val+=sign*sum(nd[z]*((1,0,2)[r] if e else r) for z,e in pp[1])
      entry_phase.append(val)
     ex['entry_sensitive_code_phases']=entry_phase
     if ex['entry_token_balance']:assert entry_phase==[0,0,0]
   cohorts[name]=dict(path_count=len(paths),distinct_support_count=len(set(map(support,paths))),path_lengths=dict(Counter(len(p[1]) for p in paths)),
    paths=paths,selection_cap_per_stratum=cap,strata=counts,pair_readouts=rr,collective_availability=coll)
  ages[str(target)]=dict(reached=True,checkpoint=cp,recent_selection=rs,cohorts=cohorts)
  print('ANALYZED',H,seed,target,cp['n'],{k:v['path_count'] for k,v in cohorts.items()},flush=True)
 result=dict(metadata=meta,epochs=ages,event_count=len(events),analysis_seconds=time.perf_counter()-started)
 outdir.mkdir(exist_ok=True,parents=True)
 with gzip.open(outdir/(prefix.name+'.analysis.json.gz'),'wt') as f:json.dump(jsonable(result),f,separators=(',',':'))
 # Complete parentage permits independent path and selected-burst validation; no endpoint labels are particles.
 import shutil
 shutil.copy2(prefix.with_suffix('.nodes.tsv'),outdir/(prefix.name+'.nodes.tsv'))
 shutil.copy2(prefix.with_suffix('.audit.tsv'),outdir/(prefix.name+'.audit.tsv'))
 return result

def aggregate(results):
 out={}
 for H in (0,6):
  worlds=[r for r in results if r['metadata']['H']==H];erows={}
  for target in EPOCHS:
   ews=[r['epochs'][str(target)] for r in worlds if r['epochs'][str(target)]['reached']]
   row=dict(reached_worlds=len(ews),checkpoints=[x['checkpoint'] for x in ews],cohorts={})
   for label in ('short_L3','short_L4','recent_full'):
    tags=sorted({k for e in ews for k in e['cohorts'][label]['strata']})
    ss={}
    for tag in tags:
     info=next(e['cohorts'][label]['strata'][tag] for e in ews if tag in e['cohorts'][label]['strata'])
     z=dict(length=info['length'],delta=info['delta'],modes={})
     for mode,windows in (('global_windows',GLOBAL_WINDOWS),('local_windows',LOCAL_WINDOWS)):
      z['modes'][mode]={}
      for window in windows:
       wrs=[];npairs=0;ncens=0
       for e in ews:
        rs=[p[mode][str(window)] for p in e['cohorts'][label]['pair_readouts'] if p['stratum']==tag]
        if not rs:continue
        npairs+=len(rs);ncens+=sum(not p['complete'] for p in rs)
        # Complete-only quantities are explicitly separate; bounds include all selected pairs.
        complete=[p for p in rs if p['complete']]
        rec=dict(selected=len(rs),complete=len(complete),
         correct_lower={str(Q):sum(p['correct_word'][str(Q)] for p in complete)/len(rs) for Q in (27,81)},
         correct_upper={str(Q):(sum(p['correct_word'][str(Q)] for p in complete)+len(rs)-len(complete))/len(rs) for Q in (27,81)})
        if complete:
         rec.update({key:sum(p[key] for p in complete)/len(complete) for key in ('elapsed_tau','elapsed_steps','union_activity','private_activity','shared_activity','squared_difference','predictable_variance','multi_hit_bursts')})
         rec['correct_word']={str(Q):sum(p['correct_word'][str(Q)] for p in complete)/len(complete) for Q in (27,81)}
        wrs.append(rec)
       def meanse(v):
        return dict(n=len(v),mean=float(np.mean(v)) if v else None,se=float(np.std(v,ddof=1)/math.sqrt(len(v))) if len(v)>1 else None)
       summary=dict(selected_pairs=npairs,censored_pairs=ncens,per_world=wrs,world_count=len(wrs))
       for k in ('elapsed_tau','elapsed_steps','union_activity','private_activity','shared_activity','squared_difference','predictable_variance','multi_hit_bursts'):
        summary[k]=meanse([r[k] for r in wrs if k in r])
       for k in ('correct_word','correct_lower','correct_upper'):
        summary[k]={str(Q):meanse([r[k][str(Q)] for r in wrs if k in r]) for Q in (27,81)}
       z['modes'][mode][str(window)]=summary
     ss[tag]=z
    row['cohorts'][label]=dict(strata=ss,path_counts=[e['cohorts'][label]['path_count'] for e in ews])
   erows[str(target)]=row
  out[str(H)]=erows
 return out

def analyse_campaign(workdir,outdir):
 workdir=Path(workdir);outdir=Path(outdir);results=[]
 for H in (0,6):
  for seed in range(20350920,20350924):
   prefix=workdir/f'H{H}_s{seed}';cached=outdir/(prefix.name+'.analysis.json.gz')
   if cached.exists():
    with gzip.open(cached,'rt') as f:results.append(json.load(f))
    print('REUSING COMPLETE ANALYSIS',H,seed,flush=True)
   else:results.append(analyze_history(prefix,outdir))
 agg=aggregate(results)
 result=dict(protocol='Cell35 fixed epoch and activity-matched persistence; separate four-path extension',
  epoch_targets=EPOCHS,global_windows=GLOBAL_WINDOWS,union_service_windows=LOCAL_WINDOWS,
  primary_phase_modulus=27,control_phase_modulus=81,source_parameters=dict(Gamma=3,m=0,H=(0,6)),
  trajectories=results,summary=agg,total_native_iterations=sum(r['metadata']['steps'] for r in results),
  sum_of_maintenance_across_worlds=sum(r['metadata']['tau'] for r in results),
  chronological_iterations_are_not_ticks=True,physical_maturity_epoch_calibrated=False,
  native_counting_law_changed=False,four_path_fanout_is_separate_extension=True,
  quantum_phase_protection_is_not_a_necessary_and_sufficient_matter_test=True)
 with gzip.open(outdir/'DCU_Mass_Cell_35_RESULTS.json.gz','wt') as f:json.dump(jsonable(result),f,separators=(',',':'))
 return result
