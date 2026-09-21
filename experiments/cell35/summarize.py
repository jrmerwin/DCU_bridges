from pathlib import Path
import gzip,json,numpy as np,math
ROOT=Path(__file__).resolve().parent
with gzip.open(ROOT/'results/DCU_Mass_Cell_35_RESULTS.json.gz','rt') as f:r=json.load(f)
rep=json.loads((ROOT/'replication/replication_results.json').read_text())
print('MAIN iterations',r['total_native_iterations'],'tau sum',r['sum_of_maintenance_across_worlds'])
print('EPOCH TABLE')
for target in r['epoch_targets']:
 cps=[w['epochs'][str(target)]['checkpoint'] for w in r['trajectories']]
 print(target, 'n',min(c['n'] for c in cps),max(c['n'] for c in cps),'t',min(c['step'] for c in cps),max(c['step'] for c in cps),'F/n',min(c['queue_ratio'] for c in cps),max(c['queue_ratio'] for c in cps))
print('PRIMARY new selected L3 delta2, completed-only plus censor counts')
for H in (0,6):
 for target in r['epoch_targets']:
  q=r['summary'][str(H)][str(target)]['cohorts']['short_L3']['strata']['L3_delta2']['modes']
  rows=[]
  for mode,w in [('global_windows','8192'),('local_windows','32'),('local_windows','128')]:
   z=q[mode][w];rows.append((z['correct_word']['27']['mean'],z['correct_word']['27']['se'],z['censored_pairs'],z['selected_pairs']))
  print(H,target,rows)
print('REPLICATES conditional on four starting worlds')
rep_summary={}
for H in (0,6):
 rep_summary[str(H)]={}
 for age in (128,4194304):
  rep_summary[str(H)][str(age)]={}
  tasks=[t for t in rep['tasks'] if t['H']==H and t['epoch']==age]
  for mode in ('global8192','local32'):
   arrays=[np.array([q['correct_word']['27'] for q in t['readouts'] if q['mode']==mode]) for t in tasks]
   na=[np.array([q['union_services'] for q in t['readouts'] if q['mode']==mode]) for t in tasks]
   nt=[np.array([q['elapsed_tau'] for q in t['readouts'] if q['mode']==mode]) for t in tasks]
   ns=[np.array([q['elapsed_steps'] for q in t['readouts'] if q['mode']==mode]) for t in tasks]
   d2=[np.array([q['difference']**2 for q in t['readouts'] if q['mode']==mode]) for t in tasks]
   n=len(arrays);mean=float(np.mean([a.mean() for a in arrays]));se=float(np.sqrt(sum(a.var(ddof=1)/len(a) for a in arrays))/n)
   out=dict(mean=mean,conditional_MC_SE=se,world_mean_SE=float(np.std([a.mean() for a in arrays],ddof=1)/np.sqrt(n)),
    per_world_means=[float(a.mean()) for a in arrays],mean_activity=float(np.mean(na)),mean_elapsed_tau=float(np.mean(nt)),mean_elapsed_steps=float(np.mean(ns)),mean_D_squared=float(np.mean(d2)),
    completed=sum(q['complete'] for t in tasks for q in t['readouts'] if q['mode']==mode),total=sum(q['mode']==mode for t in tasks for q in t['readouts']))
   rep_summary[str(H)][str(age)][mode]=out;print(H,age,mode,out)
 # Paired MC difference, conditional on same source-world panel.
 dif=[]
 for seed in range(20350920,20350924):
  a=next(t for t in rep['tasks'] if (t['H'],t['epoch'],t['source_seed'])==(H,128,seed));b=next(t for t in rep['tasks'] if (t['H'],t['epoch'],t['source_seed'])==(H,4194304,seed))
  aa=np.array([q['correct_word']['27'] for q in a['readouts'] if q['mode']=='local32']);bb=np.array([q['correct_word']['27'] for q in b['readouts'] if q['mode']=='local32']);dif.append(bb-aa)
 rep_summary[str(H)]['late_minus_early_local32']=dict(mean=float(np.mean(dif)),conditional_paired_MC_SE=float(np.sqrt(sum(a.var(ddof=1)/len(a) for a in dif))/len(dif)))
print('COLLECTIVE CENSUS')
for target in r['epoch_targets']:
 counts=[];entries=[];trunc=0;recent=[]
 for w in r['trajectories']:
  ep=w['epochs'][str(target)];ns=ne=0
  for label in ('short_L3','short_L4'):
   for x in ep['cohorts'][label]['collective_availability'].values():
    if not x['complete']:trunc+=1;continue
    ns+=x['register_balance_relations'];ne+=x['entry_token_balance_relations']
  counts.append(ns);entries.append(ne)
  co=ep['cohorts']['recent_full'];recent.extend([int(l) for l in co['path_lengths']])
 print(target,'register',counts,'entry',entries,'incomplete panels',trunc,'recent lengths',min(recent),max(recent))
allrows=[p for w in r['trajectories'] for ep in w['epochs'].values() if ep['reached'] for co in ep['cohorts'].values() for p in co['pair_readouts']]
paths=sum(co['path_count'] for w in r['trajectories'] for ep in w['epochs'].values() if ep['reached'] for co in ep['cohorts'].values())
print('SELECTED PAIRS',len(allrows),'PATH APPEARANCES',paths)
print('eligibility pairs',sum(z['eligible_pair_count'] for w in r['trajectories'] for ep in w['epochs'].values() if ep['reached'] for co in ep['cohorts'].values() for z in co['strata'].values()))
print('fork_steps',sum(max(q['elapsed_steps'] for q in t['readouts'] if q['rep']==i) for t in rep['tasks'] for i in range(32)))
(ROOT/'replication/replication_summary.json').write_text(json.dumps(rep_summary,indent=2))
