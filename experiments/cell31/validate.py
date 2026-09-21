from pathlib import Path
from collections import Counter
from fractions import Fraction
from itertools import combinations
from random import Random
from math import sqrt
import gzip,json,sys,statistics,time
from support import helpers,clean

ROOT=Path(__file__).resolve().parent
start=time.monotonic()
env,h=helpers()
result=json.load(gzip.open(ROOT/'DCU_Mass_Cell_31_RESULTS.json.gz','rt'))
counts=Counter()

# Reconstruct every private frame using only its private cards and permitted
# readout operations; crosswalk and full parentage are validation-only.
for cohort,states in result['cohorts'].items():
 for ename,cases in states.items():
  for arm,panel in cases.items():
   for hist in panel['histories']:
    parentage=hist['final_parents'];anc=[]
    for z,pair in enumerate(parentage):
     anc.append({z} if pair is None else {z}|anc[pair[0]]|anc[pair[1]])
    for policy,variants in hist['policies'].items():
     for b,proto in variants.items():
      for role,rd in enumerate(proto['readers'],1):
       back={v:int(k) for k,v in rd['auditor_handle_map'].items()}
       assert len(back)==len(rd['auditor_handle_map'])
       known={a for a,p in rd['initial_local_cards']}
       graph={a:None if p is None else tuple(p) for a,p in rd['initial_local_cards']}
       active=rd['initial_active_handle'];track=[back[active]]
       rng=Random(hist['seed']*1000003+1009*role+rd['origin'])
       hops=ties=0
       for label,pair in graph.items():
        g=back[label]
        assert (None if pair is None else tuple(back[p] for p in pair)) == (None if parentage[g] is None else tuple(parentage[g]))
        if pair:assert set(pair)<=known
       for t,card in enumerate(rd['local_log'],1):
        assert set(card)=={'local_tick','carrier_before','carrier_after','new_children','partners','newly_known_parent_cards'}
        assert card['local_tick']==t and card['carrier_before']==active
        old_known=set(known)
        gained={a for a,p in card['newly_known_parent_cards']}
        assert not (gained&known)
        known|=gained
        for a,pair in card['newly_known_parent_cards']:
         g=back[a]; pp=None if pair is None else tuple(pair)
         assert (None if pp is None else tuple(back[p] for p in pp))==(None if parentage[g] is None else tuple(parentage[g]))
         if pp:assert set(pp)<=known
         graph[a]=pp
        own=card['new_children'];old_id=back[active]
        expected_gained=set()
        for child in own:
         assert active in graph[child]
         expected_gained|=anc[back[child]]
        assert {back[a] for a in gained}==expected_gained-{back[a] for a in old_known}
        assert card['partners']==[next(z for z in graph[c] if z!=active) for c in own]
        after=active
        if policy=='follow_unique' and len(own)==1:after=own[0]
        if policy=='follow_any' and own:
         after=own[0] if len(own)==1 else own[rng.randrange(len(own))]
        assert card['carrier_after']==after
        ties+=len(own)>1
        if after!=active:
         assert old_id in parentage[back[after]] and after in own
         track.append(back[after]);hops+=1
        active=after;counts['private_event_cards']+=1
       assert track==rd['auditor_track'] and back[active]==rd['active']
       assert hops==rd['hops'] and ties==rd['ambiguous_services']
       assert rd['ticks']==len(rd['local_log']) and hops<=rd['ticks']
       assert {back[a] for a in known}==set(rd['final_known'])
       counts['private_reader_frames']+=1
      rda,rdb=proto['readers']
      for receipt in proto['all_receipt_cards']:
       z=receipt['child'];b0=receipt['active_before'][1]
       assert b0 in parentage[z] and panel['marker'] in anc[z]
       c=rdb['local_log'][receipt['receiver_tick']-1]
       assert rdb['auditor_handle_map'][str(z)] in c['new_children']
       assert c['carrier_before']==rdb['auditor_handle_map'][str(b0)]
       assert receipt['jointly_witnessed_with_source']==(receipt['active_before'][0] in parentage[z])
       counts['receipt_certificates']+=1
      returned=proto['returned_remote_receipt']
      if returned:
       r=returned['receipt_child'];w=returned['child']
       assert r in anc[w] and r<w
       rr=next(c for c in proto['all_receipt_cards'] if c['child']==r)
       assert not rr['jointly_witnessed_with_source'] and rr['auditor_step']<returned['auditor_step']
       lc=rda['local_log'][returned['source_tick']-1]
       assert rda['auditor_handle_map'][str(w)] in lc['new_children']
       assert returned['active_before'][0] in parentage[w]
       counts['returned_certificates']+=1
      merge=proto['coalescence']
      if merge:
       z=merge['shared_carrier']
       assert set(parentage[z])==set(merge['active_before'])
       assert rda['active']==rdb['active']==z and proto['stop_step']==merge['auditor_step']
       counts['coalescences']+=1
      else:assert rda['active']!=rdb['active'] and proto['stop_step']==4096
      counts['pair_protocols']+=1
   # Independent arithmetic for complete finite ensemble summaries.
   for policy,variants in panel['summary']['by_policy'].items():
    for b,s in variants.items():
     rows=[v['policies'][policy][b] for v in panel['histories']]
     D=[r['readers'][0]['ticks']-r['readers'][1]['ticks'] for r in rows]
     expect=statistics.mean(D)
     se=statistics.stdev(D)/sqrt(len(D))
     assert abs(expect-s['stopped_clock_difference']['mean'])<1e-14
     assert abs(se-s['stopped_clock_difference']['se'])<1e-14
     for k,field in [('all_receipts','receipt'),('all_remote_returns','returned_remote_receipt'),('coalescences','coalescence')]:
      assert s[k]==sum(r[field] is not None for r in rows)
     counts['summary_rows']+=1
print('PASS private logs/choices/certificates:',dict(counts))

class Independent:
 def __init__(self,parents):
  self.parents=[None,None];self.bits=[1,2];self.ch=[1,1];self.R=0;self.pairs={}
  self.add([tuple(p) for p in parents[2:]],sequential_prepare=True)
 def add(self,pairs,sequential_prepare=False):
  if sequential_prepare:
   for p in pairs:self.add([p])
   return
  cost=0;recorded=self.R
  # Independent sequential first/repeat formulation, instead of event hit counts.
  for a,b in pairs:
   bits=self.bits[a]|self.bits[b]
   while bits:
    bit=bits&-bits;z=bit.bit_length()-1;bits-=bit
    cost+=(2 if recorded&bit else 11)*(2*self.ch[z]-2)
    recorded|=bit
  born=[]
  for a,b in pairs:
   z=len(self.parents);self.parents.append((a,b));self.bits.append((1<<z)|self.bits[a]|self.bits[b]);self.ch.append(self.ch[a]+self.ch[b]);self.pairs[a,b]=z;born.append(z)
  self.R=recorded
  return born,cost
 def step(self,F,P,rng,H):
  n=len(self.parents);total=F+n;q=min(3,total)
  selected=sorted(x for x in rng.sample(range(total),q) if x<n)
  f=q-len(selected);fm=F-f;pm=P+2*f
  quota=max(1,pm//6)
  if quota%2:quota+=1
  v=min(quota,H,fm,pm) if fm>=3 and pm>=6 else 0
  pairs=[p for p in combinations(selected,2) if p not in self.pairs]
  born,cost=self.add(pairs)
  return fm-v+cost,pm-v,dict(n=n,F=F,P=P,total=total,q=q,served=tuple(selected),forced_served=f,relief_removed=v,new_pairs=tuple(pairs),born=tuple(born),cost=cost)

step=h['native_step']
for cohort,states in result['cohorts'].items():
 for ename,cases in states.items():
  for arm,panel in cases.items():
   for hist in (panel['histories'][0],panel['histories'][-1]):
    independent=Independent(panel['initial_parentage'])
    st=h['construct'](panel['initial_parentage'])
    F=P=None; F,P=panel['initial_F'],panel['initial_P'];fi,pi=F,P
    r1,r2=Random(hist['seed']),Random(hist['seed'])
    for t in range(4096):
     F,P,e=step(st,F,P,r1,panel['H'])
     fi,pi,ei=independent.step(fi,pi,r2,panel['H'])
     assert e==ei and (F,P)==(fi,pi)
     assert st.parents==independent.parents and st.chains==independent.ch
     assert st.recorded=={z for z in range(len(st)) if independent.R>>z&1}
     if hist['global_auditor_trace']:
      archived=hist['global_auditor_trace'][t]
      assert clean(dict(e,step=t+1,F_after=F,P_after=P))==archived
     counts['independent_native_transitions']+=1
    assert (len(st),F,P)==(hist['final_n'],hist['final_F'],hist['final_P'])
    assert clean(st.parents)==hist['final_parents']
    counts['observer_free_replays']+=1
print('PASS independent native transitions and observer-free final states:',dict(counts))
exact=h['exact_audits']()
assert clean(exact)==result['exact']
counts['exact_child_laws']=exact['conditional_child_law_cases']
counts['merger_laws']=len(exact['two_childless_carrier_merger_checks'])
ce=exact['identity_counterexample']
assert ce['source_local_log_identical'] and ce['tracked_receipt_exists']==(False,True)
assert ce['tracked_return_exists']==(False,True)
# Rebuild every fork-counterexample transition with independent native code, too.
initial=result['cohorts']['replayed']['no_relief']['recorded_pair']
st=Independent(ce['initial_parentage']);F,P=initial['initial_F'],initial['initial_P']
class Draw:
 def __init__(self,e):self.e=e
 def sample(self,population,k):
  return list(self.e['served'])+list(range(self.e['n'],self.e['n']+self.e['forced_served']))
for e in ce['conditioned_native_events']:
 F,P,ee=st.step(F,P,Draw(e),0)
 assert clean(ee)==clean({k:v for k,v in e.items() if k not in ('step','F_after','P_after')})
 counts['conditioned_fork_events']+=1
assert clean(st.parents)==clean(ce['final_parentage'])
print('PASS exact sampling, merger law, and identical-source-log identity counterexample.')
print('ALL CHECKS:',json.dumps(dict(counts),sort_keys=True))
print('Validation elapsed seconds',time.monotonic()-start)
print('No statistical fluctuation was used as a pass/fail filter; no physical SI, mass, or spacetime claim is checked.')
