# CELL 35 — later-epoch persistence and activity-matched controls.
# SELF-CONTAINED: independent native trajectories; does not modify earlier cells.
# Python 3.10+, NumPy and a C++17 compiler (clang++/g++/c++) are required for a full run.
# ~2.3 GB scratch for full event logs. No private neutron-search files are accessed.
# The native law is unchanged; the quantum wiring remains a declared attachment.
# Paired record storage is not a necessary-and-sufficient definition of matter.

def _run_dcu_mass_35(workdir=None, workers=2, include_replication=True):

    from concurrent.futures import ThreadPoolExecutor
    import subprocess, shutil, sys
    NATIVE_CPP = r'''
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using U=uint64_t;
struct RNG {
 U x; explicit RNG(U seed):x(seed){}
 U next(){ U z=(x+=0x9e3779b97f4a7c15ULL); z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;
 z=(z^(z>>27))*0x94d049bb133111ebULL;return z^(z>>31); }
 U below(U n){if(!n)throw std::runtime_error("zero range"); U t=-n%n,z;do{z=next();}while(z<t);return z%n;}
};
U add_checked(U a,U b){if(a>std::numeric_limits<U>::max()-b)throw std::runtime_error("integer overflow");return a+b;}
U mul_checked(U a,U b){if(b && a>std::numeric_limits<U>::max()/b)throw std::runtime_error("integer overflow");return a*b;}
struct Node {int a,b; U chains,birth_tau,birth_step,write_tau=0,write_step=0; std::vector<U> anc;};
struct State {
 std::vector<Node> nodes; std::unordered_map<U,uint32_t> pairs;
 U F=0,P=0,tau=0,t=0,work=0,forced=0,relieved=0; int H;
 explicit State(int h):H(h){nodes.push_back({-1,-1,1,0,0,0,0,{1}});nodes.push_back({-1,-1,1,0,0,0,0,{2}});}
 U weight(size_t i)const{return mul_checked(2,nodes[i].chains)-2;}
 void step(const std::vector<U>& drawn,std::vector<uint32_t>& W,U &C,U &v){
  auto n=nodes.size(); U q=std::min<U>(3,F+n); if(drawn.size()!=q)throw std::runtime_error("bad q");
  W.clear(); for(U z:drawn){if(z>=F+n)throw std::runtime_error("bad draw"); if(z<n)W.push_back(z);}
  auto sd=drawn;std::sort(sd.begin(),sd.end());if(std::adjacent_find(sd.begin(),sd.end())!=sd.end())throw std::runtime_error("duplicate draw");
  std::sort(W.begin(),W.end());U f=q-W.size(),fm=F-f,pm=add_checked(P,2*f);
  U quota=2*((std::max<U>(1,pm/6)+1)/2); v=(fm>=3 && pm>=6)?std::min({quota,U(H),fm,pm}):0;
  std::vector<std::array<uint32_t,2>> batch;
  for(size_t i=0;i<W.size();++i)for(size_t j=i+1;j<W.size();++j){U key=(U(W[i])<<32)|W[j];if(!pairs.count(key))batch.push_back({W[i],W[j]});}
  // Sequential first/repeat charging equals complete-burst formula; all parents preexist.
  C=0; U newtau=add_checked(tau,W.size()),newt=t+1;
  for(auto uv:batch){uint32_t a=uv[0],b=uv[1],z=nodes.size();auto &va=nodes[a].anc;auto &vb=nodes[b].anc;
   size_t nw=z/64+1;std::vector<U> an(nw,0);
   for(size_t k=0;k<va.size();++k){an[k]|=va[k];} for(size_t k=0;k<vb.size();++k){an[k]|=vb[k];}
   U cost=0;
   for(size_t k=0;k<nw;++k){U mask=an[k];while(mask){int bit=__builtin_ctzll(mask);size_t w=k*64+bit;mask&=mask-1;
    cost=add_checked(cost,mul_checked(nodes[w].write_step?2:11,weight(w)));}}
   C=add_checked(C,cost);
   for(uint32_t p:uv)if(!nodes[p].write_step){nodes[p].write_step=newt;nodes[p].write_tau=newtau;}
   U chains=add_checked(nodes[a].chains,nodes[b].chains);an[z/64]|=U(1)<<(z%64);
   nodes.push_back({int(a),int(b),chains,newtau,newt,0,0,std::move(an)});pairs[(U(a)<<32)|b]=z;
  }
  F=add_checked(fm-v,C);P=pm-v;t=newt;tau=newtau;
  work=add_checked(work,C);forced=add_checked(forced,f);relieved=add_checked(relieved,v);
  if(F!=work-forced-relieved || P!=2*forced-relieved || tau+forced!=3*t-1)throw std::runtime_error("budget");
 }
};
// Event records: six uint64 (step,tau,n_before,F_before,P_before,cost),
// four uint32 (served count,served0/1/2), for 64 bytes. Little-endian hosts only.
int main_standard(int argc,char**argv){try{
 if(argc<6)throw std::runtime_error("usage: native35 prefix H seed target_tau max_steps [audit_steps]");
 std::string prefix=argv[1];int H=std::stoi(argv[2]);U seed=std::stoull(argv[3]),target=std::stoull(argv[4]),maxsteps=std::stoull(argv[5]);
 U audit=argc>6?std::stoull(argv[6]):0; if(H!=0&&H!=6)throw std::runtime_error("H must 0 or6");
 unsigned endian=1;if(*(char*)&endian!=1)throw std::runtime_error("little endian required");
 std::ofstream events(prefix+".events.bin",std::ios::binary),out(prefix+".nodes.tsv"),au(prefix+".audit.tsv");
 if(!events||!out||!au)throw std::runtime_error("output paths");
 State st(H);RNG rng(seed);std::vector<uint32_t> W;std::vector<U> draw;U lastTau=0;
 auto start=std::chrono::steady_clock::now();
 while(st.tau<target && st.t<maxsteps && st.nodes.size()<50000){
  U n=st.nodes.size(),F=st.F,P=st.P,R=n+F,q=std::min<U>(3,R);draw.clear();
  for(U j=0;j<q;++j){U x;do{x=rng.below(R);}while(std::find(draw.begin(),draw.end(),x)!=draw.end());draw.push_back(x);}
  U C,v;st.step(draw,W,C,v);
  if(!W.empty()){
   U fields[6]={st.t,st.tau,n,F,P,C};uint32_t w[4]={uint32_t(W.size()),UINT32_MAX,UINT32_MAX,UINT32_MAX};
   for(size_t k=0;k<W.size();++k){w[k+1]=W[k];} events.write((char*)fields,sizeof(fields));events.write((char*)w,sizeof(w));
  }
  if(st.t<=audit){au<<st.t<<'\t'<<n<<'\t'<<F<<'\t'<<P;for(U x:draw)au<<'\t'<<x;for(U j=q;j<3;++j)au<<"\t-1";
   au<<'\t'<<st.nodes.size()<<'\t'<<st.F<<'\t'<<st.P<<'\t'<<st.tau<<'\t'<<C<<'\t'<<v<<'\n';}
  if(st.tau/65536!=lastTau/65536){std::cerr<<prefix<<" tau="<<st.tau<<" n="<<st.nodes.size()<<" steps="<<st.t<<'\n';lastTau=st.tau;}
 }
 for(size_t i=0;i<st.nodes.size();++i){auto &a=st.nodes[i];out<<i<<'\t'<<a.a<<'\t'<<a.b<<'\t'<<a.chains<<'\t'<<a.birth_tau<<'\t'<<a.birth_step<<'\t'<<a.write_tau<<'\t'<<a.write_step<<'\n';}
 std::ofstream meta(prefix+".meta.json");double secs=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 meta<<"{\"seed\":"<<seed<<",\"H\":"<<H<<",\"target_tau\":"<<target<<",\"tau\":"<<st.tau<<",\"steps\":"<<st.t<<",\"n\":"<<st.nodes.size()<<",\"F\":"<<st.F<<",\"P\":"<<st.P<<",\"total_work\":"<<st.work<<",\"forced_served\":"<<st.forced<<",\"relief_removed\":"<<st.relieved<<",\"rng_state\":"<<rng.x<<",\"complete\":"<<(st.tau>=target?"true":"false")<<",\"runtime_seconds\":"<<secs<<"}\n";
 std::cout<<"tau="<<st.tau<<" t="<<st.t<<" n="<<st.nodes.size()<<" seconds="<<secs<<'\n';return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}

int fork_main(int argc,char**argv){
 if(argc!=7)throw std::runtime_error("--fork nodes_header output seed count max_steps");
 std::ifstream input(argv[2]);std::string output=argv[3];U seed=std::stoull(argv[4]),reps=std::stoull(argv[5]),maxsteps=std::stoull(argv[6]);
 int H;U F,P,t,tau,n;input>>H>>F>>P>>t>>tau>>n;
 State base(H);base.nodes.clear();base.pairs.clear();base.F=F;base.P=P;base.t=t;base.tau=tau;
 base.forced=3*t-1-tau;base.relieved=2*base.forced-P;base.work=F+base.forced+base.relieved;
 for(U i=0;i<n;++i){Node a;input>>a.a>>a.b>>a.chains>>a.birth_tau>>a.birth_step>>a.write_tau>>a.write_step;
  a.anc.assign(i/64+1,0);if(i>=2){for(int p:{a.a,a.b}){for(size_t k=0;k<base.nodes[p].anc.size();++k)a.anc[k]|=base.nodes[p].anc[k];}base.pairs[(U(a.a)<<32)|a.b]=i;}
  a.anc[i/64]|=U(1)<<(i%64);base.nodes.push_back(std::move(a));
 }
 int ka,kb;input>>ka;std::vector<int>A(ka);for(auto&x:A)input>>x;input>>kb;std::vector<int>B(kb);for(auto&x:B)input>>x;if(!input)throw std::runtime_error("fork input malformed");
 std::ofstream out(output);if(!out)throw std::runtime_error("fork output");
 for(U k=0;k<reps;++k){State st=base;RNG rng(seed+k);U na=0,nb=0,un=0;bool gotg=false,gotl=false;
  auto save=[&](const std::string&mode,bool done){out<<k<<'\t'<<seed+k<<'\t'<<mode<<'\t'<<int(done)<<'\t'<<st.t-t<<'\t'<<st.tau-tau<<'\t'<<na<<'\t'<<nb<<'\t'<<un<<'\t'<<st.nodes.size()<<'\t'<<st.F<<'\n';};
  while((!gotg||!gotl)&&st.t-t<maxsteps&&st.tau-tau<262144&&st.nodes.size()<50000){
   U R=st.F+st.nodes.size(),q=std::min<U>(3,R);std::vector<U>draw;
   for(U j=0;j<q;++j){U x;do{x=rng.below(R);}while(std::find(draw.begin(),draw.end(),x)!=draw.end());draw.push_back(x);}
   U C,v;std::vector<uint32_t>W;st.step(draw,W,C,v);
   for(int z:W){bool a=std::find(A.begin(),A.end(),z)!=A.end(),b=std::find(B.begin(),B.end(),z)!=B.end();na+=a;nb+=b;un+=a||b;}
   if(!gotg&&st.tau-tau>=8192){save("global8192",true);gotg=true;}
   if(!gotl&&un>=32){save("local32",true);gotl=true;}
  }
  if(!gotg){save("global8192",false);} if(!gotl){save("local32",false);}
 }
 return 0;
}
int main(int argc,char**argv){try{if(argc>1&&std::string(argv[1])=="--fork")return fork_main(argc,argv);return main_standard(argc,argv);}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}

'''
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

    def confirm_fixed_support(root):
        out=root/'replication'
        out.mkdir(exist_ok=True)
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
        with ThreadPoolExecutor(max_workers=workers) as ex:results=list(ex.map(run,tasks))
        answer=dict(protocol=protocol,tasks=results,wall_seconds=time.perf_counter()-start)
        (out/'replication_results.json').write_text(json.dumps(jsonable(answer),indent=2))
        return answer

    root = Path(workdir) if workdir is not None else Path.cwd()/"DCU_Cell35_run"
    root.mkdir(parents=True,exist_ok=True)
    outdir=root/"results"
    compiler=next((shutil.which(x) for x in ("c++","clang++","g++") if shutil.which(x)),None)
    if compiler is None:
        raise RuntimeError("A C++17 compiler is required for the million-tick generation. "
                           "The bundled replay_readout.py needs only NumPy and does not regenerate histories.")
    source=root/"native35.cpp";binary=root/"native35"
    expected=hashlib.sha256(NATIVE_CPP.encode()).hexdigest()
    pin=root/"NATIVE_SOURCE_SHA256.txt"
    if pin.exists() and pin.read_text().strip()!=expected:
        raise RuntimeError("This run directory contains a different implementation. Use a new directory.")
    source.write_text(NATIVE_CPP,encoding="utf-8")
    if not binary.exists():
        subprocess.run([compiler,"-O3","-std=c++17",str(source),"-o",str(binary)],check=True)
    pin.write_text(expected+"\n")
    if type(workers) is not int or not 1<=workers<=8:raise ValueError("workers must be 1..8")

    def generate(task):
        H,seed=task;prefix=root/f"H{H}_s{seed}"
        meta=prefix.with_suffix(".meta.json")
        if meta.exists() and prefix.with_suffix(".events.bin").exists():
            r=json.loads(meta.read_text())
            if (r['seed'],r['H'],r['target_tau'])!=(seed,H,4456448):raise ValueError("Cached protocol changed")
            return r
        with prefix.with_suffix(".runlog.txt").open("w") as log:
            subprocess.run([str(binary),str(prefix),str(H),str(seed),"4456448","4000000000","16384"],
                           check=True,stdout=log,stderr=log)
        return json.loads(meta.read_text())

    print("CELL 35: actual maintenance-age targets",EPOCHS)
    print("Generating eight native histories; this can take several minutes. No tick-rate extrapolation.")
    with ThreadPoolExecutor(max_workers=workers) as ex:
        metadata=list(ex.map(generate,[(H,s) for H in (0,6) for s in range(20350920,20350924)]))
    result=analyse_campaign(root,outdir)
    print("Reached ages and finite local readouts are saved in",outdir)
    print("PRIMARY Q=27: L3 supports sharing two registers, fresh selection at each epoch")
    for H in (0,6):
        for target in EPOCHS:
            e=result['summary'][str(H)][str(target)]
            if not e['reached_worlds']:continue
            g=e['cohorts']['short_L3']['strata'].get('L3_delta2')
            if g is None:continue
            a=g['modes']['global_windows']['8192'];b=g['modes']['local_windows']['32']
            print(H,target,"P at 8192 global ticks:",a['correct_word']['27']['mean'],
                  "P at 32 local services:",b['correct_word']['27']['mean'],
                  "local censoring:",b['censored_pairs'])
    if include_replication:
        result['fixed_support_replication']=confirm_fixed_support(root)
        print("Fixed-support replication complete: 32 fresh continuations at each of 16 source checkpoints.")
    print("Neither a native particle identity nor a calibrated matter epoch follows from this test.")
    print("Primary native iterations:",result['total_native_iterations'])
    return result


dcu_mass_35 = _run_dcu_mass_35()
