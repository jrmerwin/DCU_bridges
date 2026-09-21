# CELL 39 — infer native service activity from finite quantum outcome counts.
# Self-contained Python + NumPy. No target mass, SI calibration, or space model.
# The estimator receives only the fixed apparatus and its synthetic measurement counts.
# Ideal independent re-preparations/readouts are not a compiled autonomous DCU device.
import json as _j39
from pathlib import Path as _Path39
from collections import Counter
class DCUStructure39:
    """Unordered parent pairs, exact ancestry, and the paper's recording cost."""

    def __init__(self):
        self.parents = [None, None]
        self.chains = [1, 1]
        self.ancestors = [frozenset({0}), frozenset({1})]
        self.recorded = set()
        self.pair_to_id = {}

    def __len__(self):
        return len(self.parents)

    def path_weight(self, node):
        return 2 * self.chains[node] - 2

    def _checked_pairs(self, pairs):
        result = []
        for u, v in pairs:
            if type(u) is not int or type(v) is not int:
                raise TypeError('Parents must be integer IDs in this structure.')
            if not (0 <= u < len(self) and 0 <= v < len(self)):
                raise ValueError('Both parents must already exist before the burst.')
            if u == v:
                raise ValueError('A parent pair must contain two distinct objects.')
            pair = (min(u, v), max(u, v))
            if pair in self.pair_to_id:
                raise ValueError('This pair already exists; it cannot create a duplicate.')
            result.append(pair)
        if len(result) != len(set(result)):
            raise ValueError('A burst cannot contain the same pair twice.')
        return sorted(result)

    def _cost(self, pairs):
        hits = Counter()
        for u, v in pairs:
            hits.update(self.ancestors[u] | self.ancestors[v])
        repeat_part = sum((2 * self.path_weight(z) * count for z, count in hits.items()))
        first_use_extra = sum((9 * self.path_weight(z) for z in hits if z not in self.recorded))
        return repeat_part + first_use_extra

    def recording_cost(self, pairs):
        """Quote the cost without changing the structure."""
        return self._cost(self._checked_pairs(pairs))

    def add_batch(self, pairs):
        """Commit a specified batch. All parents belong to the OLD structure."""
        pairs = self._checked_pairs(pairs)
        cost = self._cost(pairs)
        newborns = []
        for u, v in pairs:
            node = len(self)
            self.parents.append((u, v))
            self.chains.append(self.chains[u] + self.chains[v])
            self.ancestors.append(self.ancestors[u] | self.ancestors[v] | frozenset({node}))
            self.pair_to_id[u, v] = node
            self.recorded.update((u, v))
            newborns.append(node)
        return (newborns, cost)

def _run_dcu_mass_39():
    import numpy as np
    from random import Random
    from itertools import combinations
    from collections import Counter
    from fractions import Fraction
    from math import comb
    checkpoints=(16,64,256)
    shots_per_setting=32768
    repetitions=64
    depths=(27,81)
    settings=np.arange(9,dtype=float)/9
    native_steps=0
    states=[]; raw_native=[]

    def prepare():
        st=DCUStructure39()
        b,c=st.add_batch([(0,1)]);assert b==[2] and c==0
        b,c=st.add_batch([(0,2),(1,2)]);assert b==[3,4] and c==26
        return st,26,0,5

    def native_step(st,F,P,rng,H):
        n=len(st);R=n+F;q=min(3,R)
        selected=tuple(sorted(j for j in rng.sample(range(R),q) if j<n))
        f=q-len(selected);Fm=F-f;Pm=P+2*f
        quota=2*((max(1,Pm//6)+1)//2)
        relief=min(quota,H,Fm,Pm) if Fm>=3 and Pm>=6 else 0
        pairs=[e for e in combinations(selected,2) if e not in st.pair_to_id]
        # Independent sequential first/repeat computation before original batch commit.
        hits=Counter(z for a,b in pairs for z in st.ancestors[a]|st.ancestors[b])
        new_regs=sorted(z for z in hits if z>=2 and z not in st.recorded)
        repeat_work=sum(2*st.path_weight(z)*v for z,v in hits.items())
        first_premium=sum(9*st.path_weight(z) for z in new_regs)
        rec=set(st.recorded);expected=0
        for a,b in sorted(pairs):
            for z in st.ancestors[a]|st.ancestors[b]:
                expected+=(2 if z in rec else 11)*st.path_weight(z)
            rec.update(st.ancestors[a]|st.ancestors[b])
        born,cost=st.add_batch(pairs);assert cost==expected==repeat_work+first_premium
        assert st.recorded=={v for pp in st.parents[2:] for v in pp}
        return Fm-relief+cost,Pm-relief,len(selected),dict(n=n,F=F,P=P,served=list(selected),
                      forced_served=f,relief_removed=relief,pairs=pairs,born=born,cost=cost,
                      repeat_work=repeat_work,first_use_premium=first_premium,
                      newly_recorded_nonprimitive=new_regs,
                      new_record_information_trits=len(new_regs))

    st,F,P,tau=prepare()
    states.append(dict(label='preparation_control',H=None,step=0,n=len(st),F=F,P=P,R=len(st)+F,tau=tau,parents=st.parents.copy()))
    for H in (0,6):
        st,F,P,tau=prepare();rng=Random(20380301)
        for t in range(1,max(checkpoints)+1):
            F,P,s,event=native_step(st,F,P,rng,H);tau+=s;native_steps+=1
            raw_native.append(dict(H=H,step=t,**event))
            if t in checkpoints:
                states.append(dict(label=f'H{H}_t{t}',H=H,step=t,n=len(st),F=F,P=P,R=len(st)+F,tau=tau,
                                   parents=st.parents.copy()))
    def exact_law(R):
        p=Fraction(3,R);p2=Fraction(6,R*(R-1))
        return (1-2*p+p2,2*(p-p2),p2)

    def apparatus(Q):
        x=settings[:,None,None]+np.arange(3)[None,None,:]/Q-np.arange(3)[None,:,None]/3
        a=((1+2*np.cos(2*np.pi*x))**2)/9
        a[a<1e-28]=0
        assert np.max(abs(a.sum(axis=1)-1))<2e-14
        return a.reshape(27,3)

    def weights(p):
        p=np.asarray(p);p2=2*p*p/(3-p)
        return np.stack([1-2*p+p2,2*(p-p2),p2],axis=-1)

    def make_estimator(A):
        # No actual n,F,R,K, true rate, state name, or data-generation RNG enters.
        grid=np.r_[0.,np.geomspace(1e-8,0.6,2048)]
        kernel=(A@weights(grid).T).T
        logs=np.log(np.maximum(kernel,1e-300))
        def estimate(counts):
            counts=np.asarray(counts,dtype=np.int64).reshape(-1)
            assert counts.shape==(27,) and np.all(counts>=0)
            coarse=logs@counts
            i=int(np.argmax(coarse))
            left=grid[max(0,i-1)];right=grid[min(len(grid)-1,i+1)]
            def loglike(p):
                probs=A@weights(p)
                if np.any((probs<=0)&(counts>0)):return -np.inf
                good=counts>0
                return float(np.dot(counts[good],np.log(probs[good])))
            # Bounded golden section inside the globally scanned maximum's bracket.
            ratio=(np.sqrt(5)-1)/2
            a,b=left,right;c=b-ratio*(b-a);d=a+ratio*(b-a)
            fc,fd=loglike(c),loglike(d)
            for _ in range(72):
                if fc>fd:
                    b,d,fd=d,c,fc;c=b-ratio*(b-a);fc=loglike(c)
                else:
                    a,c,fc=c,d,fd;d=a+ratio*(b-a);fd=loglike(d)
            candidates=[(a+b)/2,grid[i],left,right,0.,0.6]
            phat=max(candidates,key=loglike)
            return float(phat),loglike(phat)
        return estimate

    outputs={'protocol':'Cell39 inverse activity v1','native_law':{'Gamma':3,'m':0,'H':[0,6]},
             'states':states,'native_steps_generated':native_steps,'native_events':raw_native,
             'shots_per_setting':shots_per_setting,'settings':[f'{i}/9' for i in range(9)],
             'repetitions_per_state_depth':repetitions,'primary_phase_modulus':27,
             'measurement_domain':'independent one-step re-preparations at each fixed saved state',
             'new_physical_measurements':False,'F_n_separately_identifiable':False,
             'unknown_parameter':'p=3/(F+n)','depths':{}}
    checks=dict(maintenance_subsets=0,weighted_token_combinations=0,direct_tensor_probabilities=0,
                estimator_runs=0,native_preparation_steps=native_steps)
    print('CELL 39: INVERSE NATIVE ACTIVITY FROM QUANTUM OUTCOMES')
    print('Nine frozen offsets; 32,768 ideal shots per offset; 64 synthetic datasets per state/depth.')
    print('Estimator gets counts and fixed apparatus only; no hidden native rates or service histories.')
    for state in states:
        n,F,R=state['n'],state['F'],state['R'];tally=Counter()
        for size in range(4):
            if size>n or 3-size>F:continue
            factor=comb(F,3-size)
            for ss in combinations(range(n),size):
                tally[int(3 in ss)+int(4 in ss)]+=factor
                checks['maintenance_subsets']+=1
        assert sum(tally.values())==comb(R,3)
        assert tuple(Fraction(tally[k],comb(R,3)) for k in range(3))==exact_law(R)
        checks['weighted_token_combinations']+=comb(R,3)
    for Q in depths:
        A=apparatus(Q)
        singular=np.linalg.svd(A,compute_uv=False)
        assert len(singular)==3 and singular[-1]>1e-8
        # Direct two-factor amplitudes, independent of the cosine-kernel implementation.
        for setting in settings:
            for K in range(3):
                direct=np.zeros(3)
                for a in range(3):
                    for b in range(3):
                        amp=sum(np.exp(2j*np.pi*r*(setting+K/Q+(a-b)/3)) for r in range(3))/np.sqrt(27)
                        direct[(b-a)%3]+=abs(amp)**2;checks['direct_tensor_probabilities']+=1
                ref=A.reshape(9,3,3)[int(round(setting*9)),:,K]
                assert np.max(abs(direct-ref))<2e-14
        estimator=make_estimator(A)
        qout={'matrix_singular_values':singular.tolist(),'matrix_condition_number':float(singular[0]/singular[-1]),
              'rows':[]}
        print(f'  Q={Q}: mixture matrix rank 3, condition number {qout["matrix_condition_number"]:.3f}')
        for state_index,state in enumerate(states):
            truth=float(Fraction(3,state['R']))
            probs=(A@np.array(list(map(float,exact_law(state['R']))))).reshape(9,3)
            probs=np.clip(probs,0,1);probs/=probs.sum(axis=1,keepdims=True)
            rr=np.random.Generator(np.random.PCG64(20380303+1000*Q+state_index))
            estimates=[];counts_all=[];ll=[]
            for rep in range(repetitions):
                counts=np.array([rr.multinomial(shots_per_setting,p) for p in probs])
                p_est,l=estimator(counts)
                estimates.append(p_est);counts_all.append(counts.tolist());ll.append(l);checks['estimator_runs']+=1
            e=np.array(estimates);quant=np.quantile(e,[.025,.5,.975])
            rec={'label':state['label'],'true_R':state['R'],'true_p':truth,
                 'median_estimated_p':float(quant[1]),'p_empirical_95percent_range':[float(quant[0]),float(quant[2])],
                 'mean_bias_p':float(np.mean(e)-truth),'RMSE_p':float(np.sqrt(np.mean((e-truth)**2))),
                 'relative_RMSE_p':float(np.sqrt(np.mean((e-truth)**2))/truth),
                 'median_implied_R':None if quant[1]==0 else float(3/quant[1]),
                 'boundary_zero_count':int(sum(e==0)),'estimates_p':estimates,'counts':counts_all,'maximized_loglikelihood':ll}
            qout['rows'].append(rec)
            print(f'    {state["label"]:20s} R={state["R"]:5d}  p={truth:.6f} '
                  f'median p_hat={quant[1]:.6f}  relative RMSE={100*rec["relative_RMSE_p"]:.2f}%  zeros={rec["boundary_zero_count"]}')
        outputs['depths'][str(Q)]=qout
    # Same total pool, different split -> exactly the same marginal instrument law.
    outputs['identifiability']={'mixture_weights_identifiable_from_exact_distributions':True,
        'proof':'Fourier harmonics reveal C0,C1,C2; their three K=0,1,2 phase columns form a nonsingular Vandermonde matrix.',
        'F_n_separation':'Only R=F+n enters. Fixing n additionally would be extra information, not read from these outcomes.',
        'noiseless_versus_finite_shot':'Exact identifiability does not imply useful precision at finite counts; no phase/budget refit.'}
    # Secondary explanatory audit added after the first information/metrology run.
    # All events retained; a trit is written once, not redrawn on every service.
    productive=[e for e in raw_native if e['born']]
    work_groups={}
    for v in sorted(set(e['new_record_information_trits'] for e in productive)):
        cases=[e for e in productive if e['new_record_information_trits']==v]
        work_groups[str(v)]={'events':len(cases),'min_charge':min(e['cost'] for e in cases),
                             'max_charge':max(e['cost'] for e in cases),
                             'sum_charge':sum(e['cost'] for e in cases)}
    outputs['work_information_audit']={'status':'secondary exact bookkeeping on all new native preparation events',
        'positive_charge_no_new_trit_events':sum(e['cost']>0 and not e['newly_recorded_nonprimitive'] for e in productive),
        'productive_events':len(productive),'groups_by_new_information_trits':work_groups,
        'new_record_entropy_bits':sum(e['new_record_information_trits'] for e in productive)*float(np.log2(3)),
        'interpretation':'Conditional ensemble entropy of newly written uniform trits, accessible if their values are read. Not thermodynamic entropy production or work-to-joule conversion.'}
    print('Native recording charge versus new record entropy:',work_groups)
    print('Positive work without any new trit:',outputs['work_information_audit']['positive_charge_no_new_trit_events'])
    outputs['validation']=checks
    print('PASS:',checks)
    print('Sampling ranges describe repeated synthetic datasets, not laboratory intervals. Smaller p is harder to estimate.')
    print('No SI clock, mass, registry-specific service weight, local redshift, or free autonomous apparatus supplied.')
    return outputs

dcu_mass_39 = _run_dcu_mass_39()
_out39=_Path39.cwd()/'DCU_Cell_39_outputs';_out39.mkdir(exist_ok=True)
(_out39/'DCU_Mass_Cell_39_RESULTS.json').write_text(_j39.dumps(dcu_mass_39,indent=2)+'\n',encoding='utf-8')
