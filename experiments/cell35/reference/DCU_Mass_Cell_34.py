# CELL 34 — exact protected-state audit of the existing record-path architecture.
# Paste after Cell26 (or any later cell). NumPy + standard library only.
# PRIMARY: unchanged Cell26 +/- support-count coupling to its TWO-QUTRIT apparatus.
# Paths are classical control inputs, NOT newly asserted quantum tensor factors.
# No echoes, fitted noise, native-law changes, particle assignments or SI calibration.
from collections import Counter, defaultdict
from copy import deepcopy
from fractions import Fraction
from itertools import combinations, product
from math import comb, sqrt


def _run_dcu_mass_34(saved=None):
    import numpy as np
    if saved is None:
        if 'dcu_mass_26' not in globals():
            raise RuntimeError('Run Cell26 first, or use the standalone Cell34 reproduce.py.')
        saved = dcu_mass_26
    before = deepcopy(saved)
    def check(ok, message):
        if not ok:
            raise AssertionError(message)
    def get(d, key):
        return d[key] if key in d else d[str(key)]
    def stats(v):
        a=np.asarray(v,dtype=float)
        return dict(n=len(a),mean=float(a.mean()),
                    se=float(a.std(ddof=1)/sqrt(len(a))) if len(a)>1 else None)
    def rat(q):
        return str(q.numerator) if q.denominator==1 else f'{q.numerator}/{q.denominator}'
    check(tuple(saved['record_lengths'])==(3,4), 'Original record lengths changed')
    check(tuple(saved['phase_depths'])==(3,4), 'Original phase moduli changed')
    check(saved['freeze_step']==512 and saved['horizon']==2048, 'Original windows changed')
    times=tuple(saved['checkpoints']); seeds=tuple(saved['seeds'])
    check(len(seeds)==256, 'Do not select a subset of the archived histories')
    check(not saved['entry_reflection_used_as_phase_sign'], 'Do not replace entry reflection by a phase sign')
    basis=tuple(product(range(3),repeat=2))
    refl=(1,0,2)
    counts=Counter(); errors=defaultdict(float)

    def signature_blocks(path_a,path_b,Q,entry_sensitive=False):
        # One integer signature component per ACTUAL register; exact modulo Q.
        # G=diag(0,1,2). Secondary model uses S^e G S^e, not (-1)^e G.
        ea=dict(path_a[1]); eb=dict(path_b[1])
        regs=tuple(sorted(set(ea)|set(eb)))
        blocks=defaultdict(list)
        for i,(r,s) in enumerate(basis):
            sig=[]
            for z in regs:
                left=(refl[r] if entry_sensitive and ea[z] else r) if z in ea else 0
                right=(refl[s] if entry_sensitive and eb[z] else s) if z in eb else 0
                sig.append((left-right)%Q)
            blocks[tuple(sig)].append(i)
        return tuple(sorted(tuple(v) for v in blocks.values()))

    def enumerate_paths(parentage,length):
        recorded={v for pair in parentage if pair is not None for v in pair}
        cache={0:[(0,())],1:[(1,())]}
        for z,pair in enumerate(parentage[2:],2):
            cache[z]=[(root,tokens+((z,e),)) for e,p in enumerate(pair)
                      for root,tokens in cache[p] if len(tokens)<length]
        return tuple(sorted((root,tokens) for z in sorted(recorded) if z>=2
                            for root,tokens in cache[z] if len(tokens)==length))

    rr=np.arange(3)
    omega=np.zeros(9,complex);omega[[0,4,8]]=1/sqrt(3)
    # Three source phase words; D(k/3) on Alice BEFORE storage, no in-storage controls.
    codewords=np.array([omega*np.kron(np.exp(2j*np.pi*rr*k/3),np.ones(3)) for k in range(3)])
    check(np.max(abs(codewords@codewords.conj().T-np.eye(3)))<2e-14,'Orthogonal phase words')
    projectors=np.empty((3,3,9,9),complex)
    for a,b in product(range(3),repeat=2):
        av=np.exp(-2j*np.pi*rr*a/3)/sqrt(3)
        bv=np.exp(+2j*np.pi*rr*b/3)/sqrt(3)
        v=np.kron(av,bv);projectors[a,b]=np.outer(v,v.conj())
    check(np.max(abs(projectors.sum(axis=(0,1))-np.eye(9)))<2e-14,'Terminal POVM completeness')
    def decode(rho):
        p=np.real(np.einsum('ij,abji->ab',rho,projectors))
        check(p.min()>-3e-14 and abs(p.sum()-1)<3e-14,'Decoded probability validity')
        return np.array([sum(p[a,b] for a,b in product(range(3),repeat=2) if (b-a)%3==k)
                         for k in range(3)])
    ideal_confusion=np.array([decode(np.outer(v,v.conj())) for v in codewords])
    erased=sum(np.diag(abs(v)**2) for v in codewords)/3
    check(np.max(abs(ideal_confusion-np.eye(3)))<2e-14,'Codewords are readable')
    check(np.max(abs(decode(erased)-1/3))<2e-14,'Readout must actually use coherence')
    # Check a spanning matrix-unit basis: protection is not just classical populations.
    logical_matrix_units=tuple((r,s) for r,s in product(range(3),repeat=2))

    out={}
    for name in ('no_relief','relief_H6'):
        env=saved['environments'][name]
        parentage=tuple(None if p is None else tuple(p) for p in env['parentage'])
        n=len(parentage);F=int(env['checkpoint']['F']);R=F+n;q=3
        check(F>=3 and n>=5, 'Identity/single-register service witnesses require this pool')
        zero_prob=Fraction(comb(F,3),comb(R,3))
        singleton_prob=Fraction(comb(F,2),comb(R,3))
        check(zero_prob>0 and singleton_prob>0,'Every register generator must be a legal native event')
        hist=env['histories']
        check(tuple(h['seed'] for h in hist)==seeds,'History selection or ordering changed')
        service=np.asarray([h['counts'] for h in hist],dtype=np.int64)
        check(service.shape==(256,len(times),n),'Unexpected per-object service data')
        check(np.all(service>=0) and np.all(np.diff(service,axis=1)>=0),'Invalid cumulative counts')
        eout=dict(checkpoint=deepcopy(env['checkpoint']),identity_service_probability=rat(zero_prob),
                  singleton_register_service_probability=rat(singleton_prob),panels={})
        for length in (3,4):
            panel=get(env['panels'],length)
            paths=tuple((int(root),tuple(tuple(map(int,tok)) for tok in tokens)) for root,tokens in panel['paths'])
            check(paths==enumerate_paths(parentage,length),'Archived paths do not equal complete native panel')
            supports=tuple(tuple(z for z,e in tokens) for root,tokens in paths)
            check(supports==tuple(tuple(s) for s in panel['supports']),'Support mapping changed')
            counts['record_paths']+=len(paths)
            pairlist=tuple(combinations(range(len(paths)),2))
            rows=[];grouped=defaultdict(list)
            for i,j in pairlist:
                same=supports[i]==supports[j]
                if same:
                    check(paths[i][0]!=paths[j][0] and paths[i][1][0][0]==2
                          and paths[j][1][0][0]==2 and paths[i][1][1:]==paths[j][1][1:],
                          'Equal support must not hide a different composite chain')
                    counts['genesis_entry_alias_pairs']+=1
                delta=len(set(supports[i])^set(supports[j]))
                blocks=signature_blocks(paths[i],paths[j],27)
                secondary=signature_blocks(paths[i],paths[j],27,True)
                dims=sorted(map(len,blocks),reverse=True)
                check(dims==([3,2,2,1,1] if same else [1]*9),'Complete pair-level classification')
                for Q in (27,81):
                    check(signature_blocks(paths[i],paths[j],Q)==blocks,'Unexpected refinement alias')
                    check(signature_blocks(paths[i],paths[j],Q,True)==secondary,'Entry model refinement changed')
                    counts['exact_noise_classifications']+=2
                check(max(map(len,secondary))==1,'Secondary full-entry model unexpectedly degenerate')
                grouped[delta].append(len(rows))
                rows.append(dict(pair=(i,j),shared_registers=length-delta//2,symmetric_difference=delta,
                    same_support=same,protected_block_dimensions=dims,blocks=blocks,
                    secondary_entry_sensitive_max_dimension=max(map(len,secondary)),terminal={}))
                counts['path_pairs']+=1;counts['protected_pairs']+=int(same)
            incidence=np.zeros((len(paths),n),dtype=np.int64)
            for i,ss in enumerate(supports):incidence[i,list(ss)]=1
            pc=service@incidence.T
            ii,jj=np.asarray(pairlist).T
            difference=pc[:,:,ii]-pc[:,:,jj]
            positive=grouped.get(0,[])
            if positive:
                check(np.all(difference[:,:,positive]==0),'Shared support clocks are not identical pathwise')
                counts['exact_zero_difference_snapshots']+=len(positive)*256*len(times)
            curves={}
            for depth in (3,4):
                Q=3**depth
                phase=(difference%Q)/Q
                zz=np.exp(2j*np.pi*phase)
                correct=(1+2*np.cos(2*np.pi*phase))**2/9
                for k,row in enumerate(rows):
                    z=zz[:,-1,k];c1=complex(z.mean());c2=complex((z*z).mean())
                    row['terminal'][depth]=dict(correct_phase_word=stats(correct[:,-1,k]),
                        C1=c1,C2=c2,negativity=float((2*abs(c1)+abs(c2))/3))
                gcurve={}
                for delta,indices in sorted(grouped.items()):
                    v=correct[:,:,indices].mean(axis=2)
                    prior=get(get(panel['groups'],delta)['phase_results'],depth)
                    oldmean=np.asarray(prior['mean_P0'])
                    error=float(np.max(abs(v.mean(axis=0)-oldmean)))
                    errors['old_cell26_P0_regression']=max(errors['old_cell26_P0_regression'],error)
                    check(error<1e-12,'Original Cell26 readout changed')
                    counts['prior_group_checkpoint_regressions']+=len(times)
                    gcurve[delta]=dict(pair_count=len(indices),correct_phase_word_by_time=tuple(stats(v[:,k]) for k in range(len(times))))
                curves[depth]=gcurve
                # Direct 9D check at the first pair of EVERY overlap class, endpoints of seed blocks.
                for delta,indices in sorted(grouped.items()):
                    i,j=pairlist[indices[0]]
                    for h in (0,127,128,255):
                        for t in (0,len(times)//2,len(times)-1):
                            a,b=int(pc[h,t,i]),int(pc[h,t,j])
                            u=np.kron(np.exp(2j*np.pi*rr*(a%Q)/Q),np.exp(-2j*np.pi*rr*(b%Q)/Q))
                            for k,v in enumerate(codewords):
                                noise=u*v;rho=np.outer(noise,noise.conj());decoded=decode(rho)
                                theory=np.array([(1+2*np.cos(2*np.pi*((a-b)/Q+(k-word)/3)))**2/9 for word in range(3)])
                                err=float(np.max(abs(decoded-theory)))
                                errors['dense_decode']=max(errors['dense_decode'],err)
                                check(err<2e-13,'Dense source/readout mismatch')
                                counts['dense_encoded_state_checks']+=1
                            if delta==0:
                                for r,s in logical_matrix_units:
                                    err=abs(u[4*r]*u[4*s].conjugate()-1)
                                    errors['protected_matrix_units']=max(errors['protected_matrix_units'],float(err))
                                    check(err<2e-14,'Protected coherence changed')
                                    counts['protected_matrix_unit_checks']+=1
            # Informational entry alias, and the explicitly DIFFERENT coupling's response.
            example=None
            if positive:
                k=positive[0];i,j=pairlist[k];S=supports[i]
                ea=dict(paths[i][1]);eb=dict(paths[j][1]);agree=0
                for values in product(range(3),repeat=length):
                    g=dict(zip(S,values))
                    read=lambda path:tuple(refl[g[z]] if e else g[z] for z,e in path[1])
                    agree+=int(read(paths[i])==read(paths[j]))
                rates={}
                for Q in (27,81):
                    # Diagonal-source phase on |r,r> for all retained histories.
                    phase=np.zeros((256,len(times),3),dtype=np.int64)
                    for z in S:
                        shifts=np.array([(refl[r] if ea[z] else r)-(refl[r] if eb[z] else r) for r in range(3)])
                        phase+=service[:,:,z,None]*shifts
                    fidelity=abs(np.exp(2j*np.pi*(phase%Q)/Q).mean(axis=2))**2
                    rates[Q]=dict(secondary_correct_word_terminal=stats(fidelity[:,-1]),
                                  primary_correct_word_terminal=1.0)
                example=dict(pair=(i,j),paths=(paths[i],paths[j]),support=S,
                    address_agreement_probability=rat(Fraction(agree,3**length)),
                    secondary_noise_is_new_hypothesis=True,secondary_not_a_phase_sign=True,
                    readout_sensitivity=rates)
            eout['panels'][length]=dict(path_count=len(paths),paths=paths,supports=supports,
                distinct_support_count=len(set(supports)),path_pair_count=len(pairlist),
                protected_pair_count=len(positive),unique_support_pair_count=comb(len(set(supports)),2),
                protected_pairs_after_unique_support_quotient=0,
                secondary_entry_sensitive_protected_pairs=0,pair_results=rows,curves=curves,
                representative_alias=example)
        out[name]=eout
    check(saved==before,'The supplied data were modified')
    print('CELL 34 — ARCHITECTURE -> PROTECTED AND READABLE STATE (EXACT PAIR TEST)')
    print('Primary: original two-qutrit apparatus, original +/- path-support coupling; no storage pulses.')
    print('Native input: all four Cell26 panels and all 512 saved histories; no new native simulation.')
    print('\n environment       length paths supports pairs protected(max dim 3) unique-support protected')
    for name,e in out.items():
        for length,p in e['panels'].items():
            print(f' {name:17s} {length:3d} {p["path_count"]:5d} {p["distinct_support_count"]:8d}'
                  f' {p["path_pair_count"]:5d} {p["protected_pair_count"]:20d} {0:24d}')
    print('\nExact classification: identical supports -> blocks [3,2,2,1,1]; distinct equal-length supports -> nine singletons.')
    print('All 36 equal-support pairs differ only in primitive entry into register 2; their composite chains coincide.')
    print('All-time protection follows from integer character equality; initial legal singleton draws prove necessity.')
    print('No approximate zero, finite-ensemble eigenvalue cutoff, or target mass is used.')
    print('Protected diagonal code: span{|00>,|11>,|22>}. Three phase codewords decode perfectly via b-a modulo 3.')
    print('Codeword populations and local marginals are identical. Full phase erasure reduces decoding success to 1/3.')
    print('\nPRIMARY length-3 terminal correct-word probability, grouped by shared registers:')
    for name in out:
        p=out[name]['panels'][3]
        for delta,g in p['curves'][3].items():
            z=g['correct_phase_word_by_time'][-1]
            print(f' {name:17s} shared={3-delta//2}: {z["mean"]:.9f} +/- {z["se"]:.9f} MC SE')
    print('\nENTRY-SENSITIVE SENSITIVITY (a DIFFERENT coupling, not a revision to the primary):')
    for name in out:
        p=out[name]['panels'][3];ex=p['representative_alias']
        z=ex['readout_sensitivity'][27]['secondary_correct_word_terminal']
        print(f' {name:17s} first alias {ex["pair"]}, support={ex["support"]},'
              f' address agreement={ex["address_agreement_probability"]}; secondary success={z["mean"]:.9f}')
    print('Actual S=(0 1) conjugation is NOT sign inversion. Every distinct path pair is nondegenerate in this secondary model.')
    print('\nPASS:',dict(counts))
    print('Max numerical errors:',dict(errors))
    print('Verdict: usable conditional shared-input quantum memory; NO extra two-path architecture-specific protection.')
    print('Not a particle, mass mechanism, autonomous memory, readout-cost derivation, SI map or spatial result.')
    return dict(protocol='Cell34 v1 exact pairwise DFS classification and readable phase-word transfer',
        environments=out,checkpoints=times,seeds=seeds,phase_depths=(3,4),primary_phase_depth=3,
        codeword_definition='Omega_k=sum_r exp(2pi*i*k*r/3)|rr>/sqrt(3), k=0,1,2',
        logical_code_basis=(0,4,8),ideal_decode_confusion=ideal_confusion,
        erased_decode_distribution=decode(erased),exact_counts=dict(counts),max_numerical_errors=dict(errors),
        new_native_updates=0,reused_native_histories=512,
        primary_coupling_changed=False,noise_entry_sensitive_variant_is_added=True,
        arbitrary_record_paths_are_quantum_factors=False,apparatus_is_existing_two_sibling_instrument=True,
        accessible_three_word_encoding_is_ideal_declared_instrument=True,full_logical_control_compiled=False,
        protected_against='fixed support-count diagonal phase operations only',
        protection_beyond_identical_driving_supports=False,storage_feedback=False,
        native_gate_memory_cost_derived=False,neutron_search_modified=False,
        physical_calibration=None,particle_identification=None)


dcu_mass_34 = _run_dcu_mass_34()
