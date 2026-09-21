# CELL 24 — finite record phases, collective readout, and energy-map selection.
# Paste after Cell 23. Requires NumPy and the _Native constructor from Cell 1.
# No network/data downloads. No mass, screening, service, or clock calibration fit.
# Source basis: revised DCU manuscript §§5.1,5.4,7.1 and Cell23's declared phases.
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from itertools import combinations
from math import cos, gcd, pi, sqrt


def _run_dcu_mass_24():
    try:
        import numpy as np
    except ImportError as exc:
        raise ImportError('Cell 24 needs NumPy in the notebook kernel.') from exc
    needed = ('dcu_mass_23', '_reference', '_Native')
    missing = [k for k in needed if k not in globals()]
    if missing:
        raise RuntimeError('Run Cell 23 first. Missing: ' + ', '.join(missing))
    old = dcu_mass_23
    before = deepcopy((old, _reference))
    Q = Fraction
    if (len(_reference['rows']), Counter(r['sector'] for r in _reference['rows'])) != \
            (137, Counter(S=81, I=40, G=16)):
        raise ValueError('Original registry changed; no substitute selected.')
    if old['quantum']['Kn_gap_choices'] != {
            n: dict(signed=1+Q(1,n-1), magnitude=1-Q(1,n-1)) for n in (3,4,5)}:
        raise ValueError('Cell23 energy-map controls changed.')

    def triadic(x):
        """Exact element of Z[1/3]/Z; accepts Fractions/integers, NOT rounded floats."""
        if isinstance(x, (float, np.floating)):
            raise TypeError('Use an exact Fraction for a finite native phase setting.')
        x = Q(x) % 1
        denominator = x.denominator
        while denominator % 3 == 0:
            denominator //= 3
        if denominator != 1:
            raise ValueError('This setting is not in a finite triadic catalog.')
        return x

    def phases(charges, x):
        x = triadic(x)
        return np.array([np.exp(2j*pi*float((int(q)*x) % 1)) for q in charges])

    labels = np.arange(3)
    F = np.exp(2j*pi*np.outer(labels, labels)/3)/sqrt(3)
    Fd = F.conj().T
    initial = np.array([1,0,0], dtype=complex)

    def ramsey(x, charges=(0,1,2)):
        return np.abs(Fd @ (phases(charges, x)*(F @ initial)))**2

    def kernel(x):
        # This is the source's L, evaluated only at exact triadic settings.
        x = triadic(x)
        return (1+2*cos(2*pi*float(x)))**2/9

    # A. Exact congruences underpin finite-clock composition and refinement.
    clock_checks = 0
    group_checks = 0
    for depth in range(1,6):
        q = 3**depth
        for k in range(q):
            x = Q(k,q)
            assert x == Q(3*k,3*q)  # SAME coarse tick, NOT one new fine tick.
            p = phases((0,1,2),x)
            assert np.max(np.abs(p-phases((0,1,2),Q(3*k,3*q)))) < 1e-14
            analytic = np.array([kernel(x-Q(a,3)) for a in range(3)])
            assert np.max(np.abs(ramsey(x)-analytic)) < 3e-14
            assert abs(sum(analytic)-1) < 3e-14
            clock_checks += 1
        if depth <= 4:
            for k in range(q):
                for l in range(q):
                    assert all((r*Q(k,q)+r*Q(l,q)-r*Q((k+l)%q,q)).denominator == 1
                               for r in (0,1,2))
                    group_checks += 1
    assert np.max(np.abs(ramsey(Q(1,9))-old['quantum']['interference']['1/9'])) < 1e-14
    # Both values are genuinely different elements: refinement does NOT slow an old tick.
    assert triadic(Q(1,3)) == triadic(Q(3,9)) != triadic(Q(1,9))

    # One-depth logarithm ambiguity; the full source phase family resolves it.
    # These are GENERATOR LABELS in D(x)=exp(+2pi i x G), not physical energies.
    main_labels, alias_labels = (0,1,2), (0,1,11)
    coarse_error = 0.0
    for k in range(9):
        assert all(Q((a-b)*k,9).denominator == 1 for a,b in zip(main_labels,alias_labels))
        coarse_error = max(coarse_error, float(np.max(np.abs(
            phases(main_labels,Q(k,9))-phases(alias_labels,Q(k,9))))))
    refined_x = Q(1,27)  # Fixed first new fine setting, not chosen for maximum contrast.
    p_main, p_alias = ramsey(refined_x,main_labels), ramsey(refined_x,alias_labels)
    alias = dict(coarse_denominator=9, primary_labels=main_labels, alternative_labels=alias_labels,
        coarse_max_gate_error=coarse_error, refined_setting=refined_x,
        primary_probabilities=tuple(map(float,p_main)), alias_probabilities=tuple(map(float,p_alias)),
        refined_total_variation=float(np.sum(np.abs(p_main-p_alias))/2),
        alternative_obeys_entire_source_phase_family=False)
    assert alias['refined_total_variation'] > .1
    # Fourier harmonics of the existing local readout, NOT a fitted spectrum.
    phase_grid = np.arange(81)/81
    trace = np.array([kernel(Q(k,81)) for k in range(81)])
    coeff = np.fft.fft(trace)/81
    expected = np.zeros(81,dtype=complex)
    expected[0]=1/3
    expected[1]=expected[-1]=2/9
    expected[2]=expected[-2]=1/9
    assert np.max(np.abs(coeff-expected)) < 2e-14

    # B. The declared shared-source state + factorwise copies of the same setting/readout.
    # j is NUMBER OF SIBLINGS, not K_j graph order or particle identity.
    # Tiny native fixtures establish legal pair preparation, not a service-history sample.
    def prepare_siblings(j):
        state = _Native()
        setup = 0
        while len(state) < j+1:
            _, work = state.add_batch([(0,len(state)-1)])
            setup += work
        source = len(state)-1
        assert source not in state.recorded
        partners = tuple(range(j))
        served = tuple(range(len(state)))
        # TG1/CO1: the burst includes EVERY absent pair of the served old objects,
        # not just the desired siblings. Preserve the extra children too.
        events = [e for e in combinations(served,2) if e not in state.pair_to_id]
        quoted = state.recording_cost(events)
        born, work = state.add_batch(events)
        siblings = tuple(v for v in born if source in state.parents[v])
        assert work == quoted and len(siblings) == j
        assert len({p for pair in events for p in pair}) == j+1
        assert all(e in state.pair_to_id for e in combinations(served,2))
        return dict(parentage=tuple(state.parents), common_source=source, siblings=siblings,
            other_children=tuple(v for v in born if v not in siblings), co_served_set=served,
            setup_work=setup, entire_burst_work=work, minimum_co_service_capacity=j+1,
            all_enabled_pairs_constructed=True, service_history_simulated=False,
            quantum_state_is_declared_instrument=True)

    def sibling_readout(settings, coherent=True):
        """Local phases, local inverse Fourier measurements, coarse outcome sum mod3."""
        xs = tuple(triadic(x) for x in settings)
        j = len(xs)
        if not 2 <= j <= 5:
            raise ValueError('This finite audit uses 2..5 sibling factors.')
        shape = (3,)*j
        if coherent:
            state = np.zeros(shape,dtype=complex)
            for r in range(3):
                phase = 1.0+0j
                for x in xs:
                    phase *= phases((0,1,2),x)[r]
                state[(r,)*j] = phase/sqrt(3)
            for axis in range(j):
                state = np.moveaxis(np.tensordot(Fd,state,axes=(1,axis)),0,axis)
            probabilities = np.abs(state)**2
        else:
            # Independently transform each diagonal-mixture component. This erases
            # only inter-label coherence, not the classical equal-label correlations.
            probabilities = np.zeros(shape)
            for r in range(3):
                vec = np.array([1+0j])
                for x in xs:
                    vec = np.kron(vec,Fd[:,r]*phases((0,1,2),x)[r])
                probabilities += np.abs(vec.reshape(shape))**2/3
        total = sum(xs,Q(0)) % 1
        sums = np.indices(shape).sum(axis=0) % 3
        marginal = np.array([probabilities[sums==a].sum() for a in range(3)])
        target = np.array([kernel(total-Q(a,3)) for a in range(3)]) if coherent else np.ones(3)/3
        assert np.max(np.abs(marginal-target)) < 5e-14 and abs(probabilities.sum()-1) < 5e-14
        if coherent:
            joint_target = np.array([target[a] for a in sums.ravel()]).reshape(shape)/(3**(j-1))
            assert np.max(np.abs(probabilities-joint_target)) < 5e-14
        for axis in range(j):
            local = probabilities.sum(axis=tuple(k for k in range(j) if k != axis))
            assert np.max(np.abs(local-1/3)) < 5e-14  # Local readouts carry no remote-setting signal.
        return dict(phase_sum=total, probabilities=tuple(map(float,marginal)),
                    local_marginals_uniform=True, coherent=coherent)

    def allocations(total, parts):
        if parts == 1:
            yield (total,)
        else:
            for first in range(total+1):
                for rest in allocations(total-first,parts-1):
                    yield (first,)+rest

    sibling_results, allocation_checks, tensor_grid_checks = {}, 0, 0
    for j in (2,3,4,5):
        fixture = prepare_siblings(j)
        one = sibling_readout((refined_x,)+(Q(0),)*(j-1))
        all_ports = sibling_readout((refined_x,)*j)
        same_budget = sibling_readout((j*refined_x,)+(Q(0),)*(j-1))
        erased = sibling_readout((refined_x,)*j,coherent=False)
        assert np.max(np.abs(np.array(all_ports['probabilities'])-same_budget['probabilities'])) < 5e-14
        for k in range(27):
            sibling_readout((Q(k,27),)*j)
            tensor_grid_checks += 1
        # Equal TOTAL number of finite local setting increments; no nontriadic x/j.
        for weights in allocations(7,j):
            tested = sibling_readout(tuple(Q(w,27) for w in weights))
            assert tested['phase_sum'] == Q(7,27)
            allocation_checks += 1
        sibling_results[j] = dict(fixture=fixture,one_port=one,one_unit_on_each=all_ports,
            same_total_on_one_port=same_budget,phase_erased=erased,
            total_local_setting_units=j,setting_depth=3,
            new_recording_work_for_settings=None)
    # Direct match to manuscript Eq56 after the explicit outcome/orientation relabeling.
    x,y = Q(1,9),Q(2,27)
    m = sibling_readout((x,-y))
    for a in range(3):
        for b in range(3):
            assert abs(m['probabilities'][(b-a)%3]/3-kernel(x-y+Q(a-b,3))/3) < 2e-14

    # C. Does finite phase compatibility select S versus spectral |S|? No:
    # integer mode charges make BOTH represent the same phase-clock group.
    # beta is an algebraic period normalization, NOT a duration/energy calibration.
    beta = 1
    for n in (3,4,5):
        beta = beta*(n-1)//gcd(beta,n-1)
    assert beta == 12
    cavities = {}
    representation_checks = 0
    for n in (3,4,5):
        P0 = np.ones((n,n))/n
        Pc = np.eye(n)-P0
        S = -(np.ones((n,n))-np.eye(n))/(n-1)
        eig, vec = np.linalg.eigh(S)
        magnitude = (vec*np.abs(eig))@vec.T
        assert np.max(np.abs(S-(-P0+Pc/(n-1)))) < 2e-14
        assert np.max(np.abs(magnitude-(P0+Pc/(n-1)))) < 2e-14
        result = {}
        for name, common in (('signed',-beta),('magnitude',beta)):
            contrast = beta//(n-1)
            def gate(x):
                a,b = phases((common,contrast),x)
                return a*P0+b*Pc
            charges = (common,)+(contrast,)*(n-1)
            gap = abs(common-contrast)
            for depth in range(1,5):
                q = 3**depth
                for k in range(q):
                    z=Q(k,q); U=gate(z)
                    assert np.max(np.abs(U.conj().T@U-np.eye(n))) < 5e-14
                    assert np.max(np.abs(U-gate(Q(3*k,3*q)))) < 2e-14
                    assert np.max(np.abs(gate(z)@gate(Q(1,q))-gate(z+Q(1,q)))) < 5e-14
                    representation_checks += 1
            V=gate(refined_x)
            returned=float(abs(V[0,0])**2)
            formula=1-4*(n-1)/n**2*(np.sin(pi*float((gap*refined_x)%1))**2)
            assert abs(returned-formula) < 3e-14
            # Candidate permutation symmetry is retained by either map.
            for k in range(n-1):
                order=list(range(n)); order[k],order[k+1]=order[k+1],order[k]
                assert np.max(np.abs(V[np.ix_(order,order)]-V)) < 2e-14
            result[name]=dict(integer_mode_charges=charges,phase_gap=gap,
                unscaled_gap=Q(gap,beta),vertex_return_at_1_27=returned,
                triadic_refinement_compatible=True,compiled_native_gate_word=False)
        cavities[n]=result
    assert cavities[5]['signed']['unscaled_gap']/cavities[3]['signed']['unscaled_gap'] == Q(5,6)
    assert cavities[5]['magnitude']['unscaled_gap']/cavities[3]['magnitude']['unscaled_gap'] == Q(3,2)

    # Same statement applies to the actual registry's rational eigenvalue inventory.
    exact_eigenvalues=(Q(-1),-Q(62,73),Q(383,4964),Q(1,73),Q(1,136))
    beta_registry=1
    for v in exact_eigenvalues:
        beta_registry=beta_registry*v.denominator//gcd(beta_registry,v.denominator)
    assert beta_registry == 9928
    registry_charges={name:tuple(int(beta_registry*v) for v in values) for name,values in (
        ('signed',exact_eigenvalues),('magnitude',tuple(abs(v) for v in exact_eigenvalues)))}
    saved_eigenvalues=sorted(old['quantum']['symmetric_generator_spectrum'])
    expected_eigenvalues=sorted([float(exact_eigenvalues[0]),float(exact_eigenvalues[1]),
        float(exact_eigenvalues[2])]+[1/73]*124+[1/136]*10)
    assert np.max(np.abs(np.array(saved_eigenvalues)-expected_eigenvalues)) < 3e-14
    for charges in registry_charges.values():
        for k in range(81):
            assert all(Q(q*k,81) == Q(q*3*k,243) for q in charges)
    assert (old,_reference) == before

    print('CELL 24 — FINITE RECORD PHASES, NOT A PHYSICAL MASS/ENERGY CALIBRATION')
    print(f'Local source phase D(x)=diag(1,e^(2pi ix),e^(4pi ix)); {clock_checks} grid settings,')
    print(f'  {group_checks} exact composition checks; same coarse tick survives refinement.')
    print('Fixed local return: L(x)=1/3+(4/9)cos(2pi x)+(2/9)cos(4pi x).')
    print('\nCoarse-generator alias: (0,1,2) and (0,1,11) agree on EVERY k/9 setting.')
    print('At the first finer setting x=1/27:')
    print('  source labels:',alias['primary_probabilities'])
    print('  aliased labels:',alias['alias_probabilities'])
    print(f"  total variation={alias['refined_total_variation']:.12f}; finer source family rejects alias.")
    print('\nSHARED-SOURCE SIBLINGS: outcome sum mod3, same local inverse-Fourier readout.')
    print(' j   min Gamma    one setting unit P(sum=0)   one unit on EACH P(sum=0)')
    for j,row in sibling_results.items():
        print(f" {j} {j+1:11d} {row['one_port']['probabilities'][0]:30.12f} "
              f"{row['one_unit_on_each']['probabilities'][0]:29.12f}")
    print('Each-unit and all-j-units-on-one-port results agree: response depends on SUM of phases.')
    print(f'{allocation_checks} allocations of seven total units and {tensor_grid_checks} tensor/grid cases checked.')
    print('Phase erasure yields (1/3,1/3,1/3); every one-factor marginal stays uniform.')
    print('j counts siblings, NOT K_j order. These j>=2 burst fixtures do NOT occur at Gamma=2.')
    print('Birth work is stored; phase settings have NO supplied native workload or time assignment.')
    print('\nCLOCK-ADMISSIBLE CAVITY COUNTERMODELS: common algebraic normalization beta=12.')
    print(' n    signed gap    magnitude gap    P_return signed/magnitude at x=1/27')
    for n,row in cavities.items():
        print(f" {n} {str(row['signed']['unscaled_gap']):>13s} {str(row['magnitude']['unscaled_gap']):>16s} "
              f"{row['signed']['vertex_return_at_1_27']:.12f} / {row['magnitude']['vertex_return_at_1_27']:.12f}")
    print('K5/K3 gaps: signed 5/6; magnitude 3/2. BOTH pass phase/refinement/permutation checks.')
    print(f'{representation_checks} cavity unitary/refinement/composition checks. No gate synthesis claimed.')
    print(f'Actual registry eigenvalue clock normalizer={beta_registry}; charges={registry_charges}.')
    print('A finite clock element has torsion: no nonzero additive map of that group into real elapsed time.')
    print('An unwrapped event count or a separately specified timing law is additional information.')
    print('PASS: prior Cell23/registry unchanged, legal sibling fixtures, exact phase identities, tensor Born readout.')
    print('Selected: phase law as a function of supplied settings. NOT selected: S vs |S| as physical energy.')
    return dict(protocol='declared finite phase family + sibling tensor readout + generator identifiability audit',
        local_phase_labels=(0,1,2),local_probability_harmonics=dict(constant=Q(1,3),cos1=Q(4,9),cos2=Q(2,9)),
        clock_settings_checked=clock_checks,exact_composition_checks=group_checks,coarse_alias=alias,
        siblings=sibling_results,allocation_checks=allocation_checks,tensor_grid_checks=tensor_grid_checks,
        cavity_normalizer=beta,cavities=cavities,cavity_representation_checks=representation_checks,
        registry_normalizer=beta_registry,registry_integer_charges=registry_charges,
        helpers=dict(triadic=triadic,phases=phases,kernel=kernel,ramsey=ramsey,sibling_readout=sibling_readout),
        selected_physical_energy_generator=None,selected_phase_to_service_map=None,physical_calibration=None,
        K_n_identified_with_n_siblings=False,native_dynamics_changed=False,quantum_postulates_derived=False,
        empirical_data_tested=False,source_parameters_changed=False)


dcu_mass_24 = _run_dcu_mass_24()
