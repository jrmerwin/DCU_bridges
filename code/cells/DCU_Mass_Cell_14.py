# CELL 14 — exact global spectrum and conditional cross-sector formula closure.
# Run after Cell 1 (_reference). Cells 6–13 are NOT required.
# Standard library only. No geometry, dynamics, calibration or mass fit is changed.
from collections import Counter
from fractions import Fraction
from copy import deepcopy


def _run_dcu_mass_14():
    if '_reference' not in globals():
        raise RuntimeError('Run Cell 1 first to construct the complete _reference.')
    ref = _reference
    before = deepcopy(ref)
    rows, ancestors = ref['rows'], ref['ancestors']
    adjacency = tuple(frozenset(a) for a in ref['adjacency'])
    N = len(rows)
    counts = Counter(r['sector'] for r in rows)
    degree = tuple(map(len, adjacency))
    assert N == 137 and counts == Counter(S=81, I=40, G=16)

    # Partition by actual ancestry, NOT by a target eigenvalue or mass formula.
    c = ref['pairs'][(0, 1)]
    pa = ref['pairs'][tuple(sorted((0, c)))]
    pb = ref['pairs'][tuple(sorted((1, c)))]
    groups = {tag: [] for tag in ('U', 'A', 'B')}
    for i, row in enumerate(rows):
        a = ancestors[row['object_id']]
        left, right = pa in a, pb in a
        assert left or right
        groups['U' if left and right else 'A' if left else 'B'].append(i)
    tags = ('U', 'A', 'B')
    blocks = [tuple(groups[t]) for t in tags]
    u, b, b2 = map(len, blocks)
    assert (u, b, b2) == (11, 63, 63)
    block_of = {i: t for t, block in enumerate(blocks) for i in block}
    for i in range(N):
        expected = {j for j in range(N) if j != i and
                    (block_of[i] == block_of[j] or 0 in (block_of[i], block_of[j]))}
        assert adjacency[i] == expected
    dU, dB = degree[blocks[0][0]], degree[blocks[1][0]]
    assert (dU, dB) == (N-1, b-1+u)

    def transfer(vector):
        # Previously DECLARED overlap readout T=-A D^{-1}; not native service.
        return tuple(-sum((Fraction(vector[j], degree[j]) for j in adjacency[i]),
                          Fraction(0)) for i in range(N))

    def check(vector, eigenvalue):
        assert transfer(vector) == tuple(eigenvalue*x for x in vector)

    spectrum = Counter()
    checked_vectors = 0
    # Within-block zero-sum modes. These span N-3 independent directions.
    for block in blocks:
        eigenvalue = Fraction(1, degree[block[0]])
        for j in block[1:]:
            vector = [0]*N
            vector[block[0]], vector[j] = 1, -1
            check(vector, eigenvalue)
            spectrum[eigenvalue] += 1
            checked_vectors += 1

    # Three block-constant modes. The formulas use measured block sizes/degrees.
    third = 1-Fraction(u-1, dU)-Fraction(b-1, dB)
    modes = ((Fraction(-1), (dU, dB, dB)),
             (-Fraction(b-1, dB), (0, 1, -1)),
             (third, (-2*b, u, u)))
    for eigenvalue, amplitudes in modes:
        vector = [amplitudes[block_of[i]] for i in range(N)]
        check(vector, eigenvalue)
        spectrum[eigenvalue] += 1
        checked_vectors += 1
    # The three eigenvalues are distinct, hence these macro modes are independent.
    # Block-constant vectors intersect the within-block-zero-sum space only at 0.
    assert len({v for v, _ in modes}) == 3 and checked_vectors == N
    assert sum(value*mult for value, mult in spectrum.items()) == 0
    trace_square = sum((Fraction(1, degree[i]*degree[j])
                        for i in range(N) for j in adjacency[i]), Fraction(0))
    assert sum(value**2*mult for value, mult in spectrum.items()) == trace_square
    quotient = tuple(tuple(-sum((Fraction(1, degree[j])
                                  for j in block if j in adjacency[representative]), Fraction(0))
                           for block in blocks)
                     for representative in (x[0] for x in blocks))

    # Separate SOURCE PRESCRIPTIONS: NOT outputs inferred from the overlap spectrum.
    # These combine the March lepton formula and the February taxonomy formulas.
    # Conditional subset only: this does NOT assert all versions of the papers agree.
    substrate = N-1                         # Requires the proposed 1-interface dictionary.
    f3, f5 = 1-Fraction(1, 3-1), 1-Fraction(1, 5-1)
    r_mu = N*f5/f3                          # RMR lepton prescription, Eq. (9).
    omega = Fraction(counts['I'], substrate) # Taxonomy surface-to-substrate identification.
    solar = Fraction(1, 4*N*substrate**6)    # Delta m21^2 / m_e^2, not an absolute eV scale.
    r_nu = 3**2+5**2-1
    # Eliminate N. These are exact consequences, NOT independent predictions/data.
    assert omega*(2*r_mu-3) == 3*counts['I']
    assert solar*8*r_mu*(2*r_mu-3)**6 == 3**7
    assert solar == omega**7/(4*counts['I']**6*(omega+counts['I']))
    assert ref == before

    print('CELL 14 — EXACT STRUCTURE, DECLARED OPERATOR, CONDITIONAL PHYSICS KEPT SEPARATE')
    print(f'Registry N={N}; S/I/G={tuple(counts[s] for s in ("S", "I", "G"))}; '
          f'blocks U/A/B={(u, b, b2)}; edges={sum(degree)//2}.')
    print('\nT=-A D^-1 on the ACTUAL registry-overlap graph:')
    print('   eigenvalue         multiplicity')
    for value, mult in sorted(spectrum.items()):
        print(f'   {str(value):>14s} {mult:20d}')
    print(f'PASS: {checked_vectors} independent exact eigenvectors; trace and trace(T^2) checked.')
    print('This operator is a declared readout, NOT the native service transition.')
    print('\nSeparate fixed RMR/taxonomy prescriptions (not derived by this cell):')
    print(f'   muon/electron = N*f5/f3 = {r_mu}')
    print(f'   Omega_m = I/(N-1) = {omega}')
    print(f'   Delta m21^2 / m_e^2 = {solar}')
    print(f'   R_nu = 3^2+5^2-1 = {r_nu}')
    print('   A full neutrino spectrum additionally assumes normal ordering and m1=0.')
    print('\nPASS: exact N-eliminated consistency identities:')
    print(f'   Omega_m*(2*r_mu-3) = {3*counts["I"]}')
    print(f'   (Delta m21^2/m_e^2)*8*r_mu*(2*r_mu-3)^6 = {3**7}')
    print('   Delta m21^2/m_e^2 = Omega_m^7 / [4*I^6*(Omega_m+I)]')
    print('Identities are algebraic closure, NOT independent experimental confirmations.')
    print('No empirical targets, SI scale, exclusions of failed formulas, or fitted corrections used.')
    print('The 1+136 mapping, physical observable assignments and depth exponent remain hypotheses.')
    print('Previous reference, geometry results and physical calibrations are unchanged.')
    return dict(protocol='exact overlap spectrum + explicitly conditional formula closure',
                spectrum=dict(spectrum), block_quotient=quotient,
                block_sizes=dict(zip(tags, (u, b, b2))), eigenvectors_checked=checked_vectors,
                prescriptions=dict(muon_electron=r_mu, matter_fraction=omega,
                                   solar_splitting_over_electron_mass_squared=solar,
                                   neutrino_splitting_ratio=r_nu),
                empirical_comparison_performed=False,
                native_dynamics_changed=False, physical_calibration_changed=False)


dcu_mass_14 = _run_dcu_mass_14()
