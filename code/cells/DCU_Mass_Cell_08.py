# CELL 8 — metric compatibility, NOT a mass or a derivation of three-space.
# Run after Cell 5; Cells 6 and 7 are not required. Standard library only.
# Added test premise: all core-edge witnesses prescribe the SAME unit length.
# Countercontrol: the SAME graph with distances from explicit 3D coordinates.
from collections import Counter
from fractions import Fraction
from itertools import combinations


def _run_dcu_mass_8():
    needed = ('dcu_mass_5', '_panel', '_parents', '_reference')
    absent = [k for k in needed if k not in globals()]
    if absent:
        raise RuntimeError('Run Cells 1-5 first. Missing: ' + ', '.join(absent))

    def rank(matrix):
        a = [list(map(Fraction, row)) for row in matrix]
        if not a:
            return 0
        row = 0
        for col in range(len(a[0])):
            pivot = next((i for i in range(row, len(a)) if a[i][col]), None)
            if pivot is None:
                continue
            a[row], a[pivot] = a[pivot], a[row]
            value = a[row][col]
            a[row] = [x/value for x in a[row]]
            for i in range(row+1, len(a)):
                value = a[i][col]
                a[i] = [x-value*y for x, y in zip(a[i], a[row])]
            row += 1
            if row == len(a):
                break
        return row

    def determinant(matrix):
        a = [list(map(Fraction, row)) for row in matrix]
        value = Fraction(1)
        for j in range(len(a)):
            pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
            if pivot is None:
                return Fraction(0)
            if pivot != j:
                a[j], a[pivot] = a[pivot], a[j]
                value = -value
            p = a[j][j]
            value *= p
            for i in range(j+1, len(a)):
                factor = a[i][j]/p
                a[i] = [x-factor*y for x, y in zip(a[i], a[j])]
        return value

    def metric_audit(squared_distances):
        """Exact Euclidean feasibility for a complete distance table, n<=5."""
        d = [list(map(Fraction, row)) for row in squared_distances]
        n = len(d)
        if not 2 <= n <= 5 or any(len(row) != n for row in d):
            raise ValueError('Supply a square 2-to-5-point squared-distance table.')
        if any(d[i][i] for i in range(n)) or any(
                d[i][j] < 0 or d[i][j] != d[j][i] for i in range(n) for j in range(n)):
            raise ValueError('Squared distances must be symmetric, nonnegative, and zero on the diagonal.')
        means = [sum(row)/n for row in d]
        grand = sum(means)/n
        b = [[-(d[i][j]-means[i]-means[j]+grand)/2 for j in range(n)] for i in range(n)]
        # A symmetric matrix is PSD iff EVERY principal minor is nonnegative.
        psd = all(determinant([[b[i][j] for j in subset] for i in subset]) >= 0
                  for k in range(1, n+1) for subset in combinations(range(n), k))
        r = rank(b)
        cm = [[0] + [1]*n] + [[1] + row for row in d]
        return dict(gram=b, positive_semidefinite=psd, rank=r,
                    minimum_dimension=r if psd else None, fits_R3=psd and r <= 3,
                    cayley_menger_determinant=determinant(cm))

    # Actual S overlap graph. No coordinates, volumes, or lengths inferred here.
    ref = _reference
    s_roots = [r['object_id'] for r in ref['rows'] if r['sector'] == 'S']
    s_neighbors = {u: {v for v in s_roots if v != u and
                       len(ref['ancestors'][u] & ref['ancestors'][v]) > 3} for u in s_roots}
    degrees = dict(sorted(Counter(map(len, s_neighbors.values())).items()))
    n_edges = sum(map(len, s_neighbors.values()))//2
    assert len(s_roots) == 81 and degrees == {44: 72, 80: 9} and n_edges == 1944

    families = {}
    for n in (3, 4, 5):
        # Equal distances: this prescription is an INPUT, not a native length law.
        equal = [[int(i != j) for j in range(n)] for i in range(n)]
        regular = metric_audit(equal)
        assert regular['positive_semidefinite'] and regular['rank'] == n-1
        # Unequal-distance control on the 3D moment curve, not a DCU coordinate map.
        points = [(t, t*t, t*t*t) for t in range(n)]
        actual = [[sum((x-y)**2 for x, y in zip(a, b)) for b in points] for a in points]
        control = metric_audit(actual)
        assert control['fits_R3'] and control['rank'] == min(3, n-1)
        # Rigidity matrix: independent first-order squared-length constraints.
        rigidity = []
        for i, j in combinations(range(n), 2):
            row = [0]*(3*n)
            for axis in range(3):
                row[3*i+axis] = points[i][axis]-points[j][axis]
                row[3*j+axis] = -row[3*i+axis]
            rigidity.append(row)
        constraint_rank = rank(rigidity)
        assert constraint_rank == 3*n-6
        # For n>=4 no four control vertices are coplanar: nonincident straight
        # edges cannot intersect. This explicitly embeds the K5 GRAPH in 3D.
        for subset in combinations(points, 4):
            assert determinant([[b[k]-subset[0][k] for k in range(3)] for b in subset[1:]]) != 0
        families[n] = dict(equal_length=regular, unequal_3D_control=control,
                           independent_constraints=constraint_rank,
                           self_stress_dimension=len(rigidity)-constraint_rank)

    rows = []
    for record in dcu_mass_5['catalogue']:
        m = _panel['motifs'][record['archive_motif_id']]
        n = len(m['core'])
        assert m['family'] == f'K{n}'
        assert {tuple(sorted(e)) for e in m['required_edges']} == set(combinations(sorted(m['core']), 2))
        for edge, child in zip(m['required_edges'], m['edge_children']):
            assert set(_parents[child]) == set(edge)
        rows.append(dict(name=record['name'], archive_motif_id=record['archive_motif_id'],
                         core=tuple(m['core']), family_order=n))

    print('CELL 8 — EXACT CONDITIONAL GEOMETRY AUDIT; NO MASS CONVERSION')
    print(f'Actual S overlap graph: {len(s_roots)} roots; {n_edges} edges; degrees {degrees}.')
    print('Four ternary coordinates alone do NOT specify this adjacency or a 3D metric.')
    print('\n family  equal-edge dimension  fits R3?  unequal-control dimension  redundant constraints')
    for n, f in families.items():
        print(f" K{n} {f['equal_length']['rank']:20d} {str(f['equal_length']['fits_R3']):>9s}"
              f" {f['unequal_3D_control']['rank']:26d} {f['self_stress_dimension']:22d}")
    print(f'\nChecked actual core/witness edges of {len(rows)} unchanged catalogue candidates.')
    print('K5 Cayley-Menger determinant: equal unit edges =', families[5]['equal_length']['cayley_menger_determinant'],
          '; explicit 3D control =', families[5]['unequal_3D_control']['cayley_menger_determinant'])
    print('Equal-length K5 cannot fit Euclidean R3; K5 with compatible unequal lengths CAN.')
    print('One redundant constraint permits self-stress; it does not force nonzero stress.')
    print('The dimension 3 and equal preferred lengths are declared premises, NOT outputs of S.')
    print('No S-to-coordinate map, native strain cost, registry multiplier, or mass was supplied.')
    return dict(protocol='equal-edge metric premise versus explicit unequal-length 3D control',
                S_graph=dict(roots=len(s_roots), edges=n_edges, degrees=degrees),
                catalogue=rows, families=families, audit_distances=metric_audit,
                native_metric_identified=False, physical_mass_calibration=None)


dcu_mass_8 = _run_dcu_mass_8()
