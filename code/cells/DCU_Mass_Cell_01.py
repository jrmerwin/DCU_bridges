# CELL 1 — matched-support native recording response (no physical mass assigned).
# Put the original reproducibility ZIP beside the notebook, or edit ARCHIVE.
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
import ast, gzip, json, zipfile

ARCHIVE = Path('DCU_Registry_Manuscript_and_Reproducibility.zip')
if not ARCHIVE.is_file():
    raise FileNotFoundError(f'Place the original ZIP beside the notebook, or edit ARCHIVE: {ARCHIVE.resolve()}')

# Load saved observations and original definitions, WITHOUT notebook launch blocks.
_ns = {'__name__': 'dcu_mass_cell_1'}
with zipfile.ZipFile(ARCHIVE) as _zip:
    _prefixes = [n[:-len('code/native_constructor.py')] for n in _zip.namelist()
                 if n.endswith('code/native_constructor.py')]
    if len(_prefixes) != 1:
        raise ValueError('Expected one native_constructor.py in the manuscript ZIP.')
    _prefix = _prefixes[0]
    for _file in ('native_constructor.py', 'DCU_Rosetta_Calibration_F_Registry_Anatomy.py'):
        _tree = ast.parse(_zip.read(_prefix + 'code/' + _file).decode('utf-8'))
        _defs = [n for n in _tree.body if isinstance(n, (ast.Import, ast.ImportFrom,
                                                      ast.FunctionDef, ast.ClassDef))]
        exec(compile(ast.Module(body=_defs, type_ignores=[]), _file, 'exec'), _ns)
    _obs = json.loads(gzip.decompress(_zip.read(_prefix + 'results/reproduced_exact_observations.json.gz')))
    _panel = json.loads(gzip.decompress(_zip.read(_prefix + 'results/reproduced_five_family_full.json.gz')))

_Native = _ns['DCUStructure']
_reference = _ns['rosetta_dag7_reference']()  # Complete 137-entry dictionary; G multiplicity = 8.
_sector = {r['term']: r['sector'] for r in _reference['rows']}
_parents = _obs['parents']
assert _obs['observer_seed'] == _panel['source_observer_seed'] == 20260920
assert len(_obs['samples']) == 512 and len(_panel['motifs']) == 14910

@lru_cache(None)
def _term(v):
    return 'ab'[v] if v < 2 else '(' + '|'.join(sorted(_term(p) for p in _parents[v])) + ')'

# Fix selection before measuring any workload: first two by primitive-labelled
# canonical parentage in the existing U12 support-size 11 / birth-grade 4 stratum.
_selected = []
for _family in ('K3', 'K4'):
    _report = next(r for r in _panel['round_reports'] if r['round'] == 12 and r['family'] == _family)
    _ids = [i for i in _report['unique_motif_ids'] if
            (_panel['motifs'][i]['support_size'], _panel['motifs'][i]['birth_round']) == (11, 4)]
    _ids.sort(key=lambda i: tuple(sorted(_term(v) for v in _panel['motifs'][i]['core'])))
    if len(_ids) < 2:
        raise ValueError(f'Insufficient archived {_family} instances in the declared stratum.')
    _selected.extend((f'{_family}-{j+1}', i, _panel['motifs'][i]) for j, i in enumerate(_ids[:2]))

# Native accounting regression: the manuscript's one-event and two-event examples.
_fixture = _Native()
_fixture.add_batch([(0, 1)])
_fixture.add_batch([(0, 2)])
assert _fixture.recording_cost([(2, 3)]) == 48
assert _fixture.recording_cost([(0, 3), (1, 3)]) == 60


def _measure(motif, branch):
    """Fresh minimal-support preparation for each probe context; queries do not commit."""
    s, remap = _Native(), {0: 0, 1: 1}
    for old in motif['support']:
        if old >= 2:
            born, _ = s.add_batch([tuple(remap[p] for p in _parents[old])])
            remap[old] = born[0]
    support = {remap[v] for v in motif['support']}
    core_support = set().union(*(s.ancestors[remap[v]] for v in motif['core']))
    children = {remap[v] for v in motif['edge_children']}
    assert support == core_support | children and len(support) == 11
    for (u, v), w in zip(motif['required_edges'], motif['edge_children']):
        assert set(s.parents[remap[w]]) == {remap[u], remap[v]}

    # Matched probes: p0={a,b}; pj={p(j-1),a} or {p(j-1),b}, j=1,...,5.
    # Both have grade 6, support size 8 and L=12. Reuse an existing term, never duplicate it.
    p = s.pair_to_id[(0, 1)]
    for _ in range(5):
        pair = tuple(sorted((p, branch)))
        if pair in s.pair_to_id:
            p = s.pair_to_id[pair]
        else:
            born, _ = s.add_batch([pair])
            p = born[0]
    height = [0, 0]
    for pair in s.parents[2:]:
        height.append(1 + max(height[v] for v in pair))
    assert (height[p], len(s.ancestors[p]), s.path_weight(p)) == (6, 8, 12)
    assert s.recorded == {v for pair in s.parents if pair is not None for v in pair}
    assert len(s) == len(support | s.ancestors[p])  # Cache is NOT used as an ambient universe.

    # Both arms start from THIS SAME pre-probe state. The baseline is the legal
    # comparison {a,p}; every target arm is one legal {u,p}, NOT a simultaneous burst.
    baseline = s.recording_cost([(0, p)])
    reg = {remap[v]: _sector[_term(v)] for v in motif['support'] if _term(v) in _sector}
    assert tuple(list(reg.values()).count(q) for q in ('S', 'I', 'G')) == tuple(motif['full_SIG'])
    ports, parts = [], {q: 0 for q in ('S', 'I', 'G', 'nonregistry')}
    before = (len(s), frozenset(s.recorded))
    for old in motif['support']:
        if old < 2:
            continue
        u = remap[old]
        total = s.recording_cost([(u, p)])
        increment = total - baseline
        exclusive = s.ancestors[u] - s.ancestors[p]
        charges = {z: (2 if z in s.recorded else 11) * s.path_weight(z) for z in exclusive}
        assert increment == sum(charges.values()) >= 0
        for z, charge in charges.items():
            parts[reg.get(z, 'nonregistry')] += charge
        ports.append({'cache_port': old, 'total_work': total, 'baseline_work': baseline,
                      'increment': increment, 'target_recorded': u in s.recorded})
    assert before == (len(s), frozenset(s.recorded)) and len(ports) == 9
    response = Fraction(sum(v['increment'] for v in ports), len(ports))
    return {'X': response, 'parts': {q: Fraction(v, len(ports)) for q, v in parts.items()},
            'probe_support_overlap': len(support & s.ancestors[p]), 'ports': ports}


# A first; freeze ONE reporting unit before B. K4 is an operational control,
# NOT a known vacuum, a zero-mass anchor, or an assigned massive species.
_results = []
for _branch, _probe in enumerate(('A', 'B')):
    for _name, _id, _motif in _selected:
        _row = _measure(_motif, _branch)
        _row.update(name=_name, archive_motif_id=_id, probe=_probe,
                    core=_motif['core'], support=_motif['support'], full_SIG=_motif['full_SIG'])
        _results.append(_row)
    if _probe == 'A':
        _unit = next(r['X'] for r in _results if r['name'] == 'K4-1' and r['probe'] == 'A')
for _r in _results:
    _r['fixed_unit_ratio'] = _r['X'] / _unit if _unit else None

# Retain exact rational results and all individual counterfactual comparisons in memory.
dcu_mass_1 = {'protocol': 'mean extra native work per nonprimitive support port; matched minimal ideals',
              'selection': 'U12; D=11,h=4; first two K3/K4 by labelled parentage',
              'unit': _unit, 'unit_source': 'K4-1/A', 'rows': _results}
print('Matched candidates (cache IDs are labels, not physical identities):')
for _name, _id, _m in _selected:
    print(f"  {_name}: archive {_id}, core={_m['core']}, S/I/G={_m['full_SIG']}")
print('\nX = exact mean extra work over 9 ports; one comparison per fresh-state arm.')
print(f'Fixed reporting unit: K4-1/A = {_unit}; NOT a physical mass calibration.')
print(' case   probe       X exact    X/fixed unit    registry part    overlap')
for _r in _results:
    _ratio = 'undefined' if _r['fixed_unit_ratio'] is None else f"{float(_r['fixed_unit_ratio']):.6f}"
    _regpart = sum(_r['parts'][q] for q in ('S', 'I', 'G'))
    print(f" {_r['name']:6s} {_r['probe']:>3s} {str(_r['X']):>13s} {_ratio:>15s}"
          f" {str(_regpart):>16s} {_r['probe_support_overlap']:>10d}")
print('\nProbe-context change (B/A - 1), with no refitting:')
for _name, _, _ in _selected:
    _a, _b = [next(r['X'] for r in _results if r['name'] == _name and r['probe'] == p) for p in ('A', 'B')]
    print(f"  {_name}: {'undefined' if not _a else f'{100*float(_b/_a-1):+.3f}%'}")
print('\nPASS: native cost examples, legal pairs, complete supports, recorded masks,')
print('      equal probe anatomy, exact baseline subtraction, and read-only measurements.')
print('This is a conditional recording response, not yet mass or inertia.')
