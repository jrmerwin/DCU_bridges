"""Source/checkpoint loader retained from Cell30; Cell31 helper extraction only."""
from pathlib import Path
from fractions import Fraction
from collections import Counter
from itertools import combinations
import ast, contextlib, gzip, json, time
ROOT = Path(__file__).resolve().parent

def context():
    env = {'__name__': 'cell30_reproduction'}
    tree = ast.parse((ROOT / 'reference/native_constructor.py').read_text())
    keep = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef, ast.ClassDef))]
    exec(compile(ast.Module(body=keep, type_ignores=[]), 'native_constructor.py', 'exec'), env)
    env['_Native'] = env['DCUStructure']
    tree = ast.parse((ROOT / 'reference/DCU_Mass_Cell_25.py').read_text())
    run = next((n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == '_run_dcu_mass_25'))
    step = next((n for n in run.body if isinstance(n, ast.FunctionDef) and n.name == 'step'))
    scope = {'Counter': Counter, 'combinations': combinations, 'gamma': 3}
    exec(compile(ast.Module(body=[step], type_ignores=[]), 'original_Cell25_native_step', 'exec'), scope)
    env['dcu_mass_25'] = {'Gamma': 3, 'm': 0, 'carriers': (3, 4), 'primary_depth': 3, 'helpers': {'native_step': scope['step']}}
    old = json.loads((ROOT / 'reference/Cell29_emission_checkpoints.json').read_text())

    class Draw:

        def __init__(self, pair, n):
            self.draw = list(pair) + [n]

        def sample(self, population, k):
            assert k == 3 and all((z in population for z in self.draw))
            return self.draw
    for name, entry in old['environments'].items():
        cp = entry['checkpoint']
        for row in entry['arms'].values():
            st = env['_Native']()
            for p in row['post_parentage'][2:cp['n']]:
                st.add_batch([tuple(p)])
            F, P, event = scope['step'](st, cp['F'], cp['P'], Draw(row['pair'], len(st)), entry['H'])
            assert F == row['post_F'] and P == row['post_P']
            assert clean(st.parents) == row['post_parentage']
            assert clean(event) == row['event']
    env['dcu_mass_29'] = old
    return env

def clean(obj):
    if isinstance(obj, Fraction):
        return str(obj)
    if isinstance(obj, complex):
        return dict(real=obj.real, imag=obj.imag)
    if isinstance(obj, dict):
        return {str(k): clean(v) for k, v in obj.items() if not callable(v)}
    if isinstance(obj, (tuple, list)):
        return [clean(v) for v in obj]
    if isinstance(obj, (set, frozenset)):
        return [clean(v) for v in sorted(obj)]
    if hasattr(obj, 'tolist'):
        return clean(obj.tolist())
    return obj


def inputs():
    env=context()
    with gzip.open(ROOT/'reference/DCU_Mass_Cell_30_RESULTS.json.gz','rt') as f:
        env['dcu_mass_30']=json.load(f)
    return env


def helpers():
    """Compile original Cell31 definitions without launching its full campaign."""
    env=inputs()
    tree=ast.parse((ROOT/'DCU_Mass_Cell_31.py').read_text())
    run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='_run_dcu_mass_31')
    cut=next(i for i,n in enumerate(run.body) if isinstance(n,ast.Assign) and
             any(isinstance(t,ast.Name) and t.id=='result' for t in n.targets))
    run.body=run.body[:cut]+[ast.Return(value=ast.Call(func=ast.Name(id='locals',ctx=ast.Load()),args=[],keywords=[]))]
    nodes=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom,ast.FunctionDef,ast.ClassDef))]
    exec(compile(ast.fix_missing_locations(ast.Module(body=nodes,type_ignores=[])),
                 'original_Cell31_definitions','exec'),env)
    return env,env['_run_dcu_mass_31']()
