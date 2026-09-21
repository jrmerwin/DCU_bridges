#!/usr/bin/env python3
"""Offline scientific replay; original reference files are never overwritten.

Default compact: original horizontal Cells 14/17-24; channels 32-34 and
instruments 36-37; new Cells 38-39; the saved-record epoch-study analysis.
Use --suite native to add the original stochastic Cells 25/26/29-31.
Full million-tick generation is deliberately a separate opt-in command.
"""
from __future__ import annotations
import argparse, ast, contextlib, gzip, json, math, os, platform, subprocess, sys, time
from fractions import Fraction
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]

def load_json(path):
    with (gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz') else open(path, encoding='utf-8')) as f:
        return json.load(f)

def plain(v):
    if isinstance(v, Fraction): return str(v)
    if isinstance(v, np.ndarray): return plain(v.tolist())
    if isinstance(v, np.generic): return plain(v.item())
    if isinstance(v, complex): return {'real':v.real,'imag':v.imag}
    if isinstance(v, dict): return {str(k):plain(x) for k,x in v.items() if not callable(x)}
    if isinstance(v, (list,tuple)): return [plain(x) for x in v]
    if isinstance(v, (set,frozenset)): return sorted(plain(x) for x in v)
    return v

def definitions(path, env):
    tree=ast.parse(path.read_text(encoding='utf-8'))
    nodes=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom,ast.FunctionDef,ast.ClassDef))]
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),env)

def compare(expected, actual, path='root', counts=None):
    """Check exact arithmetic and numeric fields across the legacy JSON encodings."""
    if counts is None: counts={'exact_leaves':0,'floating_leaves':0,'archive_only_fields':[]}
    if isinstance(expected,dict) and set(expected).issuperset({'numerator','denominator'}) and isinstance(actual,Fraction):
        if actual!=Fraction(expected['numerator'],expected['denominator']): raise AssertionError(path)
        counts['exact_leaves']+=1
        if 'decimal' in expected: compare(expected['decimal'],float(actual),path+'.decimal',counts)
    elif isinstance(expected,dict):
        if isinstance(actual,complex): actual={'real':actual.real,'imag':actual.imag}
        if not isinstance(actual,dict): raise AssertionError((path,type(actual)))
        amap={str(k):v for k,v in actual.items() if not callable(v)}
        for k,v in expected.items():
            if k not in amap:
                if k=='validation' or k=='helpers' and not v:
                    counts['archive_only_fields'].append(path+'.'+k); continue
                raise AssertionError(('missing field',path,k))
            compare(v,amap[k],path+'.'+k,counts)
    elif isinstance(expected,list):
        if isinstance(actual,np.ndarray):actual=actual.tolist()
        if len(expected)!=len(actual):raise AssertionError(('length',path))
        for i,(a,b) in enumerate(zip(expected,actual)):compare(a,b,f'{path}[{i}]',counts)
    elif isinstance(actual,Fraction) and isinstance(expected,str):
        if actual!=Fraction(expected):raise AssertionError(('fraction',path,expected,str(actual)))
        counts['exact_leaves']+=1
    elif isinstance(expected,(float,np.floating)):
        b=float(actual)
        if not math.isclose(float(expected),b,rel_tol=3e-10,abs_tol=3e-11):
            raise AssertionError(('float',path,expected,b))
        counts['floating_leaves']+=1
    else:
        if isinstance(actual,np.generic):actual=actual.item()
        if expected!=actual:raise AssertionError(('exact',path,expected,actual))
        counts['exact_leaves']+=1
    return counts

def saved(n):
    for ext in ('.json','.json.gz'):
        p=ROOT/'data/reference'/f'DCU_Mass_Cell_{n:02}_RESULTS{ext}'
        if p.exists():return load_json(p)
    return None

def run_cell(n,env,out):
    p=ROOT/'code/cells'/f'DCU_Mass_Cell_{n:02}.py'
    start=time.perf_counter()
    with (out/f'cell_{n:02}.log').open('w',encoding='utf-8') as log,contextlib.redirect_stdout(log):
        exec(compile(p.read_text(encoding='utf-8'),str(p),'exec'),env)
    result=env[f'dcu_mass_{n:02}']; expected=saved(n)
    checks=compare(expected,result) if expected is not None else {'embedded_assertions':'passed; no separate result JSON'}
    with gzip.GzipFile(filename=str(out/f'cell_{n:02}_replayed.json.gz'),mode='wb',mtime=0) as f:
        f.write((json.dumps(plain(result),sort_keys=True,allow_nan=False,separators=(',',':'))+'\n').encode())
    item={'cell':n,'seconds':time.perf_counter()-start,'checks':checks}
    print('PASS cell',n,f"({item['seconds']:.2f} s)",flush=True)
    return item

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--suite',choices=['compact','new','horizontal','native'],default='compact')
    ap.add_argument('--output',type=Path,default=ROOT/'build/replay')
    args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();previous=Path.cwd();os.chdir(out)
    report={'suite':args.suite,'python':platform.python_version(),'numpy':np.__version__,'cells':[],
            'scope':'Offline finite calculations; archived native histories are inputs to channel replays. No late neutron corpus generated.'}
    try:
        env={'__name__':'dcu_reproduction'}
        if args.suite in ('compact','horizontal','native'):
            definitions(ROOT/'code/registry/DCU_Rosetta_Calibration_F_Registry_Anatomy.py',env)
            definitions(ROOT/'code/registry/native_constructor.py',env)
            env['_reference']=env['rosetta_dag7_reference']()
            env['_Native']=env['DCUStructure']
            for n in (14,17,18,19,20,21,22,23,24):report['cells'].append(run_cell(n,env,out))
        if args.suite=='native':
            for n in (25,26,29,30,31):report['cells'].append(run_cell(n,env,out))
        if args.suite in ('compact','native'):
            if args.suite!='native':env.update(dcu_mass_29=saved(29),dcu_mass_26=saved(26))
            for n in (32,33,34,36,37):report['cells'].append(run_cell(n,env,out))
            script=ROOT/'experiments/cell35/replay_readout.py'
            proc=subprocess.run([sys.executable,str(script)],cwd=out,text=True,capture_output=True)
            (out/'cell_35_analysis.log').write_text(proc.stdout+proc.stderr,encoding='utf-8')
            if proc.returncode: raise RuntimeError('Cell35 analysis replay failed; see log')
            report['epoch_study']='Original compact record/readout replay passed; no new million-tick generation.'
            print('PASS cell35 saved-record analysis',flush=True)
        if args.suite in ('compact','native','new'):
            for n in (38,39):report['cells'].append(run_cell(n,env,out))
        report.update(status='PASS',elapsed_seconds=time.perf_counter()-start)
    finally:
        os.chdir(previous)
        (out/'REPLAY_REPORT.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Report:',out/'REPLAY_REPORT.json')
if __name__=='__main__':main()
