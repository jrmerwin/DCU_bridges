"""Replay saved-history analysis. No network and no native evolution."""
from pathlib import Path
import ast, contextlib, gzip, json, time, math
import numpy as np
from fractions import Fraction
ROOT=Path(__file__).resolve().parent

def normalize(v):
    if isinstance(v,dict): return {str(k):normalize(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)): return [normalize(x) for x in v]
    if isinstance(v,np.ndarray): return normalize(v.tolist())
    if isinstance(v,(complex,np.complexfloating)):return {'real':float(v.real),'imag':float(v.imag)}
    if isinstance(v,Fraction): return str(v)
    if isinstance(v,np.generic):return v.item()
    return v

def definitions(path):
    tree=ast.parse(path.read_text(encoding='utf-8'))
    nodes=[v for v in tree.body if isinstance(v,(ast.Import,ast.ImportFrom,ast.FunctionDef,ast.ClassDef))]
    env={'__name__':'cell33_replay'}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),env)
    return env

def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert isinstance(b,dict) and set(a)==set(b),path
        for k in a: compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        assert isinstance(b,(int,float)) and math.isclose(a,b,rel_tol=2e-12,abs_tol=2e-12),(path,a,b)
    else:
        assert a==b,(path,a,b)

def main():
    t=time.perf_counter()
    saved=json.load(gzip.open(ROOT/'reference/DCU_Mass_Cell_29_RESULTS.json.gz','rt'))
    with (ROOT/'REPLAY_output.txt').open('w',encoding='utf-8') as f,contextlib.redirect_stdout(f):
        old=definitions(ROOT/'reference/DCU_Mass_Cell_32.py')['_run_dcu_mass_32'](saved)
        result=definitions(ROOT/'DCU_Mass_Cell_33.py')['_run_dcu_mass_33'](saved,old)
    plain=normalize(result)
    with gzip.GzipFile(filename=str(ROOT/'REPLAY_RESULTS.json.gz'),mode='wb',mtime=0) as f:
        f.write((json.dumps(plain,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode())
    print((ROOT/'REPLAY_output.txt').read_text())
    ref=ROOT/'DCU_Mass_Cell_33_RESULTS.json.gz'
    if ref.exists():
        expected=json.load(gzip.open(ref,'rt'))
        compare(expected,plain)
        print('PASS: all exact fields match; floating fields agree to numerical tolerance.')
        print('Bit-identical normalized scientific data:',expected==plain)
    print('Runtime seconds:',time.perf_counter()-t)
    return result
if __name__=='__main__': main()
