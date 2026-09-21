"""Standalone Cell34 replay. Python3.10+ and NumPy; no network or notebook required."""
from pathlib import Path
import contextlib,gzip,hashlib,json,runpy,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parent

def plain(x):
    if isinstance(x,np.ndarray):return plain(x.tolist())
    if isinstance(x,np.generic):return plain(x.item())
    if isinstance(x,complex):return {'real':x.real,'imag':x.imag}
    if isinstance(x,dict):return {str(k):plain(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [plain(v) for v in x]
    return x

def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(u,v) in enumerate(zip(a,b)):compare(u,v,path+f'[{i}]')
    elif isinstance(a,float):
        assert abs(a-b)<=3e-12*max(1,abs(a),abs(b)),(path,a,b)
    else:assert a==b,(path,a,b)

def main():
    start=time.perf_counter()
    with gzip.open(ROOT/'reference/DCU_Mass_Cell_26_RESULTS.json.gz','rt') as f:inputs=json.load(f)
    with open(ROOT/'REPLAY_output.txt','w',encoding='utf-8') as log,contextlib.redirect_stdout(log):
        env=runpy.run_path(str(ROOT/'DCU_Mass_Cell_34.py'),init_globals={'dcu_mass_26':inputs})
    result=plain(env['dcu_mass_34'])
    with gzip.open(ROOT/'REPLAY_RESULTS.json.gz','wt',encoding='utf-8') as f:json.dump(result,f,separators=(',',':'),allow_nan=False)
    expected=ROOT/'DCU_Mass_Cell_34_RESULTS.json.gz'
    if expected.exists():
        with gzip.open(expected,'rt') as f:compare(result,json.load(f))
        print('PASS: every serialized field reproduced (platform numerical tolerance only).')
    print(f'Cell34 elapsed: {time.perf_counter()-start:.3f} s. Input unchanged; no native simulations run.')
    return result
if __name__=='__main__':main()
