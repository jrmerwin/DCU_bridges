"""Run Cell36 and compare every saved scientific field. NumPy + standard library."""
from pathlib import Path
import contextlib,gzip,json,math,runpy,time
ROOT=Path(__file__).resolve().parent
start=time.perf_counter()
with (ROOT/'Cell36_REPLAY_output.txt').open('w') as f,contextlib.redirect_stdout(f):
    result=runpy.run_path(str(ROOT/'DCU_Mass_Cell_36.py'))['dcu_mass_36']
reference=json.load(gzip.open(ROOT/'DCU_Mass_Cell_36_RESULTS.json.gz','rt'))
counts={'exact_fields':0,'float_fields':0};maxerror=[0.]
def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,f'{path}[{i}]')
    elif isinstance(a,float):
        assert isinstance(b,(float,int)) and math.isfinite(a) and math.isfinite(b),path
        assert math.isclose(a,b,rel_tol=2e-12,abs_tol=2e-12),(path,a,b)
        maxerror[0]=max(maxerror[0],abs(a-b));counts['float_fields']+=1
    else:
        assert type(a)==type(b) and a==b,(path,a,b);counts['exact_fields']+=1
compare(result,reference)
with gzip.open(ROOT/'Cell36_REPLAY_RESULTS.json.gz','wt') as f:json.dump(result,f,separators=(',',':'))
print(json.dumps(dict(status='PASS',counts=counts,max_float_difference=maxerror[0],wall_seconds=time.perf_counter()-start),indent=2))
