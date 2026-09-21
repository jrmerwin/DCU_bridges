"""Replay Cell31 in a new directory, preserving all supplied outputs. Offline."""
from pathlib import Path
import contextlib,gzip,json,time,math
from support import inputs,clean
ROOT=Path(__file__).resolve().parent

def main():
    env=inputs(); start=time.monotonic()
    with (ROOT/'REPLAY_output.txt').open('w',encoding='utf-8') as out,contextlib.redirect_stdout(out):
        exec(compile((ROOT/'DCU_Mass_Cell_31.py').read_text(),str(ROOT/'DCU_Mass_Cell_31.py'),'exec'),env)
    result=clean(env['dcu_mass_31'])
    with gzip.open(ROOT/'REPLAY_RESULTS.json.gz','wt',encoding='utf-8') as f:
        json.dump(result,f,separators=(',',':'))
    with gzip.open(ROOT/'DCU_Mass_Cell_31_RESULTS.json.gz','rt',encoding='utf-8') as f:
        expected=json.load(f)
    last_bit_differences=[]
    def compare(a,b,path='root'):
        if isinstance(a,float) and isinstance(b,float):
            if a!=b:
                if not math.isclose(a,b,rel_tol=2e-12,abs_tol=1e-12):
                    raise AssertionError(f'Scientific float mismatch at {path}: {a} versus {b}')
                last_bit_differences.append(path)
        elif isinstance(a,dict) and isinstance(b,dict):
            if a.keys()!=b.keys():raise AssertionError(f'Key mismatch at {path}')
            for k in a:compare(a[k],b[k],path+'/'+k)
        elif isinstance(a,list) and isinstance(b,list):
            if len(a)!=len(b):raise AssertionError(f'Length mismatch at {path}')
            for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i))
        elif type(a)!=type(b) or a!=b:
            raise AssertionError(f'Exact-data mismatch at {path}')
    compare(result,expected)
    print(f'PASS: entire Cell31 result reproduced in {time.monotonic()-start:.1f}s.')
    print(f'Within-tolerance floating summation differences: {len(last_bit_differences)}; discrete data exact.')
    print('Saved REPLAY_output.txt and REPLAY_RESULTS.json.gz; supplied outputs unchanged.')

if __name__=='__main__':main()
