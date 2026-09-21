"""Replay the one-cell calculation. Requires Python3.10+ and NumPy only."""
from pathlib import Path
import contextlib,gzip,io,json,runpy,time,sys
ROOT=Path(__file__).resolve().parent

def replay(save_reference=False):
    stream=io.StringIO();start=time.perf_counter()
    with contextlib.redirect_stdout(stream):
        ns=runpy.run_path(str(ROOT/'DCU_Mass_Cell_37.py'))
    result=ns['dcu_mass_37']
    normalized=json.loads(json.dumps(result))
    name='DCU_Mass_Cell_37' if save_reference else 'REPLAY'
    (ROOT/f'{name}_RESULTS.json.gz').write_bytes(gzip.compress(json.dumps(normalized,sort_keys=True,separators=(',',':')).encode(),mtime=0))
    (ROOT/f'{name}_output.txt').write_text(stream.getvalue())
    print(stream.getvalue(),end='')
    if not save_reference:
        reference=json.load(gzip.open(ROOT/'DCU_Mass_Cell_37_RESULTS.json.gz','rt'))
        def compare(a,b,path='root'):
            if isinstance(a,dict):
                assert set(a)==set(b),path
                for k in a:compare(a[k],b[k],path+'.'+k)
            elif isinstance(a,list):
                assert len(a)==len(b),path
                for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
            elif isinstance(a,float):
                assert abs(a-b)<=2e-12*max(abs(a),abs(b),1),(path,a,b)
            else:assert a==b,(path,a,b)
        compare(reference,normalized)
        print('PASS: reference replay, exact discrete fields and numerical floating tolerance.')
    print(f'Elapsed seconds: {time.perf_counter()-start:.3f}')
    return normalized

if __name__=='__main__':replay('--save-reference' in sys.argv)
