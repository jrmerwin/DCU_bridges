"""Run the single cell with the exact frozen Cell29 history file, no old paths."""
from pathlib import Path
from fractions import Fraction
import contextlib,gzip,json,runpy,sys,time
ROOT=Path(__file__).resolve().parent

def clean(x):
    if isinstance(x,Fraction):return str(x)
    if isinstance(x,complex):return {'real':x.real,'imag':x.imag}
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items() if not callable(v)}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if hasattr(x,'tolist'):return clean(x.tolist())
    return x

def main():
    inp=json.load(gzip.open(ROOT/'reference/DCU_Mass_Cell_29_RESULTS.json.gz','rt'))
    outname='REPLAY' if (ROOT/'DCU_Mass_Cell_32_RESULTS.json.gz').exists() else 'DCU_Mass_Cell_32'
    start=time.perf_counter()
    with (ROOT/(outname+'_output.txt')).open('w') as f,contextlib.redirect_stdout(f):
        env=runpy.run_path(str(ROOT/'DCU_Mass_Cell_32.py'),init_globals={'dcu_mass_29':inp})
    data=clean(env['dcu_mass_32'])
    with gzip.open(ROOT/(outname+'_RESULTS.json.gz'),'wt') as f:json.dump(data,f,sort_keys=True,separators=(',',':'))
    if outname=='REPLAY':
        target=json.load(gzip.open(ROOT/'DCU_Mass_Cell_32_RESULTS.json.gz','rt'))
        if data!=target:raise AssertionError('Replay differs from saved result')
    print((ROOT/(outname+'_output.txt')).read_text())
    print(f'Replay complete in {time.perf_counter()-start:.3f} seconds; {outname} files written.')
    return data
if __name__=='__main__':main()
