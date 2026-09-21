#!/usr/bin/env python3
"""Compile the checked LaTeX source. Does not regenerate or overwrite source prose."""
import shutil,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    missing=[x for x in ('latexmk','lualatex','biber') if shutil.which(x) is None]
    if missing:raise SystemExit('Required TeX tools missing: '+', '.join(missing)+'. See README.md; the precompiled paper/main.pdf is also supplied.')
    subprocess.run(['latexmk','-lualatex','-interaction=nonstopmode','-halt-on-error','main.tex'],cwd=ROOT/'paper',check=True)
    log=(ROOT/'paper/main.log').read_text(errors='replace')
    flags=[x for x in ('Overfull \\hbox','Overfull \\vbox','undefined references','multiply-defined labels','Missing character:') if x in log]
    if flags:raise SystemExit('Review LaTeX warnings before release: '+', '.join(flags))
    print('Compiled:',ROOT/'paper/main.pdf')
if __name__=='__main__':main()
