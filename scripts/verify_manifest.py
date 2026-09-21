#!/usr/bin/env python3
"""Verify distributed bytes without modifying files or accessing the network."""
from pathlib import Path
import hashlib,sys
ROOT=Path(__file__).resolve().parents[1]
def main():
    manifest=ROOT/'SHA256SUMS.txt'
    if not manifest.exists():raise SystemExit('Missing SHA256SUMS.txt')
    checked=0;errors=[]
    for line in manifest.read_text().splitlines():
        digest,name=line.split('  ',1);p=ROOT/name
        if not p.is_file():errors.append('Missing: '+name);continue
        if hashlib.sha256(p.read_bytes()).hexdigest()!=digest:errors.append('Changed: '+name)
        checked+=1
    if errors:raise SystemExit('\n'.join(errors))
    print(f'PASS: {checked} distributed files match SHA256SUMS.txt.')
if __name__=='__main__':main()
