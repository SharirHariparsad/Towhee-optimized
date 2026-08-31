#!/usr/bin/env python3
"""
tolerance_diff.py - Phase 4 tolerance-based regression check (1e-6 relative)

Usage:
  python3 tolerance_diff.py answer_current answer_new [tol]
  exit 0 if within tolerance, 1 if exceeds tolerance, prints diff summary.
"""
import sys, re
from pathlib import Path

DEFAULT_TOL = 1e-6

def is_num(tok):
    try:
        float(tok)
        return True
    except: return False

def compare(ref, test, tol):
    rl = Path(ref).read_text(errors='ignore').splitlines()
    tl = Path(test).read_text(errors='ignore').splitlines()
    diffs=[]
    if len(rl)!=len(tl):
        diffs.append(f"line count {len(rl)} vs {len(tl)}")
    for i,(r,t) in enumerate(zip(rl,tl),1):
        if r==t: continue
        if r.split()==t.split(): continue
        rt=r.split(); tt=t.split()
        if len(rt)!=len(tt):
            diffs.append(f"line {i}: token count {len(rt)} vs {len(tt)}")
            diffs.append(f"  ref: {r[:200]}")
            diffs.append(f"  new: {t[:200]}")
            continue
        ok=True
        for a,b in zip(rt,tt):
            if is_num(a) and is_num(b):
                av,bv=float(a),float(b)
                if abs(av-bv) > tol*max(abs(av),1.0):
                    ok=False; break
            else:
                if a!=b: ok=False; break
        if not ok:
            diffs.append(f"line {i} exceeds tol {tol}:")
            diffs.append(f"  ref: {r[:200]}")
            diffs.append(f"  new: {t[:200]}")
            if len(diffs)>100:
                diffs.append("... truncated")
                break
    return diffs

if __name__=="__main__":
    if len(sys.argv)<3:
        print(f"Usage: {sys.argv[0]} ref test [tol]")
        sys.exit(2)
    ref, test = sys.argv[1], sys.argv[2]
    tol = float(sys.argv[3]) if len(sys.argv)>3 else DEFAULT_TOL
    diffs=compare(ref,test,tol)
    if not diffs:
        print(f"tolerance_diff: PASS within tol={tol}")
        sys.exit(0)
    else:
        print(f"tolerance_diff: FAIL tol={tol}, {len(diffs)} diff lines")
        for d in diffs[:40]:
            print(d)
        sys.exit(1)
