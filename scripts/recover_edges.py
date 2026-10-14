#!/usr/bin/env python3
"""Re-run the MCS for edges whose search timed out, with a much longer budget.

Macrocycle MCS is combinatorially hard and a 6 s budget leaves 30-50% of pairs unresolved in
some families.  Dropping those edges fragments the perturbation map -- MCL-1 fell into 23
components with 15 isolated compounds -- and a disconnected map has no common dG scale, so
those compounds cannot be ranked at all.  Recovering them matters more than the CPU time.
"""
import argparse, os, sys
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS
RDLogger.DisableLog("rdApp.*")

def mol(prep, slug):
    f=os.path.join(prep, slug, "ligand_v2000.sdf")
    m=Chem.MolFromMolFile(f, sanitize=True)
    if m is None: m=Chem.MolFromMolFile(f, sanitize=False)
    return m

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("mapcsv"); ap.add_argument("prep"); ap.add_argument("--timeout", type=int, default=150)
    a=ap.parse_args()
    d=pd.read_csv(a.mapcsv)
    todo=d[d.mcs_timed_out]
    if todo.empty: print(f"{a.mapcsv}: nothing to recover"); return
    print(f"{a.mapcsv}: retrying {len(todo)} timed-out edges at {a.timeout}s")
    cache={}
    def get(s):
        if s not in cache: cache[s]=mol(a.prep,s)
        return cache[s]
    fixed=0
    for i,r in todo.iterrows():
        A,B=get(r.a),get(r.b)
        if A is None or B is None: continue
        res=rdFMCS.FindMCS([A,B], timeout=a.timeout, ringMatchesRingOnly=True,
                           completeRingsOnly=False, matchValences=False,
                           atomCompare=rdFMCS.AtomCompare.CompareElements,
                           bondCompare=rdFMCS.BondCompare.CompareAny)
        delta=r.heavy_a+r.heavy_b-2*res.numAtoms
        d.loc[i,"delta_heavy"]=int(delta)
        d.loc[i,"mcs_timed_out"]=bool(res.canceled)
        if not res.canceled: fixed+=1
        sys.stdout.write(f"\r   recovered {fixed}/{len(todo)}"); sys.stdout.flush()
    d.to_csv(a.mapcsv, index=False)
    ok=d[(d.delta_heavy<=10)&(~d.mcs_timed_out)]
    print(f"\n   {a.mapcsv}: {fixed} converged on retry; usable edges now {len(ok)} of {len(d)}")

if __name__=="__main__": main()
