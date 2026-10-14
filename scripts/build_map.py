#!/usr/bin/env python3
"""Build a relative-FEP perturbation map for one FEP+ family.

Edges are scored LOMAP-style from the maximum common substructure: the fewer atoms that have
to appear or disappear, the better the edge.  Edges never cross receptor crystal forms -- the
KRAS set is built on two (7YUZ and 8JJS) and a relative free energy between compounds in
different receptor frames is not a well-defined quantity.
"""
import argparse, json, os, sys, itertools
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS, rdFingerprintGenerator
from rdkit import DataStructs
RDLogger.DisableLog("rdApp.*")

def load(prep, slugs):
    out={}
    for s in slugs:
        f=os.path.join(prep,s,"ligand_v2000.sdf")
        if not os.path.exists(f): continue
        m=Chem.MolFromMolFile(f, sanitize=True)
        if m is None: m=Chem.MolFromMolFile(f, sanitize=False)
        if m is not None: out[s]=m
    return out

def fp_rank(mols, members, keep):
    """Morgan-fingerprint pre-screen.  Exhaustive MCS over macrocycles is far too slow --
    514 pairs at up to 20 s each for KRAS alone -- so only the most similar pairs, which are
    the only ones that could give a small perturbation anyway, go to MCS."""
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps={s: gen.GetFingerprint(mols[s]) for s in members}
    pairs=[]
    for x,y in itertools.combinations(members,2):
        pairs.append(((x,y), DataStructs.TanimotoSimilarity(fps[x],fps[y])))
    pairs.sort(key=lambda kv:-kv[1])
    return [p for p,_ in pairs[:keep]]

def mcs_pair(a,b,timeout=6):
    r=rdFMCS.FindMCS([a,b], timeout=timeout, ringMatchesRingOnly=True,
                     completeRingsOnly=False, matchValences=False,
                     atomCompare=rdFMCS.AtomCompare.CompareElements,
                     bondCompare=rdFMCS.BondCompare.CompareAny)
    return r.numAtoms, r.canceled

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("manifest"); ap.add_argument("prep"); ap.add_argument("out")
    ap.add_argument("--family", required=True)
    ap.add_argument("--max-delta", type=int, default=12, help="max heavy atoms changed per edge")
    ap.add_argument("--degree", type=int, default=3, help="target edges per compound")
    a=ap.parse_args()
    import pandas as pd
    man=pd.read_csv(a.manifest)
    man=man[man.target==a.family]
    mols=load(a.prep, list(man.slug))
    print(f"{a.family}: {len(mols)} of {len(man)} ligands read")
    form={r.slug:r.receptor for r in man.itertuples()}
    heavy={s:m.GetNumHeavyAtoms() for s,m in mols.items()}
    groups={}
    for s in mols: groups.setdefault(form.get(s,"?"),[]).append(s)
    edges=[]
    for g,members in sorted(groups.items()):
        print(f"  receptor {g}: {len(members)} compounds -> {len(members)*(len(members)-1)//2} candidate pairs")
        sc={}
        cand=fp_rank(mols, members, keep=min(len(members)*6, len(members)*(len(members)-1)//2))
        print(f"     screening {len(cand)} of {len(members)*(len(members)-1)//2} pairs by MCS")
        for x,y in cand:
            n,cancel=mcs_pair(mols[x],mols[y])
            d=heavy[x]+heavy[y]-2*n
            sc[(x,y)]=(d,cancel)
        ranked=sorted(sc.items(), key=lambda kv: kv[1][0])
        deg={s:0 for s in members}; chosen=[]
        for (x,y),(d,cancel) in ranked:
            if d>a.max_delta: continue
            if deg[x]>=a.degree and deg[y]>=a.degree: continue
            chosen.append((x,y,d,cancel)); deg[x]+=1; deg[y]+=1
        # force connectivity: attach any isolated node by its best available edge
        import collections
        adj=collections.defaultdict(set)
        for x,y,_,_ in chosen: adj[x].add(y); adj[y].add(x)
        seen=set(); comps=[]
        for s in members:
            if s in seen: continue
            stack=[s]; c=set()
            while stack:
                u=stack.pop()
                if u in c: continue
                c.add(u); seen.add(u); stack+=list(adj[u]-c)
            comps.append(c)
        while len(comps)>1:
            best=None
            for i in range(1,len(comps)):
                for x in comps[0]:
                    for y in comps[i]:
                        k=(x,y) if (x,y) in sc else (y,x)
                        if k in sc and (best is None or sc[k][0]<best[2]):
                            best=(x,y,sc[k][0],sc[k][1],i)
            if best is None: break
            chosen.append(best[:4]); comps[0]|=comps[best[4]]; comps.pop(best[4])
        for x,y,d,cancel in chosen:
            edges.append(dict(family=a.family, receptor=g, a=x, b=y,
                              heavy_a=heavy[x], heavy_b=heavy[y], delta_heavy=int(d),
                              mcs_timed_out=bool(cancel)))
        print(f"     -> {len([e for e in edges if e['receptor']==g])} edges, "
              f"components {len(comps)}, isolated {sum(1 for s in members if deg[s]==0)}")
    import pandas as pd
    df=pd.DataFrame(edges); df.to_csv(a.out, index=False)
    print(f"\nTOTAL {len(df)} edges")
    print(df.delta_heavy.describe().round(1).to_string())
    print(f"timed-out MCS: {int(df.mcs_timed_out.sum())}")
    print(f"wrote {a.out}")

if __name__=="__main__": main()
