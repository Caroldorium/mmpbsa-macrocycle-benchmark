#!/usr/bin/env python3
"""Add back the edges needed to connect each receptor group's perturbation map.

A relative-FEP map yields ddG along edges, so a compound in a disconnected component has no
path to the rest and cannot be placed on the group's dG scale at all.  The quality filter
(delta <= 10 heavy atoms) removes precisely the edges the map builder added to join
components -- they are large because that is why they were needed -- leaving 64% coverage.

Here the good edges are kept as the backbone and the cheapest available edge is added back
per component join, regardless of delta, as long as the MCS converged.  Those edges are
flagged `connector=True`: they are wider perturbations, expected to converge less well, and
should carry less weight when the ddG network is solved.
"""
import argparse, collections, glob, itertools, os
import pandas as pd
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFMCS, rdFingerprintGenerator
RDLogger.DisableLog("rdApp.*")

_MOL = {}
_FP = {}
_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def mol(prep, slug):
    if slug not in _MOL:
        f = f"{prep}/{slug}/ligand_v2000.sdf"
        _MOL[slug] = Chem.MolFromMolFile(f, sanitize=True) or Chem.MolFromMolFile(f, sanitize=False)
    return _MOL[slug]


def fp(prep, slug):
    if slug not in _FP:
        m = mol(prep, slug)
        _FP[slug] = _GEN.GetFingerprint(m) if m is not None else None
    return _FP[slug]


def bridge(prep, ca, cb, timeout=60, k=3):
    """Best cross-component edge, scored on demand.

    The original map only ran MCS on the top fingerprint-ranked pairs, so most pairs spanning
    two components were never scored and cannot be recovered from the map file."""
    cand = []
    for x in ca:
        for y in cb:
            fx, fy = fp(prep, x), fp(prep, y)
            if fx is None or fy is None:
                continue
            cand.append((DataStructs.TanimotoSimilarity(fx, fy), x, y))
    cand.sort(reverse=True)
    best = None
    for _, x, y in cand[:k]:
        A, B = mol(prep, x), mol(prep, y)
        if A is None or B is None:
            continue
        r = rdFMCS.FindMCS([A, B], timeout=timeout, ringMatchesRingOnly=True,
                           completeRingsOnly=False, matchValences=False,
                           atomCompare=rdFMCS.AtomCompare.CompareElements,
                           bondCompare=rdFMCS.BondCompare.CompareAny)
        if r.canceled or r.numAtoms == 0:
            continue
        d = A.GetNumHeavyAtoms() + B.GetNumHeavyAtoms() - 2 * r.numAtoms
        if best is None or d < best[2]:
            best = (x, y, d)
    return best


def components(members, edges):
    adj = collections.defaultdict(set)
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    seen, comps = set(), []
    for s in sorted(members):
        if s in seen:
            continue
        stack, c = [s], set()
        while stack:
            u = stack.pop()
            if u in c:
                continue
            c.add(u)
            seen.add(u)
            stack += [v for v in adj[u] if v not in c]
        comps.append(c)
    return comps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest")
    ap.add_argument("--maps", default="map_*.csv")
    ap.add_argument("--max-delta", type=int, default=10)
    ap.add_argument("--out", default="edge_list_connected.txt")
    ap.add_argument("--prep", required=True)
    ap.add_argument("--bridge-timeout", type=int, default=60)
    a = ap.parse_args()
    man = pd.read_csv(a.manifest)
    keep_rows, added, unreachable = [], 0, 0

    print(f"{'family':12s} {'group':20s} {'n':>4s} {'good':>5s} {'comp':>5s} "
          f"{'added':>6s} {'final':>6s} {'covered':>8s}")
    for f in sorted(glob.glob(a.maps)):
        d = pd.read_csv(f)
        fam = d.family.iloc[0]
        conv = d[~d.mcs_timed_out]
        for g, gg in conv.groupby("receptor"):
            members = set(man[(man.target == fam) & (man.receptor == g)].slug)
            good = gg[gg.delta_heavy <= a.max_delta]
            sel = [(r.a, r.b, int(r.delta_heavy), False) for r in good.itertuples()]
            pool = sorted(((r.a, r.b, int(r.delta_heavy)) for r in gg.itertuples()
                           if r.delta_heavy > a.max_delta), key=lambda t: t[2])
            comps = components(members, [(x, y) for x, y, _, _ in sel])
            n_before = len(comps)
            # join components cheapest-first until no further join is possible
            progress = True
            while len(comps) > 1 and progress:
                progress = False
                idx = {s: i for i, c in enumerate(comps) for s in c}
                for x, y, dlt in pool:
                    if x in idx and y in idx and idx[x] != idx[y]:
                        sel.append((x, y, dlt, True))
                        added += 1
                        comps = components(members, [(p, q) for p, q, _, _ in sel])
                        progress = True
                        break
            # anything still split had no scored candidate pair; score one on demand
            guard = 0
            while len(comps) > 1 and guard < 40:
                guard += 1
                comps.sort(key=len, reverse=True)
                br = bridge(a.prep, comps[0], comps[1], timeout=a.bridge_timeout)
                if br is None:
                    joined = False
                    for j in range(2, len(comps)):
                        br = bridge(a.prep, comps[0], comps[j], timeout=a.bridge_timeout)
                        if br is not None:
                            joined = True
                            break
                    if not joined:
                        break
                sel.append((br[0], br[1], br[2], True))
                added += 1
                comps = components(members, [(p, q) for p, q, _, _ in sel])
            big = max((len(c) for c in comps), default=0)
            unreachable += len(members) - big
            for x, y, dlt, conn in sel:
                keep_rows.append(dict(family=fam, receptor=g, a=x, b=y,
                                      delta_heavy=dlt, connector=conn))
            print(f"{fam:12s} {g:20s} {len(members):4d} {len(good):5d} {n_before:5d} "
                  f"{sum(1 for s in sel if s[3]):6d} {len(comps):6d} {big:4d}/{len(members):<4d}")
    out = pd.DataFrame(keep_rows).drop_duplicates(["a", "b"])
    out.to_csv(a.out.replace(".txt", ".csv"), index=False)
    with open(a.out, "w") as fh:
        for r in out.itertuples():
            fh.write(f"{r.a} {r.b} {r.family.replace('/','_').replace(' ','_')}\n")
    print(f"\nedges: {len(out)}  ({int(out.connector.sum())} added as connectors)")
    print(f"compounds still unreachable: {unreachable}")
    print(f"wrote {a.out} and {a.out.replace('.txt', '.csv')}")


if __name__ == "__main__":
    main()
