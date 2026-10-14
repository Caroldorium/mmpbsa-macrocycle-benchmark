#!/usr/bin/env python3
"""Compare the rebuilt topologies against the authors' deposited minimisation restarts.

The repository ships `minimized/<PDB>/minimization.inpcrd` but no topology, so the
atom count is the only thing that can be checked directly -- and it is the check
that matters, because every way this build can go wrong (a missed cyclisation
bond, a fabricated chain break, a terminal variant that should not be there)
changes it.
"""
import glob, json, io, os, struct, sys
import pandas as pd

RES = sys.argv[1] if len(sys.argv) > 1 else "results_dataset1_ff03_omm"
REPO = sys.argv[2] if len(sys.argv) > 2 else "../MM_PBGBSA-CP"


def netcdf_natom(path):
    b = open(path, "rb").read(4096)
    if b[:3] != b"CDF":
        return None
    f = io.BytesIO(b); f.read(4); struct.unpack(">i", f.read(4))
    if struct.unpack(">i", f.read(4))[0] != 10:
        return None
    dims = {}
    for _ in range(struct.unpack(">i", f.read(4))[0]):
        ln = struct.unpack(">i", f.read(4))[0]
        nm = f.read(ln).decode(); f.read((4 - ln % 4) % 4)
        dims[nm] = struct.unpack(">i", f.read(4))[0]
    return dims.get("atom")


rows = []
for f in sorted(glob.glob(os.path.join(RES, "*", "result.json"))):
    r = json.load(open(f))
    ic = os.path.join(REPO, "minimized", r["pdb"], "minimization.inpcrd")
    rows.append(dict(pdb=r["pdb"], status=r.get("status"), built=r.get("natom"),
                     authors=netcdf_natom(ic) if os.path.exists(ic) else None,
                     ter_removed=(r.get("spurious_ter_removed") or {}).get("rec"),
                     cyx_fallback=r.get("term_cyx_fallback")))
d = pd.DataFrame(rows)
if d.empty:
    print("no results yet"); raise SystemExit(0)
d["diff"] = d.built - d.authors
print(f"complexes: {len(d)}   exact atom-count match with the authors: "
      f"{int((d['diff'] == 0).sum())}")
print(d.sort_values("diff", key=lambda s: s.abs(), ascending=False).head(15).to_string(index=False))
d.to_csv("zha25_build_check.csv", index=False)
print("\nspurious chain breaks removed: total %s across %s receptors"
      % (int(d.ter_removed.fillna(0).sum()), int((d.ter_removed.fillna(0) > 0).sum())))
