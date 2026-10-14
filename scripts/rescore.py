#!/usr/bin/env python3
"""Add missing method/dielectric scores to an already-built complex directory.

Used to fill in the MM/GBSA arms after the GB interior dielectric route was
fixed, without repeating the minimisations.
"""
import argparse, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zha_pipeline import score, parse_final  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--indi", default="1,2,4,6,8")
    ap.add_argument("--methods", default="gb1,gb2,gb5")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    out = os.path.abspath(a.outdir)
    rj = os.path.join(out, "result.json")
    if not os.path.exists(rj) or not os.path.exists(os.path.join(out, "min.nc")):
        print(json.dumps({"dir": out, "status": "not_built"})); return 1
    res = json.load(open(rj))
    log = os.path.join(out, "rescore.log"); open(log, "a").close()
    changed = False
    for m in a.methods.split(","):
        for indi in a.indi.split(","):
            if not indi:
                continue
            k = f"{m}_{indi}"
            if res.get(k) is not None and not a.force:
                continue
            v = score(out, m, indi, log)
            res[k] = v[0] if isinstance(v, tuple) else v
            changed = True
    if changed:
        json.dump(res, open(rj, "w"))
    print(json.dumps({"dir": os.path.basename(out),
                      "filled": sum(1 for k in res if k.startswith('gb') and res[k] is not None)}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
