#!/usr/bin/env python3
"""Sensitivity of the PB score to the PB radius set.

The released script never emits a tleap stage, so which PBRadii the published
topologies carry is not recoverable from it.  This build sets mbondi2, which is
what igb=2 requires; tleap's own default is bondi.  This rescores structures that
are already minimised, with the radii swapped, to bound the difference.
"""
import argparse, json, os, re, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zha_pipeline import run, parse_final, _cleanup


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir"); ap.add_argument("--radii", default="bondi")
    ap.add_argument("--indi", default="2"); ap.add_argument("--method", default="pb")
    a = ap.parse_args()
    out = os.path.abspath(a.outdir)
    if not os.path.exists(os.path.join(out, "min.nc")):
        print(json.dumps({"dir": os.path.basename(out), "status": "not_built"})); return 1
    log = os.path.join(out, "radii.log"); open(log, "a").close()
    tag = f"radii_{a.radii}"
    d = os.path.join(out, tag); os.makedirs(d, exist_ok=True)
    for stem in ("com", "rec", "lig"):
        open(os.path.join(d, f"{stem}.parmed"), "w").write(
            f"changeRadii {a.radii}\nparmout {stem}.prmtop\ngo\n")
        run(f"parmed -O -i {stem}.parmed ../{stem}.prmtop", d, log)
    body = (f"&general\n  interval=1, verbose=1,\n/\n"
            f"&pb\n  indi={a.indi}, exdi=80.0, inp=1,\n"
            f"  cavity_surften=0.0072, cavity_offset=0.0,\n/\n")
    open(os.path.join(d, "mmpbsa.in"), "w").write("radii sensitivity\n" + body)
    run("MMPBSA.py -O -i mmpbsa.in -cp com.prmtop -rp rec.prmtop -lp lig.prmtop "
        "-y ../min.nc -o FINAL.dat -prefix _MM_", d, log)
    _cleanup(d)
    v = parse_final(os.path.join(d, "FINAL.dat"))
    rj = os.path.join(out, "result.json")
    if v is not None and os.path.exists(rj):
        r = json.load(open(rj)); r[f"pb_{a.indi}_{a.radii}"] = v
        json.dump(r, open(rj, "w"))
    print(json.dumps({"dir": os.path.basename(out), f"pb_{a.indi}_{a.radii}": v}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
