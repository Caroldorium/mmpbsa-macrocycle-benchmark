#!/usr/bin/env bash
# Build every usable TI edge across all five FEP+ families.
#
# "Usable" means the MCS converged and at most 10 heavy atoms change.  Edges that fail to
# build are recorded rather than silently dropped -- a missing edge can disconnect the
# perturbation map, and a disconnected map has no common dG scale for the compounds in the
# orphaned component.
set -u
G=$PROJECT_ROOT/goals/zha25-mmpbsa-replication-20260924
W=$G/work/ti
PREP=$PROJECT_ROOT/goals/universal-de-novo-discovery-oracle/work/rosetta_portable/prepared
COHORT=$G/work/cohort_ff03_omm
NPROC=${NPROC:-10}
cd "$W"
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate AmberTools26

python3 - <<'PY' > edge_list.txt
import pandas as pd, glob
rows = []
for f in sorted(glob.glob("map_*.csv")):
    d = pd.read_csv(f)
    ok = d[(d.delta_heavy <= 10) & (~d.mcs_timed_out)]
    for r in ok.itertuples():
        rows.append(f"{r.a} {r.b} {r.family.replace('/','_').replace(' ','_')}")
print("\n".join(rows))
PY
N=$(wc -l < edge_list.txt)
echo "edges to build: $N"

build_one() {
  A=$1; B=$2; FAM=$3
  D="edges/${A}__${B}"
  [ -f "$D/complex.prmtop" ] && [ -f "$D/solvent.prmtop" ] && { echo "skip $A $B"; return; }
  rm -rf "$D"
  if ! python setup_edge.py "$A" "$B" "$D" --cohort "$COHORT" --prep "$PREP" >"$D.setup.log" 2>&1; then
    echo "SETUP_FAIL $A $B $FAM" >> build_failures.txt; return
  fi
  if ! python build_ti.py "$D" --prod-ns 5 >"$D.build.log" 2>&1; then
    echo "BUILD_FAIL $A $B $FAM" >> build_failures.txt; return
  fi
  echo "ok $A $B"
}
export -f build_one; export W PREP COHORT

: > build_failures.txt
mkdir -p edges
xargs -P "$NPROC" -n3 bash -c 'build_one "$@"' _ < edge_list.txt
echo "=== BUILD COMPLETE ==="
echo "built: $(ls -d edges/*/ 2>/dev/null | wc -l)"
echo "with both topologies: $(ls edges/*/complex.prmtop 2>/dev/null | wc -l)"
echo "failures: $(wc -l < build_failures.txt)"
head -20 build_failures.txt
