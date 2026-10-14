#!/usr/bin/env bash
set -u
G=$PROJECT_ROOT/goals/zha25-mmpbsa-replication-20260924
P=$PROJECT_ROOT/goals/universal-de-novo-discovery-oracle/work/rosetta_portable/prepared
OUT=${OUT:-$G/work/cohort_ff03}
NPROC=${NPROC:-10}
INDI=${INDI:-1,2,4,6,8}
METHODS=${METHODS:-pb,gb2}
MIN_ENGINE=${MIN_ENGINE:-sander}
mkdir -p "$OUT"
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate AmberTools26
# sander links a threaded BLAS and takes ~1.4 cores per process unless pinned
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
run_one() {
  IFS=, read -r slug target sid hv q dg fep <<< "$1"
  [ -z "${slug:-}" ] && return
  if grep -q '"status": "ok"' "$OUT/$slug/result.json" 2>/dev/null; then echo "skip $slug"; return; fi
  python "$G/zha_cohort.py" "$slug" "$P/$slug" "$OUT/$slug" --charge "${q:-0}" \
      --ff ff03 --indi "$INDI" --methods "$METHODS" --min-engine "$MIN_ENGINE" \
      > "$OUT/$slug.stdout" 2>&1 \
      || echo "FAIL $slug"
  echo "done $slug"
}
export -f run_one; export G P OUT INDI METHODS MIN_ENGINE
tail -n +2 "$G/cohort_manifest.csv" | xargs -P "$NPROC" -I{} bash -c 'run_one "{}"'
echo "=== COHORT DONE ==="
