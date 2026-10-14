#!/usr/bin/env bash
# Parameterise every non-canonical residue in the LptA 212-design panel.
set -u
G=$PROJECT_ROOT/goals/zha25-mmpbsa-replication-20260924
W=$G/work/ncaa
OUT=${OUT:-$W/params}
NPROC=${NPROC:-10}
mkdir -p "$OUT"
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate AmberTools26
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
one() {
  C="$1"
  if grep -q '"status": "ok"' "$OUT/$C/result.json" 2>/dev/null; then echo "skip $C"; return; fi
  timeout 3600 python "$W/ncaa_param.py" "$C" "$W/lib" "$OUT/$C" >/dev/null 2>&1 || echo "FAIL $C"
}
export -f one; export W OUT
ls "$W/lib" | xargs -P "$NPROC" -I{} bash -c 'one "{}"'
echo "=== NCAA PARAMETERISATION DONE ==="
