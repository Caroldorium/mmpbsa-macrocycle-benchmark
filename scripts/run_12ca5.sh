#!/usr/bin/env bash
set -u
G=$PROJECT_ROOT/goals/zha25-mmpbsa-replication-20260924
S=$G/work/12ca5/staged/built
OUT=$G/work/12ca5/scored
NPROC=${NPROC:-8}
mkdir -p "$OUT"
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate AmberTools26
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
run_one() {
  s=$1
  grep -q '"status": "ok"' "$OUT/$s/result.json" 2>/dev/null && { echo "skip $s"; return; }
  q=$(python - "$G/work/12ca5/poses/pose_build_report.csv" "$s" <<'PY'
import csv,sys
sid=sys.argv[2].replace("_",":",1)
for r in csv.DictReader(open(sys.argv[1])):
    if r["system_id"]==sid: print(r.get("formal_charge") or 0); break
else: print(0)
PY
)
  python "$G/zha_cohort.py" "$s" "$S/$s" "$OUT/$s" --charge "${q:-0}" --ff ff03 \
      --indi 1,2,4 --methods pb,gb1,gb2 --min-engine openmm > "$OUT/$s.stdout" 2>&1 || echo "FAIL $s"
  echo "done $s"
}
export -f run_one; export G S OUT
ls "$S" | xargs -P "$NPROC" -I{} bash -c 'run_one {}'
echo "=== 12CA5 DONE ==="
