#!/usr/bin/env bash
set -u
G=$PROJECT_ROOT/goals/zha25-mmpbsa-replication-20260924
D=$G/work/MM_PBGBSA-CP/dataset1
OUT=${OUT:-$G/work/dataset1_ff03_omm}
NPROC=${NPROC:-3}
INDI=${INDI:-1,2,4,6,8}
METHODS=${METHODS:-pb,pb2,gb1,gb2,gb5}
mkdir -p "$OUT"
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate AmberTools26
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
H2T=$(awk -F, 'NR>1 && ($7 ~ /BB/){print $1}' "$G/work/MM_PBGBSA-CP/dataset1.csv" | tr '\n' ' ')
run_one() {
  P=$1
  if grep -q '"status": "ok"' "$OUT/$P/result.json" 2>/dev/null; then echo "skip $P"; return; fi
  EXTRA=""; case " $H2T " in *" $P "*) EXTRA="--head-to-tail";; esac
  python "$G/zha_pipeline.py" "$P" "$D/$P/${P}_protein.pdb" "$D/$P/${P}_CP.pdb" \
      "$OUT/$P" --ff ff03 --indi "$INDI" --methods "$METHODS" --min-engine openmm $EXTRA \
      > "$OUT/$P.stdout" 2>&1 || echo "FAIL $P"
  echo "done $P"
}
export -f run_one; export G D OUT INDI METHODS H2T
ls "$D" | xargs -P "$NPROC" -I{} bash -c 'run_one {}'
echo "=== OMM ARM DONE ==="
