#!/usr/bin/env bash
# Zha25 Dataset I replication: one job per complex, NPROC in parallel.
set -u
G=$PROJECT_ROOT/goals/zha25-mmpbsa-replication-20260924
D=$G/work/MM_PBGBSA-CP/dataset1
OUT=${OUT:-$G/work/dataset1_ff03}
FF=${FF:-ff03}
NPROC=${NPROC:-22}
INDI=${INDI:-1,2,4,6,8}
METHODS=${METHODS:-pb,gb1,gb2,gb5}
mkdir -p "$OUT"
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate AmberTools26
# sander links a threaded BLAS and takes ~1.4 cores per process unless pinned
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# head-to-tail set taken from the paper's own Chain type column (BB / BB+SS)
H2T=$(awk -F, 'NR>1 && ($7 ~ /BB/){print $1}' "$G/work/MM_PBGBSA-CP/dataset1.csv" | tr '\n' ' ')
echo "head-to-tail complexes: $H2T"
run_one() {
  P=$1
  if grep -q '"status": "ok"' "$OUT/$P/result.json" 2>/dev/null; then echo "skip $P"; return; fi
  EXTRA=""
  case " $H2T " in *" $P "*) EXTRA="--head-to-tail";; esac
  python "$G/zha_pipeline.py" "$P" "$D/$P/${P}_protein.pdb" "$D/$P/${P}_CP.pdb" \
      "$OUT/$P" --ff "$FF" --indi "$INDI" --methods "$METHODS" $EXTRA \
      > "$OUT/$P.stdout" 2>&1 || echo "FAIL $P"
  echo "done $P"
}
export -f run_one; export G D OUT FF INDI METHODS H2T
ls "$D" | xargs -P "$NPROC" -I{} bash -c 'run_one {}'
echo "=== ALL DONE ==="
