#!/usr/bin/env bash
# Zha25 MM/PBSA on the 192 alpha-ncAA LptA designs (96 RF-parent, 96 thanatin derivative-7).
set -u
G=$PROJECT_ROOT/goals/zha25-mmpbsa-replication-20260924
IN=$G/work/lpta_ncaa_in
FF=${FF:-ff03}
NCAA=${NCAA:-$G/work/ncaa/params}
ARM=${ARM:-lpta_ncaa_omm}
OUT=$G/work/$ARM
NPROC=${NPROC:-6}
export ZHA_MAXITER=${ZHA_MAXITER:-500}
mkdir -p "$OUT"
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate AmberTools26
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
run_one() {
  IFS=, read -r cid br code cyc pos len cb <<< "$1"
  [ -z "${cid:-}" ] && return
  [ "$cid" = "candidate_id" ] && return
  if grep -q '"status": "ok"' "$OUT/$cid/result.json" 2>/dev/null; then echo "skip $cid"; return; fi
  HT=""
  [ "$cyc" = "head_to_tail_amide" ] && HT="--head-to-tail"
  python "$G/zha_pipeline.py" "$cid" "$IN/$cid/protein.pdb" "$IN/$cid/peptide.pdb" "$OUT/$cid" \
      --ff "$FF" --indi 1,2,4 --methods pb,gb2 --min-engine openmm $HT --ncaa-dir "$NCAA" \
      > "$OUT/$cid.stdout" 2>&1 || echo "FAIL $cid"
  echo "done $cid"
}
export -f run_one; export G IN NCAA OUT ZHA_MAXITER FF
tail -n +2 "$IN/ncaa_manifest.csv" | xargs -P "$NPROC" -I{} bash -c 'run_one "{}"'
echo "=== LPTA NCAA DONE ($ARM, maxiter=$ZHA_MAXITER) ==="
