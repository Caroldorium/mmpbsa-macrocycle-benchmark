#!/usr/bin/env bash
# Zha25 MM/PBSA on the 20 canonical RFpeptide LptA binders.
# These are head-to-tail cyclic peptides of canonical amino acids: the Dataset-I
# route (protein force field + explicit cyclisation bond), NOT the GAFF2 ligand route.
set -u
G=$PROJECT_ROOT/goals/zha25-mmpbsa-replication-20260924
IN=$G/work/lpta_canonical_in
FF=${FF:-ff03}
ARM=${ARM:-lpta_canonical_omm}
OUT=$G/work/$ARM
NPROC=${NPROC:-6}
INDI=${INDI:-1,2,4}
METHODS=${METHODS:-pb,gb2}
MIN_ENGINE=${MIN_ENGINE:-openmm}
export ZHA_MAXITER=${ZHA_MAXITER:-500}
mkdir -p "$OUT"
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate AmberTools26
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
run_one() {
  slug="$1"
  [ -z "${slug:-}" ] && return
  if grep -q '"status": "ok"' "$OUT/$slug/result.json" 2>/dev/null; then echo "skip $slug"; return; fi
  python "$G/zha_pipeline.py" "$slug" "$IN/$slug/protein.pdb" "$IN/$slug/peptide.pdb" "$OUT/$slug" \
      --ff "$FF" --indi "$INDI" --methods "$METHODS" --min-engine "$MIN_ENGINE" --head-to-tail \
      > "$OUT/$slug.stdout" 2>&1 || echo "FAIL $slug"
  echo "done $slug"
}
export -f run_one; export G IN OUT INDI METHODS MIN_ENGINE ZHA_MAXITER FF
cat "$IN/manifest.txt" | xargs -P "$NPROC" -I{} bash -c 'run_one "{}"'
echo "=== LPTA CANONICAL DONE ($ARM, maxiter=$ZHA_MAXITER) ==="
