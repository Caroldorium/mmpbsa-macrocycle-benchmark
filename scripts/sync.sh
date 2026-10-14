#!/usr/bin/env bash
# Pull only the small result files back from Brown.  The local copy is rebuilt
# from scratch each time: arms have been renamed mid-run and a stale local file
# from a superseded protocol must never survive into an analysis.
set -u
R=$GPU_HOST
G=$PROJECT_ROOT/goals/zha25-mmpbsa-replication-20260924
K=$HOME/.ssh/id_ed25519_labgpu
for arm in dataset1_ff03 dataset1_ff03_omm cohort_ff03_omm dataset1_ff03_omm_converged cohort_ff03_omm_converged; do
  rm -rf "results_$arm"; mkdir -p "results_$arm"
  rsync -a -e "ssh -i $K" --include='*/' --include='result.json' --include='FINAL.dat' --exclude='*' \
        "$R:$G/work/$arm/" "results_$arm/" 2>/dev/null
  echo "$arm: $(find results_$arm -name result.json 2>/dev/null | wc -l | tr -d ' ') results"
done
