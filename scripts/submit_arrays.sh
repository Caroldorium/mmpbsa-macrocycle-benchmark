#!/usr/bin/env bash
# Submit one SLURM array per family.  Walltime is set per family from the measured throughput
# so the cheap families do not queue behind MCL-1's 19 h edges.  Every task resumes: a window
# whose prod.out.gz already exists is skipped, so a job that hits the walltime can simply be
# resubmitted and will continue where it stopped.
set -u
R=$HOME/ZHA25_AMBER_TI
cd "$R"
CONC=${CONC:-20}

# 24 h for every family. The first attempt sized walltime from production throughput alone
# and ignored the per-window minimisation, 100 ps heating and 0.5 ns equilibration: 26 windows
# x ~29 min is ~12.6 h for the smallest family, so 10-18 h produced 35 TIMEOUTs. Tasks resume
# by skipping windows that already have prod.out.gz, so a timeout costs queue time, not work.
declare -A WALL=( [Cyclin_A_B]="24:00:00" [JAK2]="24:00:00" [KRAS]="24:00:00"
                  [MCL-1]="24:00:00" [PCSK9]="24:00:00" )

for FAM in "${!WALL[@]}"; do
  LIST="list_${FAM}.txt"
  [ -s "$LIST" ] || { echo "no edges for $FAM, skipping"; continue; }
  N=$(wc -l < "$LIST")
  cat > "array_${FAM}.sbatch" <<EOF
#!/bin/bash
#SBATCH --account=def-jsbrown_gpu
#SBATCH --job-name=ti_${FAM}
#SBATCH --gpus-per-node=h100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M
#SBATCH --time=${WALL[$FAM]}
#SBATCH --array=1-${N}%${CONC}
#SBATCH --output=logs/${FAM}_%A_%a.out
module --force purge
# amber/22.5's pmemd.cuda is sm_50-sm_80 and cannot run on the H100 nodes
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.6 amber-pmemd/24.3
cd "\$SLURM_SUBMIT_DIR"
EDGE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" ${LIST})
cd "edges/\$EDGE" || exit 1
echo "edge \$EDGE  task \$SLURM_ARRAY_TASK_ID"
for LEG in complex solvent; do
  for L in 0.00000 0.05000 0.10000 0.20000 0.30000 0.40000 0.50000 0.60000 0.70000 0.80000 0.90000 0.95000 1.00000; do
    D=\$LEG/lam_\$L
    [ -f "\$D/prod.out.gz" ] && continue
    # the SPFP GTI kernel overflows on a TI region this large; DPFP is the working path
    pmemd.cuda_DPFP -O -i \$D/min.in   -p \$LEG.prmtop -c \$LEG.inpcrd -o \$D/min.out   -r \$D/min.rst   || exit 2
    pmemd.cuda_DPFP -O -i \$D/heat.in  -p \$LEG.prmtop -c \$D/min.rst  -ref \$D/min.rst -o \$D/heat.out  -r \$D/heat.rst  || exit 3
    pmemd.cuda_DPFP -O -i \$D/equil.in -p \$LEG.prmtop -c \$D/heat.rst -o \$D/equil.out -r \$D/equil.rst || exit 4
    pmemd.cuda_DPFP -O -i \$D/prod.in  -p \$LEG.prmtop -c \$D/equil.rst -o \$D/prod.out -r \$D/prod.rst  || exit 5
    rm -f \$D/min.rst \$D/heat.rst \$D/equil.rst \$D/prod.rst \$D/*.info \$D/min.out \$D/heat.out
    gzip -f \$D/prod.out \$D/equil.out
  done
done
echo "TI EDGE DONE \$EDGE"
EOF
  mkdir -p logs
  JID=$(sbatch "array_${FAM}.sbatch" 2>&1 | grep -oE '[0-9]+$')
  echo "$FAM: $N edges, walltime ${WALL[$FAM]}, array job $JID"
done
