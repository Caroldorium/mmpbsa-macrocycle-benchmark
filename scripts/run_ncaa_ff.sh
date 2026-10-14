#!/usr/bin/env bash
# Build one prep set per protein force field: ff03 types the alpha carbon CT, ff14SB/ff19SB CX.
set -u
G=$PROJECT_ROOT/goals/zha25-mmpbsa-replication-20260924
W=$G/work/ncaa
FF=${FF:-ff14SB}
OUT=${OUT:-$W/params_cx}
NPROC=${NPROC:-10}
mkdir -p "$OUT"
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate AmberTools26
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
one() {
  C="$1"
  grep -q '"status": "ok' "$OUT/$C/result.json" 2>/dev/null && { echo "skip $C"; return; }
  timeout 3600 python "$W/ncaa_param.py" "$C" "$W/lib" "$OUT/$C" --protein-ff "$FF" >/dev/null 2>&1 || echo "FAIL $C"
}
export -f one; export W OUT FF
ls "$W/lib" | xargs -P "$NPROC" -I{} bash -c 'one "{}"'
# the derivative-7 frame residues: DAB and 4FO carry a side-chain NH3+, and LE1's SG is
# cross-linked so its thiol hydrogen exists only for the charge fit
for c in EU0 HYP; do timeout 3600 python "$W/ncaa_param.py" $c "$W/framelib" "$OUT/$c" --protein-ff "$FF" >/dev/null 2>&1; done
for c in DAB 4FO; do timeout 3600 python "$W/ncaa_param.py" $c "$W/framelib" "$OUT/$c" --protein-ff "$FF" --net-charge 1 >/dev/null 2>&1; done
timeout 3600 python "$W/ncaa_param.py" LE1 "$W/framelib" "$OUT/LE1" --protein-ff "$FF" --cap-h SG >/dev/null 2>&1
python3 - "$OUT" <<'PY'
import json,glob,os,sys,collections
out=sys.argv[1]
def nn(s): return s[1:]+s[0] if s and s[0].isdigit() else s
n=0
for d in sorted(glob.glob(f'{out}/*/')):
    code=os.path.basename(d.rstrip('/'))
    src=next((f'{os.path.dirname(out)}/{b}/{code}/{code}.params' for b in ('lib','framelib')
              if os.path.exists(f'{os.path.dirname(out)}/{b}/{code}/{code}.params')), None)
    if not src: continue
    names=[nn(l.split()[1]) for l in open(src,errors='ignore') if l.startswith('ATOM ')]
    json.dump(names, open(os.path.join(d,'atom_names.json'),'w')); n+=1
st=collections.Counter(); miss=0
for f in glob.glob(f'{out}/*/result.json'):
    d=json.load(open(f)); st[d.get('status')]+=1; miss+=d.get('tleap_missing_param',0) or 0
print(f'{out}: {sum(st.values())} residues  status={dict(st)}  total_missing_params={miss}  name_maps={n}')
PY
echo "=== NCAA $FF DONE ==="
