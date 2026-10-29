#!/bin/bash
# Report 21 closing chain for a converged T-released Newton-Krylov state (no tuning; fixed gates):
#   1. trace_cleanup.py : converge the trace levels (<= 1e-10 of their element) at the fixed significant state
#      (pairing $PAIRING; structural = CMFGEN's equation-variable pairing, see the script's docstring)
#   2. verify_cool_state.py : one complete native iteration at the cleaned state (spectral export), then
#      audit_cmfgen_cool.py (electron balance, global/cumulative radiation energy, last correction on all
#      levels, populations/charge, spectral material temperature)
# Verification runs only if the cleanup converged (last cleanup_<tag>_history.json record: converged and
# max trace step < its tol); otherwise the chain stops with the history preserved.
# Usage: post_convergence.sh RUN BASE X0 REPLACE [extra VADAT pairs ...]
#   RUN      driver --out directory (it00_base = x0 reference)
#   BASE     converged accepted base pass (directory name inside RUN)
#   X0       archived restart used by the run (--x0)
#   REPLACE  equation table used by the run (--replace)
# Environment (optional): EXE (default build 31, the executable of nk167_hot since its it16 resume; build 28
#   gives the same residuals to <= 3.2e-12 per row at ~2.5x the base-pass time, Report 22), CLEANUP_THREADS
#   (2: base passes are fastest at 2 threads, Report 22), VERIFY_THREADS (6), PAIRING (structural, as declared
#   in benchmarks/cmfgen_cool_validation_plan.json closing_phase; l1 = Report 21 behaviour).
# Every pass uses INC_AD=T and INC_ADV=T, as the Newton-Krylov run did. Destinations are per state and executable:
# RUN/cleanup_<BASE>_<PAIRING>_<BUILD>_NN (BUILD = the build_* directory of EXE, else the executable's name, so a
# rerun with another EXE never reuses cleanup passes of the previous one; trace_cleanup.py also refuses a reused pass
# whose executable_sha256 differs), RUN/verify_<state>[_rN] and its _audit; a rerun never overwrites.
set -euo pipefail
cd "$(dirname "$0")/.."
RUN=$1; BASE=$2; X0=$3; REPLACE=$4; shift 4
EXE=${EXE:-backends/cmfgen/build_31/cur_cmf/exe/cmfgen_dev.exe}
CLEANUP_THREADS=${CLEANUP_THREADS:-2}; VERIFY_THREADS=${VERIFY_THREADS:-6}; PAIRING=${PAIRING:-structural}
VAD=(--vadat INC_AD T --vadat INC_ADV T)
while [ $# -ge 2 ]; do VAD+=(--vadat "$1" "$2"); shift 2; done
PY=.venv/bin/python
BUILD=$(basename "$(dirname "$(dirname "$(dirname "$EXE")")")")
case "$BUILD" in build_*) ;; *) BUILD=$(basename "$EXE" .exe);; esac
TAG="${BASE}_${PAIRING}_${BUILD}"
echo "$(date '+%F %T') post-convergence chain: $RUN/$BASE (exe $EXE, pairing $PAIRING)" | tee -a "$RUN/README.txt"
$PY scripts/trace_cleanup.py "$RUN" "$RUN/$BASE" --x0 "$X0" --replace "$REPLACE" --exe "$EXE" --free-t \
    --threads "$CLEANUP_THREADS" --pairing "$PAIRING" --tag "$TAG" "${VAD[@]}"
# The cleaned state is the base of the last history record (BASE itself when sweep 0 already converged).
STATE=$($PY - "$RUN/cleanup_${TAG}_history.json" <<'EOF'
import json, sys
h = json.load(open(sys.argv[1]))[-1]
if not (h.get('converged') and h['max_trace_step'] < h['tol']):
    sys.exit(f"cleanup NOT converged (sweep {h['sweep']}: max trace step {h['max_trace_step']:.3e} >= tol "
             f"{h['tol']:g}); no verification")
print(h['base'])
EOF
)
NAME=$(basename "$STATE"); DEST="$RUN/verify_$NAME"; n=2
while [ -e "$DEST" ] || [ -e "${DEST}_audit" ]; do DEST="$RUN/verify_${NAME}_r$n"; n=$((n + 1)); done
echo "$(date '+%F %T') verification at $STATE -> $DEST ($VERIFY_THREADS threads)" | tee -a "$RUN/README.txt"
$PY scripts/verify_cool_state.py "$STATE" "$DEST" --x0 "$X0" --replace "$REPLACE" --exe "$EXE" --free-t \
    --threads "$VERIFY_THREADS" "${VAD[@]}"
echo "$(date '+%F %T') audits written to ${DEST}_audit" | tee -a "$RUN/README.txt"
