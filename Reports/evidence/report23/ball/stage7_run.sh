#!/bin/bash
# Stage 7 (Report 23) executor: one registered driver command, after the code-freeze check. Usage: stage7_run.sh ARGS...
cd /Users/abraiymbek/Projects/QStarS
if ! shasum -a 256 -c --quiet Reports/evidence/report23/ball/stage7_code_freeze.sha256; then
  echo "CODE FREEZE CHECK FAILED"; exit 99
fi
echo "=== $(date -u '+%Y-%m-%dT%H:%M:%SZ') validate_ball_track.py $*" >> Reports/evidence/report23/ball/stage7_commands.log
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 15 .venv/bin/python scripts/validate_ball_track.py "$@" > /tmp/stage7_last.out 2>&1
rc=$?
echo "    exit $rc; run dir: $(grep -o 'benchmarks/results/[^ ]*-ball-track/[^ "]*' /tmp/stage7_last.out | head -1)" >> Reports/evidence/report23/ball/stage7_commands.log
tail -4 /tmp/stage7_last.out
exit $rc
