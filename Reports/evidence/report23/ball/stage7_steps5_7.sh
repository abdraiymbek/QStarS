#!/bin/bash
# Stage 7 steps 5-7 (Report 23): relax every row-0 configuration with Amendments 3, 4, 2, 1; smoke; campaign.
# Stops (and records why) at the first outcome that differs from the runbook's expectation. bash 3.2 compatible
# (the first version used an associative array, which /bin/bash 3.2 rejects; it was stopped during the C1 relax,
# which completed on its own). Row-0 checkpoints are taken from relax directories named on the command line or run here.
# Usage: stage7_steps5_7.sh [C1_RELAX_DIR]
cd /Users/abraiymbek/Projects/QStarS
L=Reports/evidence/report23/ball/stage7_commands.log
RUN=Reports/evidence/report23/ball/stage7_run.sh
AM="--amendment benchmarks/ball2011_track_plan_amendment_3.json --amendment benchmarks/ball2011_track_plan_amendment_4.json --amendment benchmarks/ball2011_track_plan_amendment_2.json --amendment benchmarks/ball2011_track_plan_amendment_1.json"
check() {   # $1 config, $2 relax dir -> prints the row-0 checkpoint path or stops
  st=$(python3 -c "import json;r=json.load(open('$2/relax.json'));print(r['status'], 'screens' if r.get('row0_enforced_screen_errors') else 'ok', r.get('row0_checkpoint'))")
  echo "    step 5 $1: $st ($2)" >> $L
  set -- $st
  if [ "$1" != "converged_unique" ] || [ "$2" != "ok" ]; then echo "STOP: step 5 ended $st" >> $L; exit 1; fi
}
row0() { python3 -c "import json;print('$1/'+json.load(open('$1/relax.json'))['row0_checkpoint'])"; }
D_C1=$1
for r in C1 C0 C2 C5 C6 C7; do
  if [ "$r" = "C1" ] && [ -n "$D_C1" ]; then d=$D_C1; else
    $RUN $AM relax --run $r > /dev/null
    d=$(tail -1 $L | sed -n 's/.*run dir: //p'); fi
  check $r $d
  eval "D_$r=$d"
done
$RUN $AM smoke --row0 $(row0 $D_C1) > /dev/null; rc=$?
d=$(tail -1 $L | sed -n 's/.*run dir: //p')
echo "    step 6 smoke: exit $rc $(python3 -c "import json;r=json.load(open('$d/smoke.json'));print(r.get('status'), r.get('passed'), r.get('accepted_steps'), r.get('median_dt_years'), r.get('rejected_fraction'))")" >> $L
if [ $rc -ne 0 ]; then echo "STOP: step 6 smoke gate failed" >> $L; exit 1; fi
$RUN $AM campaign --row0 C0=$(row0 $D_C0) --row0 C1=$(row0 $D_C1) --row0 C2=$(row0 $D_C2) --row0 C5=$(row0 $D_C5) --row0 C6=$(row0 $D_C6) --row0 C7=$(row0 $D_C7) > /dev/null
echo "    step 7 campaign: exit $? ($(date -u '+%Y-%m-%dT%H:%M:%SZ'))" >> $L
