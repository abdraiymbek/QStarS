#!/bin/zsh
# DRYRUN of the Stage 7 command flow through the CLI, every command with --output into scratch (never benchmarks/results).
S=/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad/ball_unblock/final_fix
cd /Users/abraiymbek/Projects/QStarS
O=$S/DRYRUN_cli_flow; mkdir -p $O
A=$S/cli_amend
A3=(--amendment $A/ball2011_track_plan_amendment_3.json)
A34=($A3 --amendment $A/ball2011_track_plan_amendment_4.json)
A342=($A34 --amendment $A/ball2011_track_plan_amendment_2.json)
A3421=($A342 --amendment $A/ball2011_track_plan_amendment_1.json)
run(){ echo "\$ $*"; local t0=$(date +%s); VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 15 .venv/bin/python scripts/validate_ball_track.py "$@" > $O/last.out 2>&1; local rc=$?; tail -c 1500 $O/last.out | grep -v '^ *"' | tail -8; grep -E '"(status|passed|verdict|classification)"' $O/last.out | head -4; echo "  -> exit $rc, wall $(( $(date +%s)-t0 )) s"; }
run $A3 --output $O relax --run C1
run $A34 --output $O relax --run C1
run $A342 --output $O relax --run C1
R=$(ls -d $O/*_relax_C1_n200 | tail -1)
run $A342 --output $O smoke --row0 $R/row0_6000rsun.npz
run $A342 --fallback F1 --output $O smoke --row0 $R/row0_6000rsun.npz
SM=$(ls -d $O/*_smoke_C1_n200 | head -1)
run $A3421 --output $O smoke --row0 $R/row0_6000rsun.npz
run $A342 --output $O smoke --row0 $R/row0_3464rsun.npz
run $A342 --fallback F1 --output $O relax --run C1
run $A342 --output $O compare --run C0=$SM
ls -d benchmarks/results/*-ball-track/* > $O/results_after.txt; diff $S/results_before.txt $O/results_after.txt && echo 'benchmarks/results unchanged'
