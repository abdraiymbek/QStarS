#!/bin/bash
# Report 23: restart a running cool_newton_krylov.py run at its next base with new arguments.
# usage: switch_at_base.sh OUT NEXT_K "README note" -- <full driver arguments without --out/--start-base/--start-k>
cd /Users/abraiymbek/Projects/QStarS
O=$1; K=$2; NOTE=$3; shift 4
KK=$(printf '%02d' $K)
until grep -q "it$KK base" $O/stdout.log || grep -q "CONVERGED\|stopping for inspection" $O/stdout.log; do sleep 10; done
grep -q "it$KK base" $O/stdout.log || exit 0
B=$(grep "it$KK base" $O/stdout.log | tail -1 | sed "s/.*it$KK base \([^:]*\):.*/\1/")
pkill -f "cool_newton_krylov.py --out $O "
sleep 3
for p in $(pgrep -f cmfgen_dev.exe); do cwd=$(lsof -a -p $p -d cwd -Fn 2>/dev/null | grep ^n | cut -c2-); case "$cwd" in *$O/*) kill $p;; esac; done
sleep 3
echo "$(date '+%Y-%m-%d %H:%M'): restarted at iterate $K (base $B): $NOTE" >> $O/README.txt
.venv/bin/python scripts/cool_newton_krylov.py --out $O "$@" --start-base $O/$B --start-k $K >> $O/stdout.log 2>&1
echo "$(date '+%Y-%m-%d %H:%M'): driver exited rc=$?" >> $O/README.txt
