#!/bin/bash
# Report 27: wait until the v1_p2 driver finishes one more iterate (a new "trial l0.125: merit" line) or exits.
cd /Users/abraiymbek/Projects/QStarS
L=results/nk167_trap/v1_p2/driver.log
n0=$(grep -c "trial l0.125: merit" $L)
while [ "$(grep -c 'trial l0.125: merit' $L)" -le "$n0" ]; do
  results/nk167_trap/driver_running.sh results/nk167_trap/v1_p2 >/dev/null || { echo "DRIVER EXITED"; tail -3 $L; exit 0; }
  sleep 60
done
sleep 30
grep "trial step bounded\|trial l.*merit\|fallback\|stopping" $L | tail -5
