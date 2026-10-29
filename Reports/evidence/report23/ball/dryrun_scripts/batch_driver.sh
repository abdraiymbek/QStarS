#!/bin/zsh
# DRYRUN batch through the driver functions (scratch output root DRYRUN_driver)
export A3_DRAFT=/private/tmp/claude-501/-Users-abraiymbek-Projects-QStarS/508f230a-38f2-47d3-bc8d-313edee2339c/scratchpad/ball_unblock/amendment_3_for_dry.json
PY=/Users/abraiymbek/Projects/QStarS/.venv/bin/python
for args in "C1 A3 A2 A1" "C0 A3" "C0 A3 A2" "C2 A3" "C2 A3 A2" "C0 A3 A2 A1"; do
  VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 15 $PY dry_driver.py ${=args} 2>&1 | grep "^DRYRUN {\|Error\|Refus" | cut -c1-2500
done
echo BATCH_DONE
