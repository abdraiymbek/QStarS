#!/bin/sh
# DRYRUN batch: grey-seed relaxations at N=400 and N=800, without and with mesh cycles (scratch only)
PY=/Users/abraiymbek/Projects/QStarS/.venv/bin/python
for RUN in C0 C2; do
 for CYC in 0 40; do
  for R in 6000 3464 8485; do
   VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 15 $PY dry_relax.py $RUN $R grey $CYC > log_relax_${RUN}_${R}_grey_cyc${CYC}.txt 2>&1
  done
 done
done
echo BATCH_DONE
