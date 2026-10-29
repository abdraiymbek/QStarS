Report 21: 167-radius T-released cool closure, INC_AD=T, INC_ADV=T, build 27, T homotopy from lam=0.
Grid: results/cmfgen_ali_lambda167 (nests the 84 grid exactly). Start state: transfer/QSTARS_STEP_167.bin from
results/nk_fixT/it11_trial_l1 (shared radii copied, midpoints log-interpolated; scripts/transfer_state_167.py).
Equation table: transfer/replace_167.txt (identical to results/nk167/replace_167.txt).
it00_base: copy (cp -Rp) of results/nk167/it00_base (T-released 167 restart reference; used only for x0 populations).
2026-10-04 00:04: restarted at iterate 2 (base it02_trial_l1) with --stage-max-iters 2.
2026-10-04 01:27: switched to the residual-shift homotopy (--hom-mode shift, lam restarted at 0 with R0 = RE residual of it03_trial_l1) with --pc-t-floor 1e-3 at iterate 3.
2026-10-04 03:48: iterate 4 product pass jv02 (1 direction, build 27) deadlocked: 0 % CPU for 30 min, all threads waiting in libgomp team barriers inside lin_apply's ray-parallel region nested in qs_jv_freq's inactive direction region (stack sample in Reports/evidence/report21/omp_deadlock_sample.txt). Killed and kept as superseded/it04_jv02_omp_deadlock. Build 28 removes the nesting (threading only).
2026-10-04 03:54: restarted at iterate 4 (base it04_trial_r0.5) with build 28 (threading-only fix; the iterate-8 product pass of the 84 run reproduces build 27 to 1.2e-16, Reports/evidence/report21/build28_vs_build27_product.json). jv00/jv01 of iterate 4 (build 27) are reused as exact products.
2026-10-04 13:20: STOPPED after iterate 14 (lam=1). From iterate 12 each Newton step cooled the front d33-38 (5.5 -> 4.3 -> 3.7 kK, next step a further 14 %), the cold-branch slide of the 84-radius run (results/nk_adv84 iterates 17-21) now on the nested 167 grid: grid refinement does not remove it. Continued as results/nk167_ptc (pseudo-transient continuation on T with split directions) from it12_trial_l1.
