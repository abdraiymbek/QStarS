Report 21: 167-radius T-released cool closure with pseudo-transient continuation on T (--hom-mode ptc) and split
Krylov directions (--split), build 28, INC_AD=T, INC_ADV=T. Same grid/table as results/nk167_adv.
Reason: free Newton at lam = 1 slid down the cool branch (results/nk167_adv iterates 12-14), as at 84 radii; PTC follows
the physical relaxation direction (net heating raises T) and stops only at a root of the unmodified residual; split
directions let the least squares size the T step (84-radius test: linear residual 0.11 instead of 0.6-0.9).
Start: results/nk167_adv/it12_trial_l1 (first lam = 1 state after one Newton step; front plateau 5.4-5.5 kK).
it00_base: copy (cp -Rp) of results/nk167_adv/it00_base (x0 populations only).
2026-10-04 13:20: started (6 threads, 3 column workers, alongside the 84-radius PTC test).
2026-10-04 13:31: restarted before the first product with --component-log-cap ln(10) (see results/nk84_ptc/README.txt); completed column passes reused.
2026-10-04 14:59: at iterate 1 (base it01_trial_l1): SER decrease-only and delta reset to 0.03 (pseudo-steps ~4x larger), as for the 84-radius PTC run.
2026-10-04 15:04: the delta-reset restart stopped at once: the stopped driver's iterate-1 product pass jv00 was created after the switch script's check (kept as superseded/it01_jv00_switch_race). Relaunched.
2026-10-04 15:05: (correction) the moves above failed because superseded/ did not exist yet; the partial iterate-1 product pass is now superseded/it01_jv00_switch_race. Relaunched with --ptc-delta-reset 0.03.
2026-10-04 17:01: PAUSED after iterate 1 (accepted full step, natural 1.16 -> 0.18; front d33-38 heated 2.5-3.8 % per iterate at delta 0.03; RE net/gross 1.14e-2 -> 1.18e-2 while crossing the unstable stretch). The machine is given to the hot-start run (results/nk167_hot), which starts closer to the stable root; this run can be resumed from it02_trial_l1.
