Report 21: 167-radius T-released cool closure from a HOT initial guess (as Report 20 4.6 at 84 radii), build 28,
INC_AD=T, INC_ADV=T. Initial guess only, not a fitted parameter: the lam=1 state results/nk167_adv/it12_trial_l1 with
T raised to >= 7900 K (Report 20's local frozen-J root at 84-radius d12) at 167-depths 19-47 (= 84-depths 10-24), then
populations and n_e at all depths re-equilibrated with J frozen and T held (scripts/equilibrate_at_T.py, prep/).
Reason: at the front the coupled RE balance is net heating ~ proportional to T on the cool branch (Report 21 4.3), so free
Newton runs to T -> 0 and PTC from ~5 kK must cross a thermally unstable stretch; starting above the hot root tests
whether a stable hot root of the unmodified residual exists. The solve then uses PTC (net heating raises T) with split
directions and the unmodified convergence criterion.
it00_base: copy (cp -Rp) of results/nk167_adv/it00_base (x0 populations only).
2026-10-04 14:51: preparation started (4 threads).
2026-10-04 15:20: preparation done (8 sweeps, final local frozen-J step 1.7e-4; T = 7.90 kK plateau at d19-47). Solve started: PTC (automatic delta0 with --ptc-target-dT 0.2), split directions, component cap ln 10, 4 threads.
2026-10-04 17:50: at iterate 2 (base it02_trial_r0.0625) restarted with 8 threads / 3 column workers (warm-start PTC run paused).
2026-10-04 21:19: Report 21 written with this run still in progress (see Reports/Report_21.md 4.5). If it converges: scripts/post_convergence.sh results/nk167_hot <base> results/cmfgen_ali_lambda167 results/nk167_adv/transfer/replace_167.txt
