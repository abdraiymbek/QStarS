Passes run with --free-t (VADAT FIX_T=F) BEFORE the T-release fix (2026-10-01 23:45).
On a restart CMFGEN reads FIXED_T back from BAMATPNT ("Temperature was held fixed" in the
restart directories results/cmfgen_ali_lambda84_long_maxlam1e3 and results/cmfgen_ali_lambda167),
so FIXPOP_IN_BA_V3 zeroed the radiative-equilibrium row at all depths ("Temperature held fixed at
all depths." in OUTGEN; residual row 134 identically zero; BA T row = identity). These are
therefore fixed-T passes despite FIX_T=F. Kept for provenance; not used.
  nk_freeT/            it00_base and start_base (state it10_trial_l1 of results/nk_fixT)
  nk167_it00_base/     167-radius restart base pass (populations at the restart are still valid)
  nk_dfreq_it00_base/  doubled-frequency base pass at x0
The same inheritance applies to results/jv/oracle_t and results/jv/b21baseT: their T direction
validated the T COLUMN (other rows' response to T), not the radiative-equilibrium row.
Fixed in scripts/cool_newton_krylov.py:release_t (removes BAMATPNT/BAMAT when T is free).
