DRYRUN scratch of the plan-compliance review of the Stage 7 unblock fix round (2026-10-13 local / 2026-10-14 UTC).
Nothing here is a registered run; nothing was written under benchmarks/results. No row-0 value was compared with the published tables.
rp_relax.py            driver relax_seeds with a scratch root: registered (C1, C0, C2), A3, A3+A2 (C1)   -> DRYRUN_out/, log_registered_C0_C2.txt, log_relax_C1_A3*.txt
rp_cmp_registered.py   dry-run registered relax vs the preserved 2026-10-14 runs                           -> log_registered_compare.txt
rp_oldsolve.py         grey seeds and one step with the current and the preserved photosphere.solve        -> log_oldsolve.txt
rp_a3_oldsolve.py      Amendment 3 settings with the preserved photosphere.solve (C1, C0)                  -> DRYRUN_oldsolve/, log_a3_oldsolve.txt
rp_solve_grid.py       old vs new solve vs brute force on a 194-point (L, R) lattice (no model built)      -> log_solve_grid.txt
rp_attack_d.py, rp_reasons.py  row-0 uniqueness/signature attacks                                          -> DRYRUN_attack/, log_attack_d.txt, log_reasons.txt
rp_attack_e.py         start-row screens, injected failures, labels                                        -> DRYRUN_attack_e/, log_attack_e.txt
rp_smoke.py            registered smoke gate with and without F1                                           -> DRYRUN_smoke/, log_smoke.txt
rp_f1.py               F1 on one MESA step (forced Newton failure)                                         -> log_f1.txt
rp_verdict.py          synthetic verdict with an enforced flag at row 0                                    -> log_verdict.txt
before_copy/           new tests against the pre-change source snapshot                                    -> log_before_check.txt
targeted_nomesa.log, targeted_mesa.log, full_suite.log   pytest
