Report 21: 167-radius T-released cool closure from a second hot initial guess. Initial guess only: the hot-start state
results/nk167_hot/start_base with T raised to >= 10500 K at the net-heated front d30-40 only (RE/c > 1e-4 there),
where the frozen-J local thermal curve at d35 (results/thermal_curve_167/d35_hot, Reports/evidence/report21/
thermal_curve_167_d35_hot.json) has its stable root between 10.27 kK (RE/gross +1.4e-3) and 11.46 kK (-4.7e-3), ~10.5 kK.
Populations and n_e re-equilibrated with J frozen and T held (scripts/equilibrate_at_T.py, prep/), then PTC with split
directions; convergence judged on the unmodified residual.
it00_base: copy (cp -Rp) of results/nk167_adv/it00_base (x0 populations only).
2026-10-04 18:25: preparation started (2 threads, alongside results/nk167_hot).
2026-10-04 19:03: preparation done (7 sweeps, final local frozen-J step 9.9e-5; T = 10.5 kK at d30-40, 7.9 kK around). Solve started: PTC (automatic delta0, --ptc-target-dT 0.2), split, component cap ln 10, 4 threads, alongside results/nk167_hot (8 threads).
2026-10-04 21:19: Report 21 written with this run still in progress (see Reports/Report_21.md 4.5). If it converges: scripts/post_convergence.sh results/nk167_hot2 <base> results/cmfgen_ali_lambda167 results/nk167_adv/transfer/replace_167.txt
