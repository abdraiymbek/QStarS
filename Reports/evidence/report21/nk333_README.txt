Report 21: 333-radius T-released cool closure for the spatial resolution check (167 vs 333), INC_AD=T, INC_ADV=T, build 28.
Grid: results/cmfgen_ali_lambda333 (scripts/prepare_mesh333.py: prepare_cool factor 4 + one native bootstrap iteration,
numerical controls VADAT/IN_ITS copied from results/cmfgen_ali_lambda167; radius and velocity identical to the 167 and
84 grids at shared points, mesh333_control.json). Equation table transfer/replace_333.txt mapped from the 167 table
(scripts/transfer_state_nested.py --table-only; midpoints take their outer shared neighbour's selection).
it00_base: unperturbed T-released base pass at the 333 restart (make_base.py), used only for x0 populations.
2026-10-04 11:31: base pass started (4 threads, alongside the 167 solve).
2026-10-04 11:51: it00_base completed (1162 s, 4 threads, 407 MB after compaction).
