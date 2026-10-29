# Trapped levels: why the Newton–Krylov measure was blind (Report 23 §2.5–2.8)

All numbers are on the 167-radius cool model. Scripts here are read-only analyses of archived passes.

| File | What it shows |
|---|---|
| `n2probe_analysis.py` | Exact H I n = 2 columns (21 depths) against the fixed-J columns: diagonal ratio 3e-4–6.5e-4; least squares in the fixed-J metric cannot reduce the n = 2 row even with weight 1 and the exact columns in the span. |
| `n2local.py` | Structure of the exact n = 2 column: −4.4 on the diagonal at depth 33 (fixed-J: −1.05e4), +4.1 into n = 3 (Hα absorption), 37–54 % of the l1 response at the perturbed depth. |
| `stride_test.py` | Contamination of a coloured probe (every k-th depth perturbed in one product): 38 % at k = 1, 4.6 % at k = 4, 2.5 % at k = 6 (median, n = 2 column). |
| `mprime_test.py` | Preconditioner with the exact n = 2 column: max over significant variables of M⁻¹R rises from 7.7e-3 to 1.63 (H II, n_e at the front) at the state the old measure called nearly converged. |
| `mali_analysis.py` | One pass with the native local-ALI linearization (`QS_LIN_FAC 1`, 33 min): its n = 2 column equals the exact local column to 0.1–5 %; it lists every level whose fixed-J diagonal is more than 2× the trapped one (H I n = 1–9, He I ground and ¹P series, He II n = 1–3). As a full preconditioner its blocks are nearly singular (M⁻¹R up to 40 with the wrong sign), as Report 20 found. |
| `closure_test.py` | The coupling-significance rule tried and withdrawn: it promotes ultra-trace He II/He III levels through ionization-equation rows. |
| `it00_probe_analysis.py` | The cool start state with exact n = 1 and n = 2 columns: H I ground correction up to 366 (depth 18) where the fixed-J measure gave 1.7. |
| `b1_profile.py`, `departure_coefficients.txt` | Departure coefficients of H I n = 1 and n = 2. Cool start: b₁ = 1.0–1.3 inside depth 43 and n₂/n₁ Boltzmann to 1e-3 (memory of the LTE start of the Λ iteration). PTC states (`nk167_hot` it34, pin state): b₁ up to 6.9e10 and n₂/n₁ down to 5.7e-9 of Boltzmann at depths 20–47: T was raised by 3–6 kK while the populations stayed at their cool values. |
