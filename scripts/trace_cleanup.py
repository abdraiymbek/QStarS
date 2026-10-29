#!/usr/bin/env python3
"""Closing phase of the cool-closure solve: converge the trace levels.

The Newton-Krylov criterion covers levels holding more than 1e-10 of their element. The
native audit's last-correction gate covers ALL levels. Trace levels (<= 1e-10 of their
element) do not feed back on n_e, T or any significant level and negligibly on J (largest:
H- bound-free, <= 3.2e-4 of the electron-scattering opacity at nk167_hot it24), so with the
significant state, n_e, T and J held fixed their equations are LINEAR in the trace populations
and x(1 - s) solves them exactly (fallback: capped exp(-s) for s >= 1). Each sweep is one
complete native base pass (formal-ALI warm-started from the previous J). Significant levels are
never changed here. To first order no trace update changes the Newton (or native) correction of
the significant variables (J^-1 J[:,T] s_T has no significant component).

Which rows determine the trace variables (--pairing):
  l1 (default; Report 21): each trace variable with its own row, unless that row's l1 weight on
     significant columns is >= that on the trace columns of its depth (row and variable are then
     held fixed: 'demoted'). Each row is judged against the ORIGINAL trace set of its depth (an
     earlier in-loop numpy view made the result depend on row order: 8844 vs 8905 kept at it24).
     The test is a proxy for CMFGEN's conservation rows but also demotes ordinary SE rows of trace
     levels fed mainly by significant levels (H-, H I excited levels, He II high levels), which are
     then never converged. Per-depth diagonal blocks (no advection coupling).
  structural (opt-in): the pairing of CMFGEN's assembly (generate_full_matrix_v3.f; see
     native_correction_check.structural_pairing): own SE rows; the top ion / a replaced ground
     with the ionization-equation row C_ION linking it to a significant (or already paired)
     adjacent ion, read from the pass's QSTARS_REPLACE_TRACE; conservation rows never used;
     variables with no row of their own in the assembled matrix are held fixed and reported
     ('unpaired'). Solved as ONE sparse system over all depths with the native lo/up blocks
     (native_correction_check.structural_step; the advection coupling of He III / He II rows is not
     negligible: diagonal-only solutions differ by up to 70 % at it24), row/column equilibrated;
     componentwise backward error recorded (2e-13 at it24).
  Either way the result is the trace equilibrium AT the fixed significant state. The audit's
  native correction of the trace levels afterwards is the slaved response to the significant
  native correction, which at the cool front is amplified up to ~4e7 (He III, it24): measure it
  with scripts/native_correction_check.py (closing_phase in benchmarks/cmfgen_cool_validation_plan.json).
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import cool_newton_krylov as ck  # noqa: E402
import native_correction_check as ncc  # noqa: E402
from qstars.cmfgen_backend import sha256  # noqa: E402
from qstars.nlte_newton import read_ba_blocks, read_trace_record  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run', type=Path, help='Newton-Krylov output directory (it00_base = x0)')
    ap.add_argument('start', type=Path, help='converged base pass (directory)')
    ap.add_argument('--x0', type=Path, default=ROOT / 'results/cmfgen_ali_lambda84_long_maxlam1e3')
    ap.add_argument('--exe', type=Path, default=ROOT / 'backends/cmfgen/build_23/cur_cmf/exe/cmfgen_dev.exe')
    ap.add_argument('--replace', type=Path, default=ROOT / 'results/jv/replace_fac1_x0.txt')
    ap.add_argument('--free-t', action='store_true')
    ap.add_argument('--vadat', nargs=2, action='append', default=[], metavar=('KEY', 'VALUE'))
    ap.add_argument('--sig-frac', type=float, default=1e-10)
    ap.add_argument('--tol', type=float, default=1e-4)
    ap.add_argument('--max-sweeps', type=int, default=20)
    ap.add_argument('--cap', type=float, default=np.log(1e4))
    ap.add_argument('--threads', type=int, default=6)
    ap.add_argument('--pairing', choices=['l1', 'structural'], default='l1',
                    help='rows that determine the trace variables (module docstring); structural is opt-in')
    ap.add_argument('--tag', default=None,
                    help='names: RUN/cleanup_<tag>_NN and RUN/cleanup_<tag>_history.json (default <start>_<pairing>), '
                         'so cleanups of different states never collide')
    a = ap.parse_args()
    run = a.run.resolve()
    tag = a.tag or f'{a.start.resolve().name}_{a.pairing}'
    ck.check_one_record(run / 'it00_base' / 'QSTARS_RESID_TRACE.bin')       # stale restart record (nk333 it00_base)
    x0, _ = read_trace_record(run / 'it00_base' / 'QSTARS_RESID_TRACE.bin')
    nt, nd = x0.shape
    R = ck.Runner(a.x0, a.exe, a.replace, not a.free_t, a.threads, run / 'cleanup.log', [tuple(v) for v in a.vadat])
    base = a.start.resolve()
    hist = []
    for k in range(a.max_sweeps + 1):
        ncc.one_record(base / 'QSTARS_RESID_TRACE.bin', nt, nd)
        x, Rv = read_trace_record(base / 'QSTARS_RESID_TRACE.bin')
        lo, di, up = read_ba_blocks(base / 'QSTARS_BA_TRACE.bin', nt, nd)
        frac = ncc.fractions(x)
        active = Rv != 0.0
        rp = base / 'QSTARS_REPLACE_TRACE'
        rep = ncc.read_replace(rp if rp.exists() else a.replace, nd)
        if not np.array_equal(rep, ncc.read_replace(a.replace, nd)):
            raise SystemExit(f'{rp}: equation selection differs from the forced table {a.replace}')
        rec = dict(sweep=k, base=str(base), pairing=a.pairing, tol=a.tol, converged=False)
        if a.pairing == 'structural':
            s, solved, held, bw = ncc.structural_step(x, Rv, lo, di, up, rep, a.sig_frac)
            rec.update(trace_variables=int(solved.sum()), held_fixed_unpaired=int(held.sum()), backward_error=bw)
        else:
            solved, held = ncc.l1_split(frac, active, di, a.sig_frac)
            s = np.zeros_like(x)
            for d in range(nd):
                t = np.nonzero(solved[:, d])[0]
                if t.size:
                    s[t, d] = np.linalg.solve(di[d][np.ix_(t, t)], Rv[t, d])
            rec.update(trace_variables=int(solved.sum()), held_fixed_demoted=int(held.sum()))
        smax = float(np.max(np.abs(s))) if solved.any() else 0.0
        rec['max_trace_step'] = smax
        # native correction of the next complete iteration at this state, per class (reported only; the
        # significant part is unchanged by this phase to first order)
        cls, _ = ncc.check(base, a.sig_frac, a.replace)
        rec['native_correction'] = {n: v.get('max_abs', 0.0) for n, v in cls['classes'].items()}
        rec['native_correction_audit_normalized'] = {n: v.get('audit_normalized', 0.0) for n, v in cls['classes'].items()}
        hist.append(rec)
        nc = rec['native_correction']
        R.log(f'cleanup {tag} sweep {k}: {rec["trace_variables"]} trace variables ({a.pairing}), max|s_trace|={smax:.3e}, '
              f'native correction significant {nc["significant"]:.3e} / T {nc["T"]:.3e} / all {nc["all"]:.3e}')
        if smax < a.tol:
            rec['converged'] = True
        (run / f'cleanup_{tag}_history.json').write_text(json.dumps(hist, indent=1))
        if smax < a.tol or k == a.max_sweeps:
            break
        # With the significant state, J, n_e and T fixed, the trace equations are LINEAR in the trace
        # populations, so x(1 - s) is their exact solution; only a non-positive target (s >= 1)
        # falls back to the capped multiplicative decrease.
        f = 1.0 - s
        xn = np.where(f > 1e-4, x * f, x * np.exp(-np.minimum(s, a.cap)))
        dest = run / f'cleanup_{tag}_{k + 1:02d}'
        man = dest / 'run_manifest.json'
        if man.exists() and json.loads(man.read_text()).get('executable_sha256') != sha256(a.exe):
            # Runner.run reuses a pass by state/replay/directions only (the 167 restart relies on that)
            raise SystemExit(f'{dest} was made by another executable than {a.exe}; refusing reuse (use another --tag)')
        base = R.run(dest, 1.0 - xn / x0, warm=base / 'QSTARS_RJ.bin')
    R.log(f'cleanup {tag} finished at {base.name} ({"converged" if hist[-1]["converged"] else "NOT converged"})')


if __name__ == '__main__':
    main()
