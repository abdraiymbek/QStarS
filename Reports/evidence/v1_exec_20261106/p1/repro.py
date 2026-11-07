import sys, time
import numpy as np
from p1_common import Problem, ck
for itd in sys.argv[1:]:
    t0 = time.time()
    P = Problem(itd)
    y0 = np.linalg.lstsq(P.A, P.b, rcond=1e-12)[0]
    print(f'{itd}: {P.A.shape[1]} directions ({P.ndense} dense), base merit {P.beta:.6e}; unbounded fit {P.lin(y0):.4f}, '
          f'max|s| {np.abs(P.step(y0)[P.mask]).max():.1f}  ({time.time() - t0:.0f} s)')
    sb, mu, yb = ck.lm_step_bound((P.S, P.A, P.b), P.it.active.astype(float), P.rad, P.mask)
    print(f'   driver LM bound: mu {mu:.3e}, fit {P.lin(yb):.6f}, bound ratio {P.bound_ratio(yb):.4f}  ({time.time() - t0:.0f} s)')
    np.savez(f'repro_{itd}.npz', y0=y0, yb=yb, mu=mu)
