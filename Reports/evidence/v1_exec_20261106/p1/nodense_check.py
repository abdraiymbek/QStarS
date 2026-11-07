"""V1 plan P1 -> P2: the bounded step from the Krylov passes alone (no dense-window columns), per method. Read-only."""
import sys, copy
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from p1_common import cached, ck  # noqa: E402
from box_proj import project_lin  # noqa: E402
d = sys.argv[1]
P = cached(d, dense=False)
it = P.it
y0 = np.linalg.lstsq(P.A, P.b, rcond=1e-12)[0]
_, mu, y_ne = ck.lm_step_bound((P.S, P.A, P.b), it.active.astype(float), P.rad, P.mask)
_, mu_q, y_qr = ck.bounded_step((P.S, P.A, P.b), it.active, P.rad, P.mask, 'qr')
_, _, y_bx = ck.bounded_step((P.S, P.A, P.b), it.active, P.rad, P.mask, 'box')
E = np.stack([(project_lin(it.x, P.Sm[:, i].reshape(P.shape)) - P.Sm[:, i].reshape(P.shape)).ravel() for i in range(P.Sm.shape[1])], 1)
print(f'{d} Krylov only ({P.A.shape[1]} directions): unbounded {P.lin(y0):.4f}; LM-NE {P.lin(y_ne):.4f}; LM-QR {P.lin(y_qr):.4f}; '
      f'BOX {P.lin(y_bx):.4f} (projection inconsistency of the box step {np.abs(E @ y_bx).max():.2e}); '
      f'rank of the consistency map {np.linalg.matrix_rank(E[np.any(E != 0, 1)], 1e-10)}')
