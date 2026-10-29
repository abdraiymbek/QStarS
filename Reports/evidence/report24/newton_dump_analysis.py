#!/usr/bin/env python3
"""Report 24: existence of a nonnegative scattering fixed point at the frequencies where build 36's Newton fallback ran.

    newton_dump_analysis.py RUN_DIR [OUT.json]
Reads QSTARS_ALI_NEWTON_DUMP.bin (QSTARS_ALI_NEWTON_DUMP=1): per fallback frequency, J (Newton start), Phi(J), the
assembled A = dPhi/dJ (one-sided differences of the formal solution), CHI, CHI_SCAT, ETA. Phi is the formal solution
of TA = ETA + CHI_SCAT*J and is affine in J up to the interpolation of the emissivity, so its fixed point is
J* = J + (I - A)^-1 (Phi(J) - J). A is nonnegative; by Perron-Frobenius a nonnegative J* for a positive source exists
iff the spectral radius rho(A) < 1. Reports rho(A), the depths carrying the Perron vector, min J*/J and where J* < 0.
"""
import json, sys
from pathlib import Path
import numpy as np

run = Path(sys.argv[1]); raw = (run / 'QSTARS_ALI_NEWTON_DUMP.bin').read_bytes()
pos, out = 0, []
while pos < len(raw):
    fi, nd = np.frombuffer(raw, '<i4', 2, pos); pos += 8
    fl = float(np.frombuffer(raw, '<f8', 1, pos)[0]); pos += 8
    v = np.frombuffer(raw, '<f8', 5 * nd + nd * nd, pos); pos += 8 * v.size
    J, P = v[:nd], v[nd:2 * nd]; A = v[2 * nd:2 * nd + nd * nd].reshape((nd, nd), order='F')
    chi, scat, eta = v[2 * nd + nd * nd:].reshape(3, nd)
    w, V = np.linalg.eig(A); k = int(np.argmax(np.abs(w))); pv = np.abs(V[:, k]); pv /= pv.max()
    Js = J + np.linalg.solve(np.eye(nd) - A, P - J)
    neg = np.nonzero(Js < 0)[0]
    out.append(dict(freq_index=int(fi), nu_1e15Hz=fl, rho=float(abs(w[k])), rho_real=float(w[k].real),
                    perron_support_depths=(np.nonzero(pv > 0.1)[0] + 1).tolist(),
                    min_A=float(A.min()), depths_scat_gt_chi=(np.nonzero(scat > chi)[0] + 1).tolist(),
                    max_scat_over_chi=float((scat / chi).max()),
                    min_Jstar_over_J=float((Js / J).min()), depths_Jstar_negative=(neg + 1).tolist(),
                    fixed_point_residual_at_start=float(np.max(np.abs(P - J) / np.maximum(P, J)))))
txt = json.dumps(out, indent=1)
if len(sys.argv) > 2:
    Path(sys.argv[2]).write_text(txt + '\n')
print(txt)
