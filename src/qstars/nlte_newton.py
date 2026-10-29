"""Offline linear algebra for the CMFGEN cool-closure Newton-Krylov driver.

The native build-18 hooks dump, at one state x:
  * the raw assembled residual R(x)               (QSTARS_RESID_TRACE.bin)
  * every full Jacobian block after equation selection, FIXPOP and the native
    column scaling by POPS                        (QSTARS_BA_TRACE.bin)
  * exact radiation parts of J.(x*s) for requested fractional directions s
                                                  (QSTARS_JV_RAD.bin)
With QS_ALI_LIN=T and QS_LIN_FAC=0 the dumped blocks contain no radiation
response, so  J.(x*s) = BA_direct.s + dSTEQ_rad(s)  exactly (to native rounding).

Everything here is in the native fractional variables: a step s means
POPS -> POPS*(1 - s) (linear) or POPS*exp(-s) (positivity preserving).
"""
from __future__ import annotations

import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.linalg as sla


def read_trace_record(path, index=0):
    """One int32 (count, N, ND) + POPS + VEC record from a QSTARS *_TRACE.bin."""
    data = Path(path).read_bytes()
    off = 0
    for i in range(index + 1):
        cnt, nt, nd = struct.unpack_from('<3i', data, off)
        n = nt * nd
        if off + 12 + 16 * n > len(data):
            raise ValueError(f'{path}: truncated record {i}')
        if i == index:
            pops = np.frombuffer(data, '<f8', n, off + 12).reshape((nt, nd), order='F').copy()
            vec = np.frombuffer(data, '<f8', n, off + 12 + 8 * n).reshape((nt, nd), order='F').copy()
            return pops, vec
        off += 12 + 16 * n
    raise IndexError(index)


def read_ba_blocks(path, nt, nd, diag=2):
    """Return (lower, diagonal, upper) block arrays shaped (nd, nt, nt).

    lower[k] couples equation depth k to variable depth k-1, upper[k] to k+1
    (0-based depths). Missing boundary blocks are zero.
    """
    data = Path(path).read_bytes()
    rec = 12 + 8 * nt * nt
    if len(data) % rec:
        raise ValueError(f'{path}: partial record')
    lo = np.zeros((nd, nt, nt)); di = np.zeros((nd, nt, nt)); up = np.zeros((nd, nt, nt))
    seen = set()
    for off in range(0, len(data), rec):
        k, b, n = struct.unpack_from('<3i', data, off)
        if n != nt:
            raise ValueError('block dimension mismatch')
        blk = np.frombuffer(data, '<f8', nt * nt, off + 12).reshape((nt, nt), order='F')
        key = (k, b)
        if key in seen:
            raise ValueError(f'duplicate block {key}: dump must hold exactly one assembly')
        seen.add(key)
        {diag - 1: lo, diag: di, diag + 1: up}[b][k - 1] = blk
    if len(seen) != 3 * nd - 2:
        raise ValueError(f'expected {3 * nd - 2} blocks, found {len(seen)}')
    return lo, di, up


@dataclass
class BlockTridiag:
    lo: np.ndarray
    di: np.ndarray
    up: np.ndarray

    @property
    def shape(self):
        return self.di.shape[0], self.di.shape[1]

    def matvec(self, s):
        """y[:,k] = lo[k] s[:,k-1] + di[k] s[:,k] + up[k] s[:,k+1]; s shaped (nt, nd)."""
        nd = self.di.shape[0]
        y = np.einsum('kij,jk->ik', self.di, s)
        y[:, 1:] += np.einsum('kij,jk->ik', self.lo[1:], s[:, :-1])
        y[:, :-1] += np.einsum('kij,jk->ik', self.up[:-1], s[:, 1:])
        return y

    def factor(self):
        """Block LU (Thomas). Each pivot block is row/column equilibrated."""
        nd, nt = self.shape
        self._piv = []
        self._g = np.zeros_like(self.up)   # G_k = Piv_k^{-1} up_k
        prev_g = None
        for k in range(nd):
            p = self.di[k].copy()
            if k > 0:
                p -= self.lo[k] @ prev_g
            r = 1.0 / np.maximum(np.max(np.abs(p), axis=1), 1e-300)
            ps = p * r[:, None]
            c = 1.0 / np.maximum(np.max(np.abs(ps), axis=0), 1e-300)
            ps = ps * c[None, :]
            lu = sla.lu_factor(ps, check_finite=True)
            self._piv.append((lu, r, c))
            if k < nd - 1:
                prev_g = self._solve_piv(k, self.up[k])
                self._g[k] = prev_g
        return self

    def _solve_piv(self, k, rhs):
        lu, r, c = self._piv[k]
        rr = rhs * (r[:, None] if rhs.ndim == 2 else r)
        x = sla.lu_solve(lu, rr, check_finite=False)
        return x * (c[:, None] if x.ndim == 2 else c)

    def solve(self, b):
        nd, nt = self.shape
        y = np.zeros_like(b)
        prev = None
        for k in range(nd):
            rhs = b[:, k].copy()
            if k > 0:
                rhs -= self.lo[k] @ prev
            prev = self._solve_piv(k, rhs)
            y[:, k] = prev
        x = np.zeros_like(b)
        x[:, nd - 1] = y[:, nd - 1]
        for k in range(nd - 2, -1, -1):
            x[:, k] = y[:, k] - self._g[k] @ x[:, k + 1]
        return x


def read_jv_rad(path):
    data = Path(path).read_bytes()
    nv, nt, nd = struct.unpack_from('<3i', data, 0)
    n = nt * nd * nv
    s = np.frombuffer(data, '<f8', n, 12).reshape((nt, nd, nv), order='F').copy()
    d = np.frombuffer(data, '<f8', n, 12 + 8 * n).reshape((nt, nd, nv), order='F').copy()
    return s, d


def write_fortran(path, a):
    """Raw float64 column-major array (ndarray.tofile would always write C order)."""
    Path(path).write_bytes(np.asarray(a, dtype='<f8').tobytes(order='F'))


def read_fortran(path, shape):
    return np.frombuffer(Path(path).read_bytes(), '<f8').reshape(shape, order='F').copy()


def write_directions(path, dirs):
    """dirs: (nt, nd, nv) fractional directions -> QSTARS_JV_DIRS stream."""
    dirs = np.asarray(dirs, dtype='<f8')
    nt, nd, nv = dirs.shape
    with open(path, 'wb') as f:
        f.write(struct.pack('<3i', nv, nt, nd))
        f.write(np.asfortranarray(dirs).tobytes(order='F'))


def gmres(apply, b, m, precond=None, tol=1e-12, x_in=None):
    """Right-preconditioned GMRES(m) on flattened vectors; returns (x, history).

    apply(v) -> J v,  precond(v) -> M^{-1} v. Modified Gram-Schmidt twice.
    """
    n = b.size
    beta = np.linalg.norm(b)
    V = np.zeros((n, m + 1)); H = np.zeros((m + 1, m)); Z = np.zeros((n, m))
    V[:, 0] = b / beta
    hist = [1.0]
    for j in range(m):
        z = precond(V[:, j]) if precond else V[:, j]
        Z[:, j] = z
        w = apply(z)
        for _ in range(2):
            for i in range(j + 1):
                h = V[:, i] @ w
                H[i, j] += h
                w -= h * V[:, i]
        H[j + 1, j] = np.linalg.norm(w)
        if H[j + 1, j] > 0:
            V[:, j + 1] = w / H[j + 1, j]
        e1 = np.zeros(j + 2); e1[0] = beta
        y, *_ = np.linalg.lstsq(H[:j + 2, :j + 1], e1, rcond=None)
        res = np.linalg.norm(H[:j + 2, :j + 1] @ y - e1) / beta
        hist.append(float(res))
        if res < tol or H[j + 1, j] == 0:
            break
    x = Z[:, :len(y)] @ y
    return x, hist
