"""Damped Newton with a local, complex-step block Jacobian.

No subtraction of nearby residuals is used for differentiation. Variables
are ordered [v_1..v_N, lnT_0..lnT_{N-1}]; each equation sees neighboring
cells only. Seven colors cover the stencil, independent of zone count.
"""
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
from scipy.optimize import OptimizeResult


def local_jacobian(fun, x):
    n=len(x)//2
    rows=[];cols=[];values=[]
    for kind in (0,1):
        colors=4 if kind == 0 else 3
        for color in range(colors):
            indices=np.arange(color,n,colors)
            xx=x.astype(complex)
            xx[indices+kind*n]+=1e-25j
            derivative=fun(xx).imag/1e-25
            for j in indices:
                if kind == 0:
                    support=np.r_[np.arange(max(0,j-1),min(n,j+2)),
                                  n+np.arange(max(0,j-1),min(n,j+3))]
                else:
                    support=np.r_[np.arange(max(0,j-1),min(n,j+1)),
                                  n+np.arange(max(0,j-1),min(n,j+2))]
                rows.extend(support);cols.extend([j+kind*n]*len(support));values.extend(derivative[support])
    return csc_matrix((values,(rows,cols)),shape=(2*n,2*n))


def newton(fun,x,tolerance,max_iterations=40,jacobian=local_jacobian):
    x=x.copy();nfev=0
    for iteration in range(max_iterations):
        f=fun(x);nfev+=1
        if np.all(np.isfinite(f)) and np.max(abs(f))<=tolerance:
            return OptimizeResult(x=x,nfev=nfev,nit=iteration,success=True,message='Newton residual converged')
        jac=jacobian(fun,x);nfev+=7
        dx=spsolve(jac,-f)
        if not np.all(np.isfinite(dx)): break
        norm=np.linalg.norm(f)
        for k in range(25):
            candidate=x+dx*2.**(-k)
            try: trial=fun(candidate)
            except (ValueError,FloatingPointError): continue
            nfev+=1
            if np.all(np.isfinite(trial)) and np.linalg.norm(trial)<norm*(1-1e-4*2.**(-k)):
                x=candidate;break
        else: break
    return OptimizeResult(x=x,nfev=nfev,nit=iteration,success=False,message='Newton line search or iteration limit')


def flux_jacobian(fun,x):
    n=(len(x)+1)//3
    rows=[];cols=[];values=[]
    for kind,count,offset,colors in ((0,n,0,3),(1,n,n,2),(2,n-1,2*n,2)):
        for color in range(colors):
            indices=np.arange(color,count,colors);xx=x.astype(complex);xx[indices+offset]+=1e-25j
            derivative=fun(xx).imag/1e-25
            for j in indices:
                if kind==0:
                    support=np.r_[np.arange(max(0,j-1),min(n,j+2)),
                                  n+np.arange(j,min(n,j+2)),
                                  2*n+np.arange(max(0,j-1),min(n-1,j+2))]
                elif kind==1:
                    support=np.r_[np.arange(max(0,j-1),min(n,j+1)),n+j,
                                  2*n+np.arange(max(0,j-1),min(n-1,j+1))]
                else:
                    support=np.r_[n+j,n+j+1,2*n+j]
                rows.extend(support);cols.extend([offset+j]*len(support));values.extend(derivative[support])
    return csc_matrix((values,(rows,cols)),shape=(len(x),len(x)))
