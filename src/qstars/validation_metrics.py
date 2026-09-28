"""Profile norms on a common mass domain, independent of evolution residuals."""
import numpy as np
from scipy.interpolate import PchipInterpolator

def interpolant_norms(xa,ya,xb,yb):
    """Exact extrema and polynomial quadrature of two piecewise cubic interpolants.

    The interpolants are a representation choice, not the exact physical
    profile. Their common center domain excludes the stated boundary mass.
    """
    a=PchipInterpolator(xa,ya);b=PchipInterpolator(xb,yb)
    lo=max(xa[0],xb[0]);hi=min(xa[-1],xb[-1])
    knots=np.unique(np.r_[lo,xa[(xa>lo)&(xa<hi)],xb[(xb>lo)&(xb<hi)],hi])
    left=knots[:-1];width=np.diff(knots)
    def shifted_derivative(interp):
        i=np.clip(np.searchsorted(interp.x,left,side='right')-1,0,len(interp.x)-2)
        c=interp.derivative().c[:,i];z=left-interp.x[i]
        return np.array([c[0]*width**2,(2*c[0]*z+c[1])*width,(c[0]*z+c[1])*z+c[2]])
    A,B,C=shifted_derivative(a)-shifted_derivative(b)
    scale=np.maximum.reduce([abs(A),abs(B),abs(C),np.full_like(A,1e-300)])
    quadratic=abs(A)>1e-14*scale
    discriminant=B*B-4*A*C;valid=quadratic&(discriminant>=0)
    locations=[knots]
    for sign in (-1,1):
        root=(-B[valid]+sign*np.sqrt(discriminant[valid]))/(2*A[valid])
        inside=(root>0)&(root<1)
        locations.append(left[valid][inside]+width[valid][inside]*root[inside])
    linear=(~quadratic)&(abs(B)>1e-14*scale)
    root=-C[linear]/B[linear];inside=(root>0)&(root<1)
    locations.append(left[linear][inside]+width[linear][inside]*root[inside])
    points=np.concatenate(locations);linf=float(max(abs(a(points)-b(points))))
    nodes,weights=np.polynomial.legendre.leggauss(4)
    grid=left[:,None]+.5*width[:,None]*(nodes+1)
    integral=np.sum(.5*width[:,None]*weights*(a(grid)-b(grid))**2)
    return linf,float(np.sqrt(integral/(hi-lo)))
