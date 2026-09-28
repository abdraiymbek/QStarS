import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

def surface(x,y): return y[0]
surface.terminal=True; surface.direction=-1

def ball(x0,n=3):
 pbh=x0/(2*n); m0=pbh+2*x0**3/3
 def rhs(x,y): return [-y[1]/x**2,x**2*max(y[0],0)**n]
 s=solve_ivp(rhs,[x0,1e4],[1,m0],events=surface,rtol=1e-10,atol=1e-12,max_step=.2)
 return pbh/s.y[1,-1]
for n in (2,3,4):
 s=minimize_scalar(lambda x:-ball(x,n),bounds=(.1,2),method='bounded');print('ball',n,s.x,-s.fun)

def corr(K,N=1):
 def rhs(x,y): return [-y[1]/(K*x*x),x*x*max(y[0],0)**3]
 s=solve_ivp(rhs,[1,1e5],[6*N/K,1],events=surface,rtol=1e-10,atol=1e-12,max_step=10)
 return (1-2/3*(6*N/K)**3)/s.y[1,-1]
for N in (1/3,.5,1.,2.):
 s=minimize_scalar(lambda K:-corr(K,N),bounds=(5.242*N,50*N),method='bounded');print('corr',N,s.x,-s.fun)

def sc(x,y):
 m,z=y
 return [x*x*z**3,4*z/(3*x)-m*z**3/x**(2/3)]
x0=1e-5; coef=3*np.sqrt(3)/(5*np.sqrt(2)); z0=(coef*2.5)**(1/3)*x0**(-1/6)
def entropy_min(x,y): return 3*y[0]*y[1]**2/x**(2/3)-16/(3*x)
entropy_min.terminal=True; entropy_min.direction=1
s=solve_ivp(sc,[x0,2],[1+coef*x0**2.5,z0],events=entropy_min,rtol=1e-11,atol=1e-13,max_step=.002)
x=s.t[-1];m,z=s.y[:,-1]; print('sat',x,1/m,(x/m)**(5/3))
