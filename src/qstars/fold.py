"""Frozen-entropy hydrostatic continuation in M_BH: classification of a Ball track's end.

Plan endpoint.classification_protocol: when the termination trigger of a
Ball track fires (dt below t_dyn, or 12 consecutive retries), the
hydrostatic problem at frozen specific entropy is continued in
lambda = ln(M_BH/M_BH,terminal) from the last accepted state, with the
cavity and photosphere boundaries active (pseudo-arclength). A fold is a
sign change of det J (J the Jacobian in the state at fixed M_BH) where the
smallest singular value of J vanishes as a square root of the distance in
M_BH. Only a fold within 2% of the terminal M_BH is labelled
model_domain_exit ("no quasi-static continuation", Ball's operational
criterion); every other outcome is numerical_failure and is never relabelled.

Frozen problem (the dt -> 0 limit of the quasi-static step, whose first-law
rows become s(rho_k, T_k) = s_k). Gas keeps the specific entropy of its mass
coordinate mu above the base of the reference state: s(mu) is the PCHIP
interpolant through the reference cell centres, read at each cell's mass
centre. A continuous swept mass S (an unknown) leaves the base into the
cavity and the black hole grows out of the cavity, booked as in the step:
M_cav = M_cav,ref + S - dM_BH/(1 - eps), m_0 = M_BH + M_cav and
dm_k = (M_env,ref - S) f_k with the reference mass fractions f_k, so the
gravitating mass changes only by the mass defect eps dM_BH/(1 - eps).
Unknowns (2N+3): ln(r/r_ref) on faces 0..N, ln(T/T_ref) on cells,
ln rho_ph and S/M_env,ref. Rows: static balance at the face radii
4 pi r_j^2 (p_j-1 - p_j)/w_j = G m_j/r_j^2 on faces 1..N with p_N = p_ext (the
endpoint momentum rows of quasistatic.py); frozen entropy per cell, scaled
by c_p; the piston ln r0 = ln(r_B(M_BH, c_b)/b); the cavity mass (thesis
3.7); the photosphere row (thesis 2.16, L_N from the grey relation of the
last cell, as in the step) or, without the Ball photosphere, ln rho_ph
pinned and p_ext = boundary.pressure_cgs. An accepted quasi-static state is
a solution at lambda = 0 to its Newton tolerance.

Continuation (Arc): pseudo-arclength with the bordered tangent; at the
first sign change of t_lambda the turning point is bisected in arclength,
and sigma_min is probed at arclengths h0 2^-k behind it; its exponent in
|lambda_fold - lambda| must be 1/2 (sqrt_law). sigma_min comes from inverse
iteration with the SuperLU factors; the smallest |U_ii| pivot is reported
too but depends on the pivoting order. classify() searches +5% in M_BH,
then -2%.

For an exact Gamma = 4/3 polytropic gas the frozen problem has the
homology r -> a r, rho -> rho/a^3, T -> T/a as an exact symmetry (the piston
and the cavity mass scale the same way), so J is singular along the whole
branch and has no fold to detect: n = 3 is the degenerate limit of the
loaded-polytrope family. The identity test (spec I4) therefore uses
Gamma = 1 + 1/n with n < 3, whose fold is at M_BH/M = reference.ball_limit(n).
"""
import copy,math
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.sparse import csc_matrix,hstack,vstack
from scipy.sparse.linalg import splu
from .constants import G,PI
from .state import State
from .nonlinear import pattern_jacobian
from .mass_transport import DomainError
from .transport import luminosities
from . import cavity,photosphere

TRIAL_ERRORS=(ValueError,FloatingPointError,ZeroDivisionError,DomainError,RuntimeError)
DEFAULTS=dict(ds=2e-2,ds_min=1e-8,ds_max=.5,max_points=200,tolerance=1e-9,roundoff_residual=1e-8,max_iterations=12,
              after_fold=2,refine=16,probe_fraction=.1,probe_max=.05,search=.05,window=.02,sqrt_points=4,
              sqrt_exponent_tolerance=.1)


class Profile:
    """PCHIP interpolant of (x, y) evaluated at complex points (exact complex-step derivatives).

    The interval is chosen from the real part; outside the breakpoints the
    end polynomials extrapolate.
    """
    def __init__(self,x,y):
        p=PchipInterpolator(np.asarray(x,dtype=float),np.asarray(y,dtype=float),extrapolate=True)
        self.x,self.c=p.x,p.c
    def __call__(self,z):
        z=np.asarray(z)
        k=np.clip(np.searchsorted(self.x,np.real(z),side='right')-1,0,len(self.x)-2)
        d=z-self.x[k];c=self.c[:,k]
        return ((c[0]*d+c[1])*d+c[2])*d+c[3]


class Frozen:
    """Frozen-entropy hydrostatic residual F(x, lambda) about a reference (accepted) state."""
    def __init__(self,state,eos,opacity,config):
        from .quasistatic import settings
        s=settings(config)
        self.eos,self.opacity,self.config,self.cav=eos,opacity,config,s['cavity']
        if not self.cav['active']:raise ValueError('The frozen-entropy continuation needs the cavity boundary')
        self.n=n=len(state.dm)
        if n<3:raise ValueError('The frozen-entropy continuation needs at least three cells')
        q=eos.evaluate(state.rho,state.T)
        self.r,self.T=np.array(state.r,dtype=float),np.array(state.T,dtype=float)
        self.M,self.Mcav=float(state.bh_mass),float(state.cavity_mass)
        self.Menv=math.fsum(state.dm);self.f=state.dm/self.Menv
        centre=np.cumsum(state.dm)-.5*state.dm
        self.fc=centre/self.Menv
        self.entropy=Profile(centre,np.real(q['s']))
        self.cp=np.real(q['cp']).copy()
        self.dm0=float(state.dm[0])
        self.eps=config['engine']['efficiency']
        self.g=G*state.enclosed/state.r**2
        self.photosphere=config['boundary']['thermal']=='stars_photosphere'
        self.pext=float(config['boundary']['pressure_cgs'])
        self.bc_config=copy.deepcopy(config);self.bc_config['transport']['convection']='off'
        if self.photosphere:
            rho_ph=state.photosphere_density
            if rho_ph is None:
                L=float(np.real(luminosities(state,eos,opacity,self.bc_config,0.)[0][-1]))
                rho_ph=photosphere.solve(L,float(state.r[-1]),float(state.enclosed[-1]),eos,opacity,
                                         bracket=s['quasistatic']['photosphere_bracket'])
            self.lnrho_ph=math.log(rho_ph)
        else:self.lnrho_ph=math.log(state.photosphere_density) if state.photosphere_density else 0.
        self.pattern=pattern(n)

    def x0(self):
        x=np.zeros(2*self.n+3);x[2*self.n+1]=self.lnrho_ph
        return x

    def evaluate(self,x,lam,parts=False):
        n,eos=self.n,self.eos
        y,th=x[:n+1],x[n+1:2*n+1]
        lnrho_ph,S=x[2*n+1],x[2*n+2]*self.Menv
        M=self.M*np.exp(lam)
        accreted=(M-self.M)/(1-self.eps)
        Mcav=self.Mcav+S-accreted
        m0=M+Mcav
        dm=(self.Menv-S)*self.f
        if np.any(np.real(dm)<=0):raise ValueError('Swept mass exceeds the envelope')
        r=self.r*np.exp(y)
        if np.any(np.diff(np.real(r))<=0) or not np.real(r[0])>0:raise ValueError('Trial geometry has crossing radial faces')
        V=4*PI/3*np.diff(r)*(r[1:]**2+r[1:]*r[:-1]+r[:-1]**2)
        T=self.T*np.exp(th)
        rho=dm/V
        q=eos.evaluate(rho,T);p=q['p']
        st=State(0.,m0,dm,r,np.zeros(n+1),T,cell_volume=V)
        m=st.enclosed
        B=cavity.boundary_state(rho[:2],T[:2],p[:2],r[0],m0,dm[0],eos,self.cav['boundary_state'])
        if self.photosphere:
            LN=luminosities(st,eos,self.opacity,self.bc_config,0.)[0][-1]
            ph=photosphere.evaluate(np.exp(lnrho_ph),LN,r[-1],m[-1],eos,self.opacity)
            pext,row_ph=ph['pressure'],ph['residual']
        else:
            LN=ph=None;pext=self.pext+0.*lnrho_ph;row_ph=lnrho_ph-self.lnrho_ph
        w=.5*np.r_[dm[0],dm[:-1]+dm[1:],dm[-1]]
        dp=p-np.r_[p[1:],pext]
        momentum=(G*m[1:]/r[1:]**2-4*PI*r[1:]**2*dp/w[1:])/self.g[1:]
        centre=S+(self.Menv-S)*self.fc
        entropy=(q['s']-self.entropy(centre))/self.cp
        piston=np.log(r[0])-np.log(cavity.bondi_radius(M,B['cs'],self.cav['b']))
        row_cav=(Mcav-cavity.cavity_mass(B['rho'],r[0],self.cav['profile']))/self.dm0
        f=np.r_[momentum,entropy,piston,row_cav,row_ph]
        if not parts:return f
        return f,dict(r=r,T=T,rho=rho,dm=dm,m=m,p=p,q=q,boundary=B,photosphere=ph,LN=LN,pext=pext,
                      bh_mass=M,cavity_mass=Mcav,swept=S,total_mass=m[-1])

    def jacobian(self,x,lam,step=1e-25):
        """J = dF/dx on the declared pattern (coloured complex step) and the dense column dF/dlambda."""
        J=pattern_jacobian(lambda z:self.evaluate(z,lam),x,self.pattern)
        Fl=np.asarray(self.evaluate(np.asarray(x,dtype=complex),lam+step*1j)).imag/step
        return J,Fl

    def state(self,x,lam):
        """The frozen-problem solution as a State (diagnostics; no ledger)."""
        _,d=self.evaluate(x,lam,parts=True)
        real=lambda a:np.real(a).astype(float)
        st=State(0.,float(np.real(d['bh_mass'])),real(d['dm']),real(d['r']),np.zeros(self.n+1),real(d['T']),
                 cell_volume=4*PI/3*np.diff(real(d['r']))*(real(d['r'])[1:]**2+real(d['r'])[1:]*real(d['r'])[:-1]+real(d['r'])[:-1]**2),
                 cavity_mass=float(np.real(d['cavity_mass'])),kinetic_tracked=False)
        if self.photosphere:st.photosphere_density=float(np.exp(np.real(x[2*self.n+1])))
        return st


def pattern(n):
    """Declared sparsity of J (2N+3 square); the swept-mass column is dense (every dm_k depends on S)."""
    Y=lambda j:j;T=lambda k:n+1+k;PH,S=2*n+1,2*n+2
    base=[Y(0),Y(1),Y(2),T(0),T(1),S];outer=[Y(n-1),Y(n),T(n-1),PH,S]
    rows=[];cols=[]
    def add(i,js):
        js=sorted(set(js));rows.extend([i]*len(js));cols.extend(js)
    for j in range(1,n+1):
        add(j-1,[Y(j-1),Y(j),T(j-1),S]+([Y(j+1),T(j)] if j<n else outer))
    for k in range(n):
        add(n+k,[Y(k),Y(k+1),T(k),S])
    add(2*n,base);add(2*n+1,base);add(2*n+2,outer)
    return csc_matrix((np.ones(len(rows)),(rows,cols)),shape=(2*n+3,2*n+3))


def parity(permutation):
    """Sign of a permutation (cycle count)."""
    perm=np.asarray(permutation);seen=np.zeros(len(perm),bool);sign=1
    for i in range(len(perm)):
        if seen[i]:continue
        j=i;length=0
        while not seen[j]:
            seen[j]=True;j=perm[j];length+=1
        if length%2==0:sign=-sign
    return sign

def monitors(J,power_iterations=40):
    """Fold monitors of J: sign det J (SuperLU: U diagonal and both permutation parities),
    the smallest |U_ii| pivot, and the smallest singular value from inverse iteration
    on J^T J with the same factors."""
    lu=splu(csc_matrix(J))
    u=lu.U.diagonal()
    sign=int(np.prod(np.sign(u)))*parity(lu.perm_r)*parity(lu.perm_c)
    n=J.shape[0]
    v=np.ones(n)+.01*np.cos(np.arange(n));v/=np.linalg.norm(v);growth=1.
    for _ in range(power_iterations):
        w=lu.solve(lu.solve(v,trans='T'))
        growth=float(np.linalg.norm(w))
        if not (growth>0 and np.isfinite(growth)):break
        v=w/growth
    return dict(det_sign=sign,min_pivot=float(np.min(abs(u))),max_pivot=float(np.max(abs(u))),
                min_singular=float(1/math.sqrt(growth)) if growth>0 and np.isfinite(growth) else 0.,
                log_abs_det=float(np.sum(np.log(abs(u)))))


def solve(problem,x,lam,tolerance=None,max_iterations=None):
    """Damped Newton on F(x, lambda) = 0 at fixed lambda (SuperLU on the pattern Jacobian)."""
    tolerance=DEFAULTS['tolerance'] if tolerance is None else tolerance
    max_iterations=DEFAULTS['max_iterations']*2 if max_iterations is None else max_iterations
    x=np.array(x,dtype=float)
    f=problem.evaluate(x,lam)
    for iteration in range(max_iterations+1):
        residual=float(np.max(abs(f)))
        if residual<=tolerance:return x,dict(success=True,iterations=iteration,residual=residual)
        if iteration==max_iterations:break
        J,_=problem.jacobian(x,lam)
        dx=splu(csc_matrix(J)).solve(-f);norm=np.linalg.norm(f)
        for k in range(25):
            try:trial=problem.evaluate(x+dx*2.**-k,lam)
            except TRIAL_ERRORS:continue
            if np.all(np.isfinite(trial)) and np.linalg.norm(trial)<norm*(1-1e-4*2.**-k):
                x,f=x+dx*2.**-k,trial;break
        else:break
    return x,dict(success=False,iterations=iteration,residual=float(np.max(abs(f))))


class Arc:
    """Pseudo-arclength continuation of F(x, lambda) = 0.

    Arclength metric: |dz|^2 = |dx|^2/len(x) + dlambda^2. The tangent t solves
    the bordered system [J, F_lambda; t_prev] t = [0; 1] (normalised), so it
    stays oriented through a fold, where t_lambda changes sign with det J.
    """
    def __init__(self,problem,s=None):
        self.p=problem;self.s={**DEFAULTS,**(s or {})}
        self.w=1./(2*problem.n+3)

    def dot(self,a,b):
        return self.w*float(np.dot(a[:-1],b[:-1]))+float(a[-1]*b[-1])

    def bordered(self,J,Fl,t):
        row=csc_matrix(np.r_[self.w*t[:-1],t[-1]][None,:])
        return csc_matrix(vstack([hstack([csc_matrix(J),csc_matrix(Fl[:,None])]),row]))

    def tangent(self,z,t_prev=None):
        x,lam=z[:-1],z[-1]
        J,Fl=self.p.jacobian(x,lam)
        if t_prev is None:
            v=splu(csc_matrix(J)).solve(-Fl);t=np.r_[v,1.]
        else:
            rhs=np.zeros(len(z));rhs[-1]=1.
            t=splu(self.bordered(J,Fl,t_prev)).solve(rhs)
        t=t/math.sqrt(self.dot(t,t))
        return t,J,Fl

    def correct(self,z_pred,t):
        """Newton on [F(x, lambda); <t, z - z_pred>] = 0."""
        s=self.s;z=np.array(z_pred,dtype=float)
        def G(z):
            return np.r_[self.p.evaluate(z[:-1],z[-1]),self.dot(t,z-z_pred)]
        g=G(z)
        for iteration in range(s['max_iterations']+1):
            residual=float(np.max(abs(g)))
            if residual<=s['tolerance']:return z,iteration
            if iteration==s['max_iterations']:break
            J,Fl=self.p.jacobian(z[:-1],z[-1])
            dz=splu(self.bordered(J,Fl,t)).solve(-g);norm=np.linalg.norm(g)
            for k in range(20):
                try:trial=G(z+dz*2.**-k)
                except TRIAL_ERRORS:continue
                if np.all(np.isfinite(trial)) and np.linalg.norm(trial)<norm*(1-1e-4*2.**-k):
                    z,g=z+dz*2.**-k,trial;break
            else:
                # Roundoff exit: a stalled line search at the residual floor of the scaled rows.
                if residual<=s['roundoff_residual']:return z,iteration
                break
        raise RuntimeError(f'Arclength corrector failed (residual {float(np.max(abs(g))):.3g})')

    def point(self,z,t,J,ds,iterations):
        m=monitors(J)
        _,d=self.p.evaluate(z[:-1],z[-1],parts=True)
        return dict(lam=float(z[-1]),bh_mass=float(np.real(d['bh_mass'])),total_mass=float(np.real(d['total_mass'])),
                    bh_fraction=float(np.real(d['bh_mass']/d['total_mass'])),cavity_mass=float(np.real(d['cavity_mass'])),
                    swept=float(np.real(d['swept'])),inner_radius=float(np.real(d['r'][0])),radius=float(np.real(d['r'][-1])),
                    base_density=float(np.real(d['boundary']['rho'])),t_lambda=float(t[-1]),ds=float(ds),
                    iterations=int(iterations),**m)

    def run(self,x,lam=0.,direction=1,lam_limit=None):
        """Continue from a solution (x, lam). Returns dict(points, fold, stop).

        At the first sign change of t_lambda the fold is refined by bisection
        and probed on its approach side (probe), and the run stops after_fold
        points later.
        """
        s=self.s;z=np.r_[np.asarray(x,dtype=float),lam]
        t,J,_=self.tangent(z);t=t*direction
        points=[self.point(z,t,J,0.,0)];ds=s['ds'];fold=None;beyond=0;stop='maximum points'
        while len(points)<s['max_points']:
            try:
                znew,iterations=self.correct(z+ds*t,t)
                tnew,Jnew,_=self.tangent(znew,t)
            except TRIAL_ERRORS as exc:
                ds*=.5
                if ds<s['ds_min']:
                    stop=f'step size underflow: {type(exc).__name__}: {exc}';break
                continue
            points.append(self.point(znew,tnew,Jnew,ds,iterations))
            if fold is None and np.sign(tnew[-1])!=np.sign(t[-1]):
                fold=self.refine(z,t,ds)
                fold['bracket']=[len(points)-2,len(points)-1]
            z,t=znew,tnew
            if fold is not None:
                beyond+=1
                if beyond>=s['after_fold']:stop='fold passed';break
            if lam_limit is not None and abs(z[-1]-lam)>=lam_limit:
                stop='lambda limit';break
            ds=min(ds*1.5,s['ds_max']) if iterations<=3 else ds*.7 if iterations>6 else ds
        return dict(points=points,fold=fold,stop=stop)

    def refine(self,z,t,ds):
        """Bisect the arclength from (z, t) for t_lambda = 0 (the fold is the lambda extremum there),
        then probe the approach side at arclengths h0 2^-k for the square-root law."""
        lo,hi=0.,ds;best=None;sign=np.sign(t[-1]);message=None
        for _ in range(self.s['refine']):
            mid=.5*(lo+hi)
            try:
                zm,_=self.correct(z+mid*t,t)
                tm,Jm,_=self.tangent(zm,t)
            except TRIAL_ERRORS as exc:
                message=f'{type(exc).__name__}: {exc}';break
            if np.sign(tm[-1])==sign:lo=mid
            else:hi=mid
            best=(zm,tm,Jm)
        if best is None:
            return dict(lam=float(z[-1]),refined=False,refine_message=message,probe=[],arclength_bracket=float(ds))
        zm,tm,Jm=best
        out=self.point(zm,tm,Jm,hi-lo,0)
        out.update(refined=message is None,refine_message=message,arclength_bracket=float(hi-lo),
                   probe=self.probe(zm,tm,min(self.s['probe_fraction']*ds,self.s['probe_max'])))
        return out

    def probe(self,z,t,h0):
        """Points at arclength -h0 2^-k (k < sqrt_points) behind the fold: (lambda_fold - lambda, sigma_min)."""
        out=[]
        for k in range(self.s['sqrt_points']):
            h=h0*2.**-k
            try:
                zp,_=self.correct(z-h*t,t)
                J,_=self.p.jacobian(zp[:-1],zp[-1])
            except TRIAL_ERRORS:break
            out.append(dict(arclength=float(h),lam=float(zp[-1]),distance=float(abs(z[-1]-zp[-1])),
                            min_singular=monitors(J)['min_singular']))
        return out


def sqrt_law(fold):
    """Exponent a of sigma_min ~ |lambda_fold - lambda|^a from the probe behind the fold (1/2 at a generic fold)."""
    use=[p for p in fold.get('probe',[]) if p['distance']>0 and p['min_singular']>0]
    if len(use)<2:return None,len(use)
    x=np.log([p['distance'] for p in use]);y=np.log([p['min_singular'] for p in use])
    return float(np.polyfit(x,y,1)[0]),len(use)


def classify(state,eos,opacity,config,s=None,log=None):
    """Plan endpoint protocol from the last accepted state.

    Continues in +lambda up to `search` (ln M_BH); if no fold is found there,
    continues in -lambda up to `window`. The fold is confirmed when det J
    changes sign at a t_lambda turning point, sigma_min vanishes with
    exponent 1/2 +- sqrt_exponent_tolerance, and |M_fold/M_terminal - 1| <=
    window. The N = 400 vs 800 agreement is a campaign-level check.
    """
    s={**DEFAULTS,**(s or {})}
    problem=Frozen(state,eos,opacity,config)
    x,info=solve(problem,problem.x0(),0.)
    out=dict(terminal_bh_mass=float(state.bh_mass),zones=problem.n,start_residual=info['residual'],settings=s,
             runs=[],fold=None,fold_confirmed=False)
    if not info['success']:
        out.update(label='numerical_failure',reason=f"frozen-entropy problem does not hold at the last accepted state (residual {info['residual']:.3g})")
        return out
    arc=Arc(problem,s)
    for direction,limit in ((1,s['search']),(-1,s['window'])):
        try:run=arc.run(x,0.,direction,limit)
        except TRIAL_ERRORS as exc:
            run=dict(points=[],fold=None,stop=f'continuation could not start: {type(exc).__name__}: {exc}')
        out['runs'].append(dict(direction=direction,limit=limit,stop=run['stop'],points=run['points'],fold=run['fold']))
        if log:log(direction,run)
        if run['fold'] is not None:
            fold=dict(run['fold'])
            exponent,used=sqrt_law(run['fold'])
            pts=run['points'];b=run['fold']['bracket']
            fold.update(direction=direction,det_sign_change=pts[b[0]]['det_sign']!=pts[b[1]]['det_sign'],
                        sqrt_exponent=exponent,sqrt_points=used,
                        relative_offset=fold['bh_mass']/state.bh_mass-1)
            fold['sqrt_law']=exponent is not None and abs(exponent-.5)<=s['sqrt_exponent_tolerance']
            fold['within_window']=abs(fold['relative_offset'])<=s['window']
            out['fold']=fold
            break
    fold=out['fold']
    if fold is None:
        out.update(label='numerical_failure',reason='no fold of the frozen-entropy problem within the search range: '
                   +'; '.join(f"direction {r['direction']:+d}: {r['stop']}" for r in out['runs']))
        return out
    confirmed=bool(fold['det_sign_change'] and fold['sqrt_law'] and fold['within_window'])
    out['fold_confirmed']=confirmed
    if confirmed:out.update(label='model_domain_exit',reason='no quasi-static continuation (Ball\'s operational criterion): '
                            f"frozen-entropy fold at M_BH = {fold['bh_mass']:.6g} g, {100*fold['relative_offset']:+.3f}% of the terminal mass")
    else:
        missing=[k for k in ('det_sign_change','sqrt_law','within_window') if not fold[k]]
        out.update(label='numerical_failure',reason=f"fold found but not confirmed ({', '.join(missing)} failed)")
    return out
