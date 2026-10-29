"""EXPLORATORY (scratch, DRYRUN only, not part of the proposed patch): Newton start of a failed step by continuation
in the step's own dt on the same post-event mesh. Monkeypatches quasistatic.Problem.solve."""
import numpy as np
def install(qs):
    base_solve=qs.Problem.solve
    def solve(self):
        x,info=base_solve(self)
        if info['success'] or self.mode!='evolve' or getattr(self,'_ladder',False):return x,info
        q=self.s['quasistatic'];z=None;n=self.n;dt=self.dt
        for dtk in np.geomspace(dt*1e-5,dt,11)[:-1]:
            pred=qs.predict(self.old,dtk,self.eos,self.config,'evolve',None,self.s)
            P=qs.Problem(self.old,self.post,dtk,self.eos,self.opacity,self.config,'evolve',None,pred,self.s);P._ladder=True
            start=P.x0 if z is None else np.r_[z[:3*n+2],P.x0[3*n+2]]
            z,i=qs.newton(P.evaluate,start,P.pattern,q['tolerance'],q['max_iterations'],q['polish_iterations'],q['roundoff_residual'],q['roundoff_step'])
            if not i['success']:return x,dict(info,message=info['message']+f'; dt ladder failed at {dtk:.3g} s')
        start=np.r_[z[:3*n+2],self.x0[3*n+2]]
        x2,final=qs.newton(self.evaluate,start,self.pattern,q['tolerance'],q['max_iterations'],q['polish_iterations'],q['roundoff_residual'],q['roundoff_step'])
        if final['success']:return x2,dict(final,ladder=True,iterations=info['iterations']+final['iterations']+30)
        return x,dict(info,message=info['message']+'; dt ladder final solve failed')
    qs.Problem.solve=solve
