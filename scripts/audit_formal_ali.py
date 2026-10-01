#!/usr/bin/env python3
"""Independent audit of the Report 16 moment blocker and QS_FRM_ALI operator fixtures; never a wind."""
import argparse,json,re
from pathlib import Path
import numpy as np
from qstars.cmfgen_backend import read_spectral,sha256

ROOT=Path(__file__).resolve().parents[1]


def seed_state(path):
    """Report 16 seed frame: moment state versus the formal J the native code held at that moment."""
    x=np.load(path)['state'];R,V,sig,gam=x[:,0],x[:,1],x[:,2],x[:,6]
    formal_gJ=x[:,35]*gam*R**2;moment_gJ=x[:,34];prev_gJ=x[:,28];prev_h=x[:,29]
    face_j=.5*(prev_gJ[:-1]+prev_gJ[1:])
    return dict(file=str(path.relative_to(ROOT)),sha256=sha256(path),
        outer_depths_formal_gJ=formal_gJ[:4].tolist(),outer_depths_moment_gJ=moment_gJ[:4].tolist(),
        max_formal_over_moment=float(np.max(formal_gJ[:4]/np.abs(moment_gJ[:4]))),
        min_sigma_plus_one=float(np.min(sig+1)),velocity_monotone=bool(np.all(np.diff(V)<=0)),
        previous_face_H_over_J_min=float(np.min(prev_h[:-1]/face_j)),
        lag_weight_HS_outer=x[:7,25].tolist(),
        interpretation='Formal-ray J at the outer four depths is ~1e99 while the moment J is ~1e2: the captured '
                       'frame is downstream of an already corrupted restart state (Report 15 runaway), so its signed '
                       'forcing is a symptom, not an origin. Velocity is monotone (sigma+1>0), so per-ray blue-to-red '
                       'transport is well posed there.')


def stale_export(fixed,active):
    out={}
    for label,d in [('fixed_J_settled',fixed),('transfer_active',active)]:
        s=read_spectral(d/'QSTARS_SPECTRAL.bin');es=np.asarray(s.scattering);chi=np.asarray(s.extinction)
        ne=es[len(s.frequency)//2]/6.6524587e-25
        out[label]=dict(directory=str(d.relative_to(ROOT)),electron_density_from_scattering_min=float(ne.min()),
            electron_density_from_scattering_max=float(ne.max()),
            outer_log_chi_std_over_nu_20000_30000=float(np.std(np.log(chi[20000:30000,0]))))
    out['interpretation']=('USE_FIXED_J=T exports stale opacity buffers (constant n_e, frequency-independent chi). '
                           'That export is not a material state and must not be used for transfer analysis.')
    return out


def fixture_results(directory):
    res={}
    for p in sorted(directory.glob('**/out.txt')):
        rows={}
        for line in p.read_text().splitlines():
            if line.startswith('RESULT'):
                parts=line.split();rows[parts[1]]=[float(v) for v in parts[2:]]
        res[str(p.relative_to(ROOT))]=rows
    return res


def resolution(directory):
    P={n:np.loadtxt(directory/f'nd{n}/profile.txt') for n in [40,79,157,313]};ref=P[313];rows={}
    for n in [40,79,157]:
        a=P[n];step=(313-1)//(n-1);rr=ref[::step]
        if not np.allclose(rr[:,0],a[:,0],rtol=1e-12):raise ValueError('non-nested grids')
        rows[n]=dict(formal_ali_max_rel_error_vs_313=float(np.max(abs(a[:,1]/rr[:,1]-1))),
                     native_moment_max_rel_error_vs_313=float(np.max(abs(a[:,2]/rr[:,2]-1))))
    rows['formal_vs_moment_at_313_max_rel']=float(np.max(abs(ref[:,1]/ref[:,2]-1)))
    return rows


def jacobian(directory):
    out={}
    for p in sorted(directory.glob('jacobian_d*.txt')):
        a=np.loadtxt(p);d0=int(re.search(r'd(\d+)',p.stem).group(1))-1
        out[p.name]=dict(
            eta_max_abs_error_over_max=float(np.max(abs(a[:,3]-a[:,2]))/np.max(abs(a[:,2]))),
            eta_local_rel_error=float(abs(a[d0,3]/a[d0,2]-1)),
            chi_max_abs_error_over_max=float(np.max(abs(a[:,5]-a[:,4]))/np.max(abs(a[:,4]))),
            chi_sign_disagreements=int(np.sum(np.sign(a[:,5])!=np.sign(a[:,4]))))
    out['interpretation']=('VAR_JREL_V5 is a fixed-Eddington-factor linearization. About the formal field it '
                           'differs from the exact discrete formal-ALI response at O(1) in optically thin, '
                           'frequency-coupled layers (including sign errors in dJ/dchi inside the opacity front). '
                           'It can serve as an approximate Newton operator; it is not a 1e-6 nonlocal response.')
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    result=dict(
        seed=seed_state(ROOT/'Reports/evidence/moment_closure_20260930/seed_analysis_v2/seed.npz'),
        stale_export=stale_export(ROOT/'results/cmfgen_recovery_T9_fixedJ_ff',ROOT/'results/cmfgen_recovery_T9_refined_full_ff_bootstrap'),
        operator_fixtures=fixture_results(ROOT/'results/formal_ali_probe'),
        operator_resolution=resolution(ROOT/'results/formal_ali_probe'),
        response_fixture=jacobian(ROOT/'results/formal_ali_jacobian'),
        driver_sha256=sha256(Path(__file__)),production_certified=False)
    (a.output/'audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2)[:4000])
