"""Wind launch closure and energy demand, separate from supplied material."""
import numpy as np
from .constants import C,G,MSUN,YEAR
from .engine import eddington
from .controls import value
from .mass_transport import DomainError

NODE_ASSUMPTION=("CMFGEN-derived wind node: mdot and terminal speed are explicit (old-state) node values, held within "
    "|dln x|<=drift_tolerance of the node base data or piecewise log-linear along an ordered node sequence; no "
    "extrapolation. The envelope's own outer thermal/pressure boundary is not replaced by the node atmosphere; "
    "base-face and energy defects are diagnostics, not ledger terms.")

def selected(state,config,luminosity):
    """Node values for base data (M, R, L, rho, T) of the outer cell/face.

    One node: values held while every |dln x| <= delta. A sequence: Euclidean
    projection onto each ln-space segment (clamped to it), log-linear values
    there; the state must lie within delta (max norm) of the closest point.
    """
    from .wind_node import load,BASE
    nodes=load(config['node_path'])
    x=np.array([state.enclosed[-1],state.r[-1],luminosity,state.rho[-1],state.T[-1]],dtype=float)
    if not np.all(np.isfinite(x)) or np.any(x<=0):
        raise DomainError('Wind node base data need positive finite mass, radius, luminosity, density and temperature')
    x=np.log(x)
    P=np.log([[n[k] for k in BASE] for n in nodes])
    best=(np.max(abs(x-P[0])),0,0.)
    for i in range(len(nodes)-1):
        d=P[i+1]-P[i];s=float(np.clip(np.dot(x-P[i],d)/np.dot(d,d),0.,1.))
        best=min(best,(np.max(abs(x-P[i]-s*d)),i,s),key=lambda b:b[0])
    drift,i,s=best
    if s==1.:i,s=i+1,0.  # exact node values at a sequence vertex
    if drift>config['drift_tolerance']:
        raise DomainError(f"State left the CMFGEN wind node neighbourhood: max |dln base| {drift:.4g} > "
                          f"drift_tolerance {config['drift_tolerance']:g}")
    a,b=nodes[i],nodes[min(i+1,len(nodes)-1)]
    get=lambda k:float(np.exp((1-s)*np.log(a[k])+s*np.log(b[k]))) if s else float(a[k])
    # Log-linear in every output, so Teff=(L/(4 pi sigma R^2))^(1/4) stays exact.
    return dict(rate=get('mdot_g_s'),terminal_speed=get('v_inf_cm_s'),emergent_luminosity=get('L_emergent_erg_s'),
                photospheric_radius=get('R_phot_cm'),emergent_effective_temperature=get('Teff_emergent_K'),
                base_luminosity=get('L_base_erg_s'),max_abs_dln_drift=float(drift),
                node_id=a['node_id'] if not s else f"{a['node_id']}->{b['node_id']}@{s:.6f}")

def evaluate(state,eos,opacity,config,luminosity):
    if config["model"]=="selected_nodes": return selected(state,config,luminosity)
    # Raw config value keeps the earlier modes bitwise unchanged.
    return dict(rate=rate(state,eos,opacity,config,luminosity),terminal_speed=config["terminal_speed_cm_s"])

def state_values(state,eos,opacity,config):
    """Node values at an accepted state, from the luminosity history_row reports."""
    from .engine import rates
    from .transport import luminosities
    L,_=luminosities(state,eos,opacity,config,rates(state,eos,config["engine"])["luminosity"])
    return evaluate(state,eos,opacity,config["wind"],max(0.,float(L[-1])))

def rate(state,eos,opacity,config,luminosity):
    mode = config["model"]
    if mode == "off": return 0.
    if mode == "prescribed":
        return value(config["rate_msun_per_year"],state.time/YEAR)*MSUN/YEAR
    if mode == "selected_nodes": return selected(state,config,luminosity)["rate"]
    if mode == "dotan2011":
        k = float(opacity(state.rho[-1],state.T[-1]))
        le = eddington(state.enclosed[-1],k)
        gamma = luminosity/le
        if gamma <= 1: return 0.
        # Isothermal GAS sound speed at the assumed sonic boundary.
        cs = (eos.R*state.T[-1])**.5
        raw = config["W"]*(luminosity-le)/(C*cs)
        tiring = luminosity*state.r[-1]/(G*state.enclosed[-1])
        return min(raw, min(.2*gamma**.6,.9)*tiring)
    raise ValueError(f"Unknown wind {mode}")
