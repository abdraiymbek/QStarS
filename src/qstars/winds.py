"""Wind launch closure and energy demand, separate from supplied material."""
from .constants import C,G,MSUN,YEAR
from .engine import eddington
from .controls import value

def rate(state,eos,opacity,config,luminosity):
    mode = config["model"]
    if mode == "off": return 0.
    if mode == "prescribed":
        return value(config["rate_msun_per_year"],state.time/YEAR)*MSUN/YEAR
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
