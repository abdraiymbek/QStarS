"""Rest-mass rates are always distinct from black-hole gravitational growth."""
from .constants import G, C, PI, YEAR, MSUN
from .controls import value

def eddington(mass,kappa):
    return 4*PI*G*mass*C/kappa

def rates(state,eos,config):
    q = eos.evaluate(state.rho,state.T)
    rho,p,cs,gamma = (float(q[k][0]) if k != "rho" else float(state.rho[0])
                        for k in ("rho","p","cs","gamma1"))
    mode = config["model"]
    eps = config["efficiency"]
    eta = value(config["coupling"],state.time/YEAR)
    if mode == "off": luminosity = 0.
    elif mode == "bra2008":
        luminosity = 4*PI*G**2*eta*state.bh_mass**2*rho**1.5/p**.5
    elif mode in ("ball2011_printed", "quest2025_eq11"):
        prefactor = 4 if mode == "ball2011_printed" else 16
        luminosity = prefactor*PI*eta*(G*state.bh_mass)**2*rho/(gamma*cs)
    elif mode == "prescribed":
        luminosity = eps*value(config["rest_mass_rate_msun_per_year"],state.time/YEAR)*MSUN/YEAR*C**2
    else: raise ValueError(f"Unknown BH engine {mode}")
    if luminosity < 0: raise ValueError("BH luminosity must be nonnegative")
    rest = luminosity/(eps*C**2)
    return dict(luminosity=luminosity,rest_mass_rate=rest,bh_growth_rate=(1-eps)*rest,
                bondi_radius=2*G*state.bh_mass/cs**2)
