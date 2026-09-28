"""Independent whole-domain diagnostics, with signed non-rest-mass energies."""
from .constants import C

FIELDS = ["supplied_mass","accreted_rest_mass","wind_mass","mass_defect",
          "supplied_energy","wind_energy","advected_energy","central_generated",
          "radiated_energy","boundary_work","central_potential_work",
          "supplied_entropy","wind_entropy","advected_entropy",
          "mixing_entropy_change","remap_entropy_change","wind_mechanical_energy"]

def initialize(state,eos):
    state.ledger = {k:0. for k in FIELDS}
    state.ledger.update(initial_mass=state.bh_mass+state.dm.sum(),
                        initial_energy=state.energies(eos)["total"],
                        initial_energy_scale=state.energies(eos)["internal"])

def residuals(state,eos):
    l = state.ledger
    mass = state.bh_mass+state.dm.sum()+l["wind_mass"]+l["mass_defect"]-l["supplied_mass"]-l["initial_mass"]
    e = state.energies(eos)
    energy = (e["total"]-l["initial_energy"]+l["radiated_energy"]+l["wind_energy"]+
              l["advected_energy"]-l["supplied_energy"]-l["central_generated"]-
              l["boundary_work"]-l["central_potential_work"])
    # Separate closure on rest-mass conversion; not fitted to energy residual.
    engine = l["mass_defect"]*C**2-l["central_generated"]
    throughput = sum(abs(l[k]) for k in ("radiated_energy","wind_energy","advected_energy",
                    "supplied_energy","central_generated","boundary_work","central_potential_work"))
    return dict(mass_residual_g=float(mass),mass_residual_relative=float(mass/l["initial_mass"]),
                energy_residual_erg=float(energy),
                energy_residual_relative=float(energy/l["initial_energy_scale"]),
                energy_residual_relative_to_binding=float(energy/max(abs(l["initial_energy"]),1.)),
                energy_residual_relative_to_throughput=float(energy/max(throughput,1.)) if throughput else 0.,
                engine_conversion_residual_erg=float(engine),**e)
