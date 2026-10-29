"""Independent whole-domain diagnostics, with signed non-rest-mass energies."""
from .constants import C

FIELDS = ["supplied_mass","accreted_rest_mass","wind_mass","mass_defect",
          "supplied_energy","wind_energy","advected_energy","central_generated",
          "radiated_energy","boundary_work","central_potential_work",
          "supplied_entropy","wind_entropy","advected_entropy",
          "mixing_entropy_change","remap_entropy_change","wind_mechanical_energy","numerical_kinetic_heating"]

# Moving-cavity (Ball 2011) channels. Present only in cavity runs and read
# with .get(k,0.), so older ledgers and checkpoints keep their exact keys.
# Mass identity: swept gas stays inside bh_mass+cavity_mass, so
# cavity_swept_mass is a flow diagnostic, not a term of the mass residual.
CAVITY_FIELDS = ["cavity_swept_mass","cavity_swept_energy","cavity_swept_entropy",
                 "inner_boundary_work","remap_energy_defect","relaxation_injected_energy"]

def initialize(state,eos,cavity=False):
    state.ledger = {k:0. for k in FIELDS}
    if cavity: state.ledger.update({k:0. for k in CAVITY_FIELDS})
    state.ledger.update(initial_mass=state.bh_mass+state.cavity_mass+state.dm.sum(),
                        initial_energy=state.energies(eos)["total"],
                        initial_energy_scale=state.energies(eos)["internal"])

def residuals(state,eos):
    l = state.ledger
    mass = state.bh_mass+state.cavity_mass+state.dm.sum()+l["wind_mass"]+l["mass_defect"]-l["supplied_mass"]-l["initial_mass"]
    e = state.energies(eos)
    energy = (e["total"]-l["initial_energy"]+l["radiated_energy"]+l["wind_energy"]+
              l["advected_energy"]-l["supplied_energy"]-l["central_generated"]-
              l["boundary_work"]-l["central_potential_work"])
    cavity = any(k in l for k in CAVITY_FIELDS)
    if cavity:
        energy += (l.get("cavity_swept_energy",0.)-l.get("relaxation_injected_energy",0.)-
                   l.get("inner_boundary_work",0.)-l.get("remap_energy_defect",0.))
    # Separate closure on rest-mass conversion; not fitted to energy residual.
    engine = l["mass_defect"]*C**2-l["central_generated"]
    throughput = sum(abs(l[k]) for k in ("radiated_energy","wind_energy","advected_energy",
                    "supplied_energy","central_generated","boundary_work","central_potential_work"))
    if cavity:
        throughput += sum(abs(l.get(k,0.)) for k in ("cavity_swept_energy","relaxation_injected_energy",
                          "inner_boundary_work","remap_energy_defect"))
    return dict(mass_residual_g=float(mass),mass_residual_relative=float(mass/l["initial_mass"]),
                energy_residual_erg=float(energy),
                energy_residual_relative=float(energy/l["initial_energy_scale"]),
                energy_residual_relative_to_binding=float(energy/max(abs(l["initial_energy"]),1.)),
                energy_residual_relative_to_throughput=float(energy/max(throughput,1.)) if throughput else 0.,
                engine_conversion_residual_erg=float(engine),**e)
