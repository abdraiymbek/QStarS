"""Strict, versioned YAML configuration. Every resolved default is saved."""
import copy
from pathlib import Path
import math
import re
import yaml

class ScientificLoader(yaml.SafeLoader):
    """Safe YAML with scientific notation and only true/false booleans."""
    yaml_implicit_resolvers=copy.deepcopy(yaml.SafeLoader.yaml_implicit_resolvers)

for _char,_resolvers in ScientificLoader.yaml_implicit_resolvers.items():
    ScientificLoader.yaml_implicit_resolvers[_char]=[(tag,rx) for tag,rx in _resolvers
                                                    if tag!='tag:yaml.org,2002:bool']
ScientificLoader.add_implicit_resolver('tag:yaml.org,2002:bool',re.compile(r'^(?:true|false|True|False|TRUE|FALSE)$'),list('tTfF'))
ScientificLoader.add_implicit_resolver('tag:yaml.org,2002:float',
    re.compile(r'^[-+]?(?:[0-9][0-9_]*(?:\.[0-9_]*)?|\.[0-9_]+)[eE][-+]?[0-9]+$'),list('-+0123456789.'))

DEFAULT = dict(
    schema_version=1,name="hot_envelope",fidelity="exploratory",
    initial=dict(envelope_mass_msun=1e4,bh_mass_msun=100.,radius_cm=1e14,inner_radius_cm=1e12,
                 temperature_perturbation=0.,perturbation_mass_center=.5,perturbation_mass_width=.1),
    composition=dict(X=.75,Y=.25,Z=0.),
    eos="ideal_gas_radiation",
    mesa=dict(library_path="",mesa_dir=""),
    opacity=dict(model="electron_scattering",table_path=""),
    engine=dict(model="bra2008",efficiency=.1,coupling=.1,rest_mass_rate_msun_per_year=1e-5),
    feeding=dict(rate_msun_per_year=0.,thermodynamics="specified",density_cgs=1e-7,temperature_K=1e6),
    wind=dict(model="off",rate_msun_per_year=0.,W=10.,terminal_speed_cm_s=0.,node_path="",drift_tolerance=0.),
    transport=dict(convection="mlt_lossless",mixing_length_alpha=2.),
    boundary=dict(inner="fixed_excision",pressure_cgs=0.,thermal="grey",external_optical_depth=0.,
                  bath_temperature_K=0.,luminosity_erg_s=1e41),
    evolution=dict(end_time_years=.1,max_steps=10000,max_wall_seconds=0.),
    numerics=dict(zones=32,mesh="cosine",initial_dt_years=.001,max_dt_years=.005,min_dt_years=1e-10,
                  nonlinear_solver="newton_flux",mechanics="midpoint",boundary_remesh_fraction=.2,
                  solver_tolerance=1e-9,max_fractional_change=.02,max_retries=12,
                  conservation_tolerance=1e-7,refine_every_steps=0,max_zones=256),
    validity=dict(min_temperature_K=1e5,max_temperature_K=1e8,max_mach=.2,max_compactness=.001,
                  require_optically_thick=True),
    output=dict(profile_every_steps=10,checkpoint_every_steps=10))

PRESETS = {"exploratory":{},"research":dict(zones=96,solver_tolerance=1e-10,max_dt_years=.001,
           max_fractional_change=.01),"validation":dict(zones=192,solver_tolerance=1e-11,
           max_dt_years=.0005,max_fractional_change=.005)}

# Quasi-static mechanics (Ball et al. 2011 replication; quasistatic.py). Its sections and
# opacity.tables are absent from DEFAULT, so every other configuration resolves, and is
# saved, exactly as before; a quasi_static configuration gets all of them with every
# default (quasistatic.DEFAULTS, RELAX_DEFAULTS, REPLICATION_DEFAULTS, stepcontrol.DEFAULTS).
QUASI_STATIC_SECTIONS = ("cavity","quasistatic","stepcontrol","relaxation","replication")
# default: the opacity of opacity.model (for mesa, MESA's default handle); stars2004_switch:
# the ABI-7 STARS switch with the nearest-table-value rule off the table (published-assumption
# replication only); stars2004_switch_strict: the same switch, off-table states are errors.
OPACITY_TABLES = ("default","stars2004_switch","stars2004_switch_strict")
# Ball replication closures exist only on the quasi-static path: the evolution path's legacy
# screens apply the lossless-MLT Pe screen, and it has no implicit photosphere-density row.
QUASI_STATIC_ONLY = [("boundary","thermal","stars_photosphere"),("transport","convection","mlt_ml1"),
                     ("transport","convection","mlt_thesis_lossless")]
# Keys the quasi-static path never reads: a value other than the (preset) default is an error,
# not a silent no-op. Its tracks end on black-hole mass or Ball's termination trigger.
QUASI_STATIC_UNUSED = dict(initial=list(DEFAULT["initial"]),evolution=["end_time_years"],output=["profile_every_steps"],
    numerics=["mesh","initial_dt_years","max_dt_years","min_dt_years","solver_tolerance","max_fractional_change",
              "max_retries","boundary_remesh_fraction","refine_every_steps","max_zones"])
# The Ball photosphere sets the outer pressure itself (P_g + a Teff^4/3 at tau = 2/3, no bath).
PHOTOSPHERE_UNUSED = ["pressure_cgs","external_optical_depth","bath_temperature_K"]
# Read at the step's time without the evolution path's control breakpoints: numbers only.
QUASI_STATIC_CONSTANT = [("engine","coupling"),("engine","rest_mass_rate_msun_per_year"),("boundary","luminosity_erg_s")]

def quasi_static_defaults():
    from .quasistatic import DEFAULTS,RELAX_DEFAULTS,REPLICATION_DEFAULTS
    from .stepcontrol import DEFAULTS as STEP
    return copy.deepcopy(dict(cavity=DEFAULTS["cavity"],quasistatic=DEFAULTS["quasistatic"],stepcontrol=STEP,
                              relaxation=RELAX_DEFAULTS,replication=REPLICATION_DEFAULTS))

def merge(base,updates,path=""):
    for k,v in updates.items():
        if k not in base: raise ValueError(f"Unknown configuration key {path+k}")
        if isinstance(base[k],dict) and isinstance(v,dict): merge(base[k],v,path+k+".")
        else: base[k] = copy.deepcopy(v)
    return base

def check_control(v,label,nonnegative=True,controls=True):
    from .controls import value
    if not isinstance(v,(float,int,dict)) or isinstance(v,bool): raise ValueError(f'{label}: expected number or control mapping')
    # Only stream rates and surface luminosity are time-evaluated (with breakpoints).
    if isinstance(v,dict) and not controls: raise ValueError(f'{label}: expected a number; time controls are not evaluated here')
    if isinstance(v,dict):
        kind = v.get("kind")
        schemas = {"table":{"kind","time_years","values"},
                   "pulse":{"kind","base","amplitude","start_years","end_years"},
                   "sinusoid":{"kind","base","amplitude","period_years","phase"}}
        if kind not in schemas or set(v)-schemas[kind]: raise ValueError(f"Invalid {label} control")
        if kind == "table":
            t,x = v["time_years"],v["values"]
            if len(t)!=len(x) or len(t)<2 or any(b<=a for a,b in zip(t,t[1:])):
                raise ValueError(f"{label}: table times must increase and match values")
            samples = list(t)
        elif kind == "pulse":
            if v["end_years"] <= v["start_years"]: raise ValueError("Pulse must have positive duration")
            samples = [v["start_years"]-1,v["start_years"]]
        else:
            if v["period_years"] <= 0: raise ValueError("Sinusoid period must be positive")
            if nonnegative and v["base"] < abs(v["amplitude"]): raise ValueError(f"{label}: negative sinusoid")
            samples = [0.]
    else: samples = [0.]
    for t in samples:
        x = value(v,t)
        if not math.isfinite(x) or (nonnegative and x < 0): raise ValueError(f"Invalid {label}")

def resolve(raw):
    raw = raw or {}
    if not isinstance(raw,dict): raise ValueError("Configuration must be a mapping")
    c = copy.deepcopy(DEFAULT)
    mode = raw.get("fidelity",c["fidelity"])
    if mode not in PRESETS: raise ValueError("Unknown fidelity mode")
    c["numerics"].update(PRESETS[mode])
    numerics,opacity = raw.get("numerics"),raw.get("opacity")
    quasi = isinstance(numerics,dict) and numerics.get("mechanics")=="quasi_static"
    extra = [k for k in QUASI_STATIC_SECTIONS if k in raw]+(["opacity.tables"] if isinstance(opacity,dict) and "tables" in opacity else [])
    if extra and not quasi: raise ValueError(f"{', '.join(extra)} require numerics.mechanics quasi_static")
    if quasi:
        defaults = copy.deepcopy(c)
        c.update(quasi_static_defaults()); c["opacity"]["tables"] = "default"
        for k in QUASI_STATIC_SECTIONS:
            if not isinstance(raw.get(k,{}),dict): raise ValueError(f"{k} must be a mapping")
    merge(c,raw)
    if not isinstance(c['name'],str) or not re.fullmatch(r'[a-zA-Z0-9_-]{1,100}',c['name']):
        raise ValueError('name must use 1-100 letters, digits, hyphens, or underscores')
    if c["schema_version"] != 1: raise ValueError("Unsupported schema version")
    choices = [(c["eos"],["ideal_gas_radiation","mesa"]),
        (c['numerics']['mesh'],['cosine','core_refined','log_radius','log_radius_cosine','log_radius_surface']),
        (c['numerics']['nonlinear_solver'],['newton','newton_flux','hybr']),
        (c['numerics']['mechanics'],['midpoint','backward_euler','quasi_static']),
        (c["opacity"]["model"],["electron_scattering","table","mesa"]),
        (c["boundary"]["inner"],["fixed_excision","bondi_cavity"]),
        (c["boundary"]["thermal"],["grey","luminosity","insulated","stars_photosphere"]),
        (c["engine"]["model"],["off","bra2008","ball2011_printed","quest2025_eq11","ball2012_eq3_15","prescribed"]),
        (c["wind"]["model"],["off","prescribed","dotan2011","selected_nodes"]),
        (c["feeding"]["thermodynamics"],["specified","surface"]),
        (c["transport"]["convection"],["off","mlt_lossless","mlt_ml1","mlt_thesis_lossless"])]
    if quasi: choices.append((c["opacity"]["tables"],list(OPACITY_TABLES)))
    for v,opts in choices:
        if v not in opts: raise ValueError(f"Unsupported choice {v}; choose from {opts}")
    for section,key,choice in QUASI_STATIC_ONLY:
        if c[section][key]==choice and not quasi: raise ValueError(f"{section}.{key} {choice} requires numerics.mechanics quasi_static")
    if quasi or c["boundary"]["inner"]=="bondi_cavity":
        # The moving Bondi cavity (thesis 3.2-3.15) is the inner boundary of quasistatic.bondi_step only;
        # its swept-gas bookings assume no surface streams and the flux-unknown solve.
        if not (quasi and c["boundary"]["inner"]=="bondi_cavity"):
            raise ValueError("The Bondi cavity and quasi_static mechanics go together: set boundary.inner bondi_cavity and numerics.mechanics quasi_static")
        if c["numerics"]["nonlinear_solver"]!="newton_flux": raise ValueError("The Bondi cavity requires numerics.nonlinear_solver newton_flux")
        if c["wind"]["model"]!="off": raise ValueError("The Bondi cavity requires wind.model off")
        rate=c["feeding"]["rate_msun_per_year"]
        if isinstance(rate,(dict,bool)) or rate!=0: raise ValueError("The Bondi cavity requires zero feeding (feeding.rate_msun_per_year 0)")
        if c["engine"]["model"]=="off": raise ValueError("quasi_static mechanics needs a central engine: L_BH sets the relaxed luminosity")
    def positive(section,keys):
        for k in keys:
            x=c[section][k]
            if isinstance(x,bool) or not isinstance(x,(float,int)) or not math.isfinite(x) or x<=0:
                raise ValueError(f"{section}.{k} must be positive and finite")
    positive("initial",["envelope_mass_msun","bh_mass_msun","radius_cm","inner_radius_cm","perturbation_mass_width"])
    positive("numerics",["zones","initial_dt_years","max_dt_years","min_dt_years","solver_tolerance","max_fractional_change","max_retries","conservation_tolerance","max_zones"])
    positive("evolution",["end_time_years","max_steps"])
    wall=c['evolution']['max_wall_seconds']
    if isinstance(wall,bool) or not isinstance(wall,(int,float)) or not math.isfinite(wall) or wall<0:
        raise ValueError('evolution.max_wall_seconds must be a finite nonnegative number')
    positive("output",["profile_every_steps","checkpoint_every_steps"])
    positive("validity",["min_temperature_K","max_temperature_K","max_mach","max_compactness"])
    positive("feeding",["density_cgs"])
    positive("transport",["mixing_length_alpha"])
    for key in ["zones","max_retries","max_zones","refine_every_steps"]:
        if not isinstance(c["numerics"][key],int): raise ValueError(f"{key} must be an integer")
    if c['numerics']['refine_every_steps']<0:raise ValueError('refine_every_steps must be nonnegative')
    if not 0<=c['numerics']['boundary_remesh_fraction']<=.5:
        raise ValueError('boundary_remesh_fraction must lie between 0 (disabled) and 0.5')
    for key in ('profile_every_steps','checkpoint_every_steps'):
        if not isinstance(c['output'][key],int):raise ValueError(f'{key} must be an integer')
    if c["numerics"]["zones"]<4: raise ValueError("At least four zones required")
    floor=3000. if c['eos']=='mesa' else 1e5
    if c['validity']['min_temperature_K']<floor or c['validity']['max_temperature_K']>1e8:
        raise ValueError(f'Configured EOS limits must stay within the supported {floor:g}-1e8 K interval')
    if (c['eos']=='mesa') != (c['opacity']['model']=='mesa'):
        raise ValueError('MESA EOS and opacity must be selected together for this audited interface')
    if c['eos']=='mesa' and (not c['mesa']['library_path'] or not c['mesa']['mesa_dir']):
        raise ValueError('MESA requires explicit library_path and mesa_dir')
    if c["numerics"]["min_dt_years"]>c["numerics"]["max_dt_years"]: raise ValueError("Invalid timestep bounds")
    if not 0<c["engine"]["efficiency"]<1: raise ValueError("Efficiency must lie strictly between zero and one")
    if c["initial"]["radius_cm"]<=c["initial"]["inner_radius_cm"]: raise ValueError("Outer radius must exceed inner radius")
    if c["initial"]["temperature_perturbation"]<=-1: raise ValueError("Temperature perturbation cannot remove all heat")
    for section,key in [("engine","coupling"),("engine","rest_mass_rate_msun_per_year"),
        ("feeding","rate_msun_per_year"),("feeding","temperature_K"),("wind","rate_msun_per_year"),
        ("boundary","luminosity_erg_s")]: check_control(c[section][key],section+"."+key)
    for section,key in [("boundary","pressure_cgs"),("boundary","external_optical_depth"),
        ("boundary","bath_temperature_K"),("wind","terminal_speed_cm_s"),("wind","W")]:
        check_control(c[section][key],section+"."+key,controls=False)
    w=c['wind'];d=w['drift_tolerance']
    if not isinstance(w['node_path'],str): raise ValueError('wind.node_path must be a string')
    if isinstance(d,bool) or not isinstance(d,(int,float)) or not math.isfinite(d) or not 0<=d<1:
        raise ValueError('wind.drift_tolerance must lie in [0, 1)')
    if w['model']=='selected_nodes':
        if not w['node_path'] or d<=0: raise ValueError('selected_nodes wind requires node_path and a positive drift_tolerance')
        if w['terminal_speed_cm_s']!=0 or w['rate_msun_per_year']!=0:
            raise ValueError('selected_nodes wind takes mdot and terminal speed from the node; leave rate_msun_per_year and terminal_speed_cm_s at 0')
    if quasi: _quasi_static_checks(c,defaults)
    return c

def _quasi_static_checks(c,defaults):
    """Cross-validation of a quasi_static configuration; its sections hold every default."""
    from .quasistatic import settings,SEED_SURFACES
    from .photosphere import BRACKETS
    from .stepcontrol import settings as step_settings
    from .validity import POLICIES
    settings(c); step_settings(c)  # choices exactly as the step and the controller read them
    real=lambda x:not isinstance(x,bool) and isinstance(x,(int,float)) and math.isfinite(x)
    count=lambda x:isinstance(x,int) and not isinstance(x,bool)
    positive=lambda x:real(x) and x>0
    fraction=lambda x:real(x) and 0<x<1
    factor=lambda x:real(x) and 0<x<=1
    text=lambda x:isinstance(x,str)
    q,s,r,p=c["quasistatic"],c["stepcontrol"],c["relaxation"],c["replication"]
    rules=dict(
        quasistatic=dict(tolerance=positive,max_iterations=lambda x:count(x) and x>=1,polish_iterations=lambda x:count(x) and x>=0,
                         roundoff_residual=positive,roundoff_step=positive,max_passes=lambda x:count(x) and x>=1,pass_fraction=fraction,
                         remesh=lambda x:isinstance(x,bool),ptc=lambda x:isinstance(x,bool),
                         photosphere_bracket=lambda x:isinstance(x,str) and x in BRACKETS,split=lambda x:real(x) and x>1,merge=fraction,
                         energy_per_step=positive,mass_per_step=positive),
        cavity=dict(active=lambda x:isinstance(x,bool),b=positive),
        stepcontrol=dict(delta_target=positive,safety=factor,min_factor=factor,max_factor=lambda x:real(x) and x>=1,
                         reject_ratio=lambda x:real(x) and x>1,hold_steps=lambda x:count(x) and x>=0,
                         slow_iterations=lambda x:count(x) and x>=0,slow_factor=factor,retry_factor=fraction,
                         max_retries=lambda x:count(x) and x>=1,dt_initial_years=positive,
                         dt_max_years=lambda x:real(x) and x>=s["dt_initial_years"]),
        relaxation=dict(dt_years=positive,dt_max_years=lambda x:real(x) and x>=r["dt_years"],
                        converged_dt_years=lambda x:positive(x) and x<=r["dt_max_years"],tolerance=positive,
                        luminosity_tolerance=positive,max_steps=lambda x:count(x) and x>=1,max_retries=lambda x:count(x) and x>=0,
                        mesh_cycles=lambda x:count(x) and x>=0,mesh_tolerance=positive,cycle_dt_years=positive),
        replication=dict(plan=text,run=text,kind=text,plan_text=text,policy=lambda x:isinstance(x,str) and x in POLICIES,
                         total_mass_msun=positive,initial_bh_mass_msun=lambda x:positive(x) and x<p["total_mass_msun"],
                         seed_radii_rsun=lambda x:isinstance(x,list) and len(x)>0 and all(map(positive,x)),
                         seed_surface=lambda x:isinstance(x,str) and x in SEED_SURFACES,
                         landings_msun=lambda x:isinstance(x,list) and all(map(real,x)) and all(b>a for a,b in zip([p["initial_bh_mass_msun"]]+x,x)),
                         end_bh_mass_msun=lambda x:real(x) and (x==0 or x>p["initial_bh_mass_msun"]),
                         amendments=lambda x:isinstance(x,list) and all(map(text,x))))
    for section,checks in rules.items():
        for key,ok in checks.items():
            if not ok(c[section][key]): raise ValueError(f"Invalid {section}.{key}: {c[section][key]!r}")
    if q["merge"]*2>q["split"]: raise ValueError("quasistatic.merge must not exceed half of quasistatic.split (split/merge hysteresis)")
    if not c["cavity"]["active"]:
        raise ValueError("cavity.active false is the frozen-cavity relaxation of the identity tests; a bondi_cavity configuration moves the cavity")
    if c["opacity"]["model"]=="table": raise ValueError("quasi_static mechanics takes opacity.model electron_scattering or mesa")
    for section,key in QUASI_STATIC_CONSTANT: check_control(c[section][key],section+"."+key,controls=False)
    if c["boundary"]["thermal"]=="stars_photosphere":
        for key in PHOTOSPHERE_UNUSED:
            if c["boundary"][key]!=0: raise ValueError(f"boundary.{key} is not read with boundary.thermal stars_photosphere; leave it at 0")
    if p["seed_surface"]=="grey" and c["boundary"]["thermal"]!="stars_photosphere":
        raise ValueError("replication.seed_surface grey is the grey surface relation of the Ball photosphere; it needs boundary.thermal stars_photosphere")
    tables=c["opacity"]["tables"]
    if tables!="default" and c["opacity"]["model"]!="mesa": raise ValueError(f"opacity.tables {tables} is the MESA ABI-7 STARS switch; it needs opacity.model mesa")
    if tables=="stars2004_switch" and p["policy"]!="published_assumption":
        raise ValueError("The STARS nearest-table-value rule (opacity.tables stars2004_switch) is a published-assumption replication rule; production configurations use stars2004_switch_strict")
    for section,keys in QUASI_STATIC_UNUSED.items():
        for key in keys:
            if c[section][key]!=defaults[section][key]:
                raise ValueError(f"{section}.{key} is not read by quasi_static mechanics (seed and track end: replication; steps: stepcontrol)")

def load(path):
    with open(path) as f: config=resolve(yaml.load(f,Loader=ScientificLoader))
    if config['opacity']['model']=='table':
        p=Path(config['opacity']['table_path'])
        if not p.is_absolute():p=Path(path).resolve().parent/p
        config['opacity']['table_path']=str(p.resolve())
    if config['wind']['model']=='selected_nodes':
        p=Path(config['wind']['node_path'])
        if not p.is_absolute():p=Path(path).resolve().parent/p
        config['wind']['node_path']=str(p.resolve())
    if config['eos']=='mesa':
        for key in ('library_path','mesa_dir'):
            p=Path(config['mesa'][key])
            if not p.is_absolute():p=Path(path).resolve().parent/p
            config['mesa'][key]=str(p.resolve())
    return config
