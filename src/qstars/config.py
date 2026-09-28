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
    opacity=dict(model="electron_scattering",table_path=""),
    engine=dict(model="bra2008",efficiency=.1,coupling=.1,rest_mass_rate_msun_per_year=1e-5),
    feeding=dict(rate_msun_per_year=0.,thermodynamics="specified",density_cgs=1e-7,temperature_K=1e6),
    wind=dict(model="off",rate_msun_per_year=0.,W=10.,terminal_speed_cm_s=0.),
    transport=dict(convection="mlt_lossless",mixing_length_alpha=2.),
    boundary=dict(inner="fixed_excision",pressure_cgs=0.,thermal="grey",external_optical_depth=0.,
                  bath_temperature_K=0.,luminosity_erg_s=1e41),
    evolution=dict(end_time_years=.1,max_steps=10000,max_wall_seconds=0.),
    numerics=dict(zones=32,mesh="cosine",initial_dt_years=.001,max_dt_years=.005,min_dt_years=1e-10,
                  solver_tolerance=1e-9,max_fractional_change=.02,max_retries=12,
                  conservation_tolerance=1e-7,refine_every_steps=0,max_zones=256),
    validity=dict(min_temperature_K=1e5,max_temperature_K=1e8,max_mach=.2,max_compactness=.001,
                  require_optically_thick=True),
    output=dict(profile_every_steps=10,checkpoint_every_steps=10))

PRESETS = {"exploratory":{},"research":dict(zones=96,solver_tolerance=1e-10,max_dt_years=.001,
           max_fractional_change=.01),"validation":dict(zones=192,solver_tolerance=1e-11,
           max_dt_years=.0005,max_fractional_change=.005)}

def merge(base,updates,path=""):
    for k,v in updates.items():
        if k not in base: raise ValueError(f"Unknown configuration key {path+k}")
        if isinstance(base[k],dict) and isinstance(v,dict): merge(base[k],v,path+k+".")
        else: base[k] = copy.deepcopy(v)
    return base

def check_control(v,label,nonnegative=True):
    from .controls import value
    if not isinstance(v,(float,int,dict)) or isinstance(v,bool): raise ValueError(f'{label}: expected number or control mapping')
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
    c["numerics"].update(PRESETS[mode]); merge(c,raw)
    if not isinstance(c['name'],str) or not re.fullmatch(r'[a-zA-Z0-9_-]{1,100}',c['name']):
        raise ValueError('name must use 1-100 letters, digits, hyphens, or underscores')
    if c["schema_version"] != 1: raise ValueError("Unsupported schema version")
    choices = [(c["eos"],["ideal_gas_radiation"]),
        (c['numerics']['mesh'],['cosine','core_refined']),
        (c["opacity"]["model"],["electron_scattering","table"]),
        (c["boundary"]["inner"],["fixed_excision"]),
        (c["boundary"]["thermal"],["grey","luminosity","insulated"]),
        (c["engine"]["model"],["off","bra2008","ball2011_printed","quest2025_eq11","prescribed"]),
        (c["wind"]["model"],["off","prescribed","dotan2011"]),
        (c["feeding"]["thermodynamics"],["specified","surface"]),
        (c["transport"]["convection"],["off","mlt_lossless"])]
    for v,opts in choices:
        if v not in opts: raise ValueError(f"Unsupported choice {v}; choose from {opts}")
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
    for key in ('profile_every_steps','checkpoint_every_steps'):
        if not isinstance(c['output'][key],int):raise ValueError(f'{key} must be an integer')
    if c["numerics"]["zones"]<4: raise ValueError("At least four zones required")
    if c['validity']['min_temperature_K']<1e5 or c['validity']['max_temperature_K']>1e8:
        raise ValueError('Configured EOS limits must stay within the supported 1e5-1e8 K interval')
    if c["numerics"]["min_dt_years"]>c["numerics"]["max_dt_years"]: raise ValueError("Invalid timestep bounds")
    if not 0<c["engine"]["efficiency"]<1: raise ValueError("Efficiency must lie strictly between zero and one")
    if c["initial"]["radius_cm"]<=c["initial"]["inner_radius_cm"]: raise ValueError("Outer radius must exceed inner radius")
    if c["initial"]["temperature_perturbation"]<=-1: raise ValueError("Temperature perturbation cannot remove all heat")
    for section,key in [("engine","coupling"),("engine","rest_mass_rate_msun_per_year"),
        ("feeding","rate_msun_per_year"),("feeding","temperature_K"),("wind","rate_msun_per_year"),
        ("boundary","luminosity_erg_s")]: check_control(c[section][key],section+"."+key)
    for section,key in [("boundary","pressure_cgs"),("boundary","external_optical_depth"),
        ("boundary","bath_temperature_K"),("wind","terminal_speed_cm_s"),("wind","W")]:
        check_control(c[section][key],section+"."+key)
    return c

def load(path):
    with open(path) as f: config=resolve(yaml.load(f,Loader=ScientificLoader))
    if config['opacity']['model']=='table':
        p=Path(config['opacity']['table_path'])
        if not p.is_absolute():p=Path(path).resolve().parent/p
        config['opacity']['table_path']=str(p.resolve())
    return config
