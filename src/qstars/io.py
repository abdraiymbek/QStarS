"""Self-describing output, atomic checkpoints, and source snapshot provenance."""
import csv,hashlib,json,os,platform,subprocess,sys,zipfile
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import yaml
from . import __version__
from .state import State
from .constants import YEAR,MSUN,G,SIGMA,PI
from .transport import luminosities
from .engine import rates,eddington
from .controls import value
from .winds import evaluate as wind_evaluate
from .ledger import residuals
from .structure import hydrostatic_residual

def _package_sources():
    root=Path(__file__).parent
    return {str(p.relative_to(root.parent)):p.read_bytes() for p in sorted(root.rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts}

def _source_digest(sources):
    digest=hashlib.sha256()
    for name,data in sources.items():digest.update(name.encode()+b'\0'+data+b'\0')
    return digest.hexdigest()

# A long-running campaign imports Python modules once. Reading current disk
# files for each later run could otherwise archive code it never executed.
_PROCESS_SOURCES=_package_sources()
_PROCESS_SOURCE_DIGEST=_source_digest(_PROCESS_SOURCES)

def finish_provenance(path):
    file=Path(path)/'metadata.json';meta=json.loads(file.read_text())
    unchanged=_source_digest(_package_sources())==meta['process_package_sha256']
    verified=bool(unchanged and meta.get('source_tree_matches_import',False))
    meta.update(source_tree_matches_at_end=unchanged,source_snapshot_verified=verified)
    if not verified:meta['source_warning']='Package files changed during this process. Dynamic imports may differ from the frozen source snapshot; rerun in a fresh process with unchanged source.'
    json_write(file,meta)
    return verified

def json_write(path,data):
    Path(path).write_text(json.dumps(data,indent=2,allow_nan=False,default=lambda x:float(x))+"\n")

def create_directory(parent,name):
    stamp=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    path=Path(parent)/f"{stamp}_{name}"
    path.mkdir(parents=True,exist_ok=False)
    return path

def provenance(path,config):
    root=Path(__file__).resolve().parents[2]
    def git(*args):
        return subprocess.run(["git","-C",str(root),*args],capture_output=True,text=True).stdout.strip()
    snapshot=path/"source_snapshot.zip"
    with zipfile.ZipFile(snapshot,"w",compression=zipfile.ZIP_DEFLATED) as z:
        for name,data in _PROCESS_SOURCES.items():z.writestr(name,data)
    meta=dict(version=__version__,git_commit=git("rev-parse","HEAD"),git_status=git("status","--porcelain"),
              source_sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest(),python=sys.version,
              platform=platform.platform(),created_utc=datetime.now(timezone.utc).isoformat())
    meta['process_package_sha256']=_PROCESS_SOURCE_DIGEST
    meta['source_tree_matches_import']=_source_digest(_package_sources())==_PROCESS_SOURCE_DIGEST
    (path/"config.yaml").write_text(yaml.safe_dump(config,sort_keys=False))
    (path/"working_tree.diff").write_text(git("diff","HEAD"))
    import importlib.metadata
    meta["dependencies"]={p:importlib.metadata.version(p) for p in ("numpy","scipy","PyYAML","matplotlib")}
    json_write(path/"metadata.json",meta)
    registry=Path(__file__).parent/"data"/"sources.json"
    if registry.exists(): (path/"source_prescriptions.json").write_bytes(registry.read_bytes())

def save_checkpoint(path,state,config,schema=None):
    path=Path(path)
    # NPZ contains JSON strings, never pickled Python objects. Schema 7 adds the
    # quasi-static fields (cavity_mass, photosphere_density, kinetic_tracked, aux).
    # schema=None: a state that holds them all at their defaults (every
    # evolution-path state) is written as schema 6, byte for byte as before, and
    # schema 6 loads with exactly those defaults; any other state is schema 7.
    # schema=7 (quasistatic.save_checkpoint) always writes them.
    if schema not in (None,7): raise ValueError("Checkpoints are written as schema 7, or 6 for a state without quasi-static fields")
    extended=dict(cavity_mass=state.cavity_mass,
                  photosphere_density=np.nan if state.photosphere_density is None else state.photosphere_density,
                  kinetic_tracked=bool(state.kinetic_tracked),aux=json.dumps(state.aux))
    if schema is None and not (state.cavity_mass!=0 or state.photosphere_density is not None or not state.kinetic_tracked or state.aux):extended={}
    temporary=path.with_suffix(".tmp")
    with temporary.open("wb") as f:
        np.savez_compressed(f,schema_version=7 if extended else 6,time=state.time,bh_mass=state.bh_mass,dm=state.dm,r=state.r,v=state.v,T=state.T,
            cell_volume=state.volume,
            radius_compensation=state.radius_compensation if state.radius_compensation is not None else np.zeros_like(state.r),
            face_luminosity=state.face_luminosity if state.face_luminosity is not None else np.array([]),
            reference_mass_fractions=state.reference_mass_fractions if state.reference_mass_fractions is not None else state.dm/state.dm.sum(),
            ledger=json.dumps(state.ledger),steps=state.steps,next_dt=state.next_dt,config=json.dumps(config),**extended)
        f.flush(); os.fsync(f.fileno())
    temporary.replace(path)

def load_checkpoint(path,resolve_config=True):
    """(state, config) of a checkpoint of schema 1-7; older schemas get the defaults of later fields.

    resolve_config=False returns the stored configuration as saved
    (quasistatic.load_checkpoint). The Stage 5 quasi-static format
    ('qstars_quasistatic' 1) holds exactly the schema-7 fields.
    """
    from .config import resolve
    with np.load(path,allow_pickle=False) as f:
        schema=int(f["schema_version"])
        if 'format' in f.files:
            if str(f['format'])!='qstars_quasistatic' or schema!=1: raise ValueError("Unsupported checkpoint schema")
            schema=7
        if schema not in (1,2,3,4,5,6,7): raise ValueError("Unsupported checkpoint schema")
        state=State(float(f["time"]),float(f["bh_mass"]),f["dm"].copy(),f["r"].copy(),f["v"].copy(),f["T"].copy(),
                    json.loads(str(f["ledger"])),int(f["steps"]),float(f["next_dt"]))
        stored=json.loads(str(f["config"]))
        config=resolve(stored) if resolve_config else stored
        if schema >= 2: state.cell_volume=f['cell_volume'].copy()
        if schema >= 3 and len(f['face_luminosity']):state.face_luminosity=f['face_luminosity'].copy()
        state.reference_mass_fractions=f['reference_mass_fractions'].copy() if schema>=4 else state.dm/state.dm.sum()
        state.radius_compensation=f['radius_compensation'].copy() if schema>=5 else np.zeros_like(state.r)
        if schema<6:state.ledger.setdefault('numerical_kinetic_heating',0.)
        if schema>=7:
            state.cavity_mass=float(f['cavity_mass']);rho_ph=float(f['photosphere_density'])
            state.photosphere_density=None if np.isnan(rho_ph) else rho_ph
            state.kinetic_tracked=bool(f['kinetic_tracked']);state.aux=json.loads(str(f['aux']))
    if not (len(state.r)==len(state.dm)+1==len(state.v) and len(state.T)==len(state.dm)):
        raise ValueError("Malformed checkpoint dimensions")
    if any(not np.all(np.isfinite(x)) for x in (state.r,state.dm,state.v,state.T)):
        raise ValueError('Checkpoint contains nonfinite state arrays')
    # Quasi-static steps are controlled from aux (a t = 0 model has next_dt 0).
    if not (np.isfinite(state.time) and state.time>=0 and np.isfinite(state.bh_mass) and state.bh_mass>0
            and np.isfinite(state.next_dt) and (state.next_dt>0 or schema>=7 and state.next_dt==0) and state.steps>=0):
        raise ValueError('Checkpoint contains invalid time, central mass, or step controls')
    if not (np.isfinite(state.cavity_mass) and state.cavity_mass>=0 and isinstance(state.aux,dict) and
            (state.photosphere_density is None or np.isfinite(state.photosphere_density) and state.photosphere_density>0)):
        raise ValueError('Checkpoint contains invalid cavity, photosphere or driver state')
    from .ledger import FIELDS
    if not set(FIELDS+['initial_mass','initial_energy','initial_energy_scale'])<=state.ledger.keys():
        raise ValueError('Checkpoint ledger is incomplete')
    if any(isinstance(x,bool) or not isinstance(x,(int,float)) or not np.isfinite(x) for x in state.ledger.values()):
        raise ValueError('Checkpoint contains invalid ledger values')
    if np.any(np.diff(state.r)<=0) or np.any(state.dm<=0) or np.any(state.T<=0):
        raise ValueError("Invalid checkpoint physical state")
    geometric_volume=4*PI/3*np.diff(state.r)*(state.r[1:]**2+state.r[1:]*state.r[:-1]+state.r[:-1]**2)
    if state.volume.shape != state.dm.shape or np.any(state.volume<=0) or not np.allclose(state.volume,geometric_volume,rtol=1e-8,atol=0):
        raise ValueError('Checkpoint volume and geometry disagree')
    if state.face_luminosity is not None and (state.face_luminosity.shape!=state.r.shape or not np.all(np.isfinite(state.face_luminosity))):
        raise ValueError('Malformed checkpoint face luminosities')
    if state.radius_compensation.shape!=state.r.shape or not np.all(np.isfinite(state.radius_compensation)):
        raise ValueError('Malformed checkpoint radius compensation')
    f=state.reference_mass_fractions
    if f.shape!=state.dm.shape or np.any(f<=0) or not np.isclose(f.sum(),1.,rtol=1e-12,atol=0):
        raise ValueError('Malformed checkpoint reference mass fractions')
    return state,config

def history_row(state,eos,opacity,config):
    if config['numerics']['mechanics']=='quasi_static':
        # Ball history columns: Table 1 observables, base state, cavity, photospheric
        # radius (the tau = 2/3 face of the Ball photosphere), screens and the ledger.
        from .quasistatic import history_row as ball_row
        return ball_row(state,eos,opacity,config)
    e=rates(state,eos,config["engine"])
    L,d=luminosities(state,eos,opacity,config,e["luminosity"])
    radius=float(state.r[-1]); total=float(state.enclosed[-1])
    supply=value(config['feeding']['rate_msun_per_year'],state.time/YEAR)
    w=wind_evaluate(state,eos,opacity,config['wind'],max(0.,float(L[-1])));wind=w['rate']
    er=state.energies(eos)
    transport=supply+wind*YEAR/MSUN+e['rest_mass_rate']*YEAR/MSUN
    row=dict(step=state.steps,time_years=state.time/YEAR,bh_mass_msun=state.bh_mass/MSUN,
        envelope_mass_msun=float(state.dm.sum()/MSUN),radius_cm=radius,
        luminosity_erg_s=float(L[-1]),engine_luminosity_erg_s=e["luminosity"],
        effective_temperature_diagnostic_K=float((max(L[-1],0)/(4*PI*radius**2*SIGMA))**.25),
        rest_accretion_msun_per_year=e["rest_mass_rate"]*YEAR/MSUN,
        bh_growth_msun_per_year=e["bh_growth_rate"]*YEAR/MSUN,
        supply_msun_per_year=supply,wind_msun_per_year=wind*YEAR/MSUN,
        photospheric_radius_cm=float('nan'),  # undefined; a node's held R_phot is not this state's photosphere
        kelvin_helmholtz_time_years=abs(er['gravity'])/max(abs(float(L[-1])),1.)/YEAR,
        replacement_time_years=float(state.dm.sum()/MSUN)/transport if transport else float('inf'),
        dynamical_time_years=(radius**3/(G*total))**.5/YEAR,
        eddington_ratio=float(L[-1]/eddington(total,d["opacity"][-1])),
        max_hydrostatic_residual=float(max(abs(hydrostatic_residual(state,eos,config["boundary"]["pressure_cgs"])))),
        max_convective_saturation=float(max(d["saturation"],default=0.)),
        envelope_entropy_erg_K=float(np.dot(state.dm,eos.evaluate(state.rho,state.T)['s'])),
        **residuals(state,eos))
    if 'emergent_luminosity' in w:
        row.update(wind_terminal_speed_cm_s=w['terminal_speed'],wind_emergent_luminosity_erg_s=w['emergent_luminosity'],
                   wind_emergent_effective_temperature_K=w['emergent_effective_temperature'],
                   wind_node_photospheric_radius_cm=w['photospheric_radius'],
                   wind_node_max_abs_dln_drift=w['max_abs_dln_drift'])
    row.update({"ledger_"+k:float(v) for k,v in state.ledger.items()})
    return row

def append_history(path,row):
    exists=Path(path).exists()
    with open(path,"a",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(row))
        if not exists: writer.writeheader()
        writer.writerow(row)

def save_profile(path,state,eos,opacity,config):
    q=eos.evaluate(state.rho,state.T)
    e=rates(state,eos,config["engine"])
    L,d=luminosities(state,eos,opacity,config,e["luminosity"])
    if state.face_luminosity is not None:
        L=state.face_luminosity
        d['convective']=L-d['radiative'];d['convective'][0]=0.;d['convective'][-1]=0.
    np.savez_compressed(path,time_seconds=state.time,dm_g=state.dm,r_edges_cm=state.r,r_cm=state.rc,
        enclosed_mass_g=.5*(state.enclosed[:-1]+state.enclosed[1:]),density_cgs=state.rho,
        pressure_cgs=q["p"],temperature_K=state.T,entropy_cgs=q["s"],internal_energy_erg_g=q["u"],
        opacity_cgs=d["opacity"],luminosity_edges_erg_s=L,radiative_edges_erg_s=d["radiative"],
        convective_edges_erg_s=d["convective"],grad_ad=q["grad_ad"],velocity_edges_cm_s=state.v,
        eddington_ratio_edges=d["radiative"]*np.r_[d["opacity"],d["opacity"][-1]]/(4*PI*G*state.enclosed*2.99792458e10))
