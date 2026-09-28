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
from .winds import rate as wind_rate
from .ledger import residuals
from .structure import hydrostatic_residual

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
        for p in sorted(Path(__file__).parent.rglob("*")):
            if p.is_file() and "__pycache__" not in p.parts: z.write(p,p.relative_to(Path(__file__).parent.parent))
    meta=dict(version=__version__,git_commit=git("rev-parse","HEAD"),git_status=git("status","--porcelain"),
              source_sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest(),python=sys.version,
              platform=platform.platform(),created_utc=datetime.now(timezone.utc).isoformat())
    (path/"config.yaml").write_text(yaml.safe_dump(config,sort_keys=False))
    (path/"working_tree.diff").write_text(git("diff","HEAD"))
    import importlib.metadata
    meta["dependencies"]={p:importlib.metadata.version(p) for p in ("numpy","scipy","PyYAML","matplotlib")}
    json_write(path/"metadata.json",meta)
    registry=Path(__file__).parent/"data"/"sources.json"
    if registry.exists(): (path/"source_prescriptions.json").write_bytes(registry.read_bytes())

def save_checkpoint(path,state,config):
    path=Path(path)
    # NPZ contains JSON strings, never pickled Python objects.
    temporary=path.with_suffix(".tmp")
    with temporary.open("wb") as f:
        np.savez_compressed(f,schema_version=4,time=state.time,bh_mass=state.bh_mass,dm=state.dm,r=state.r,v=state.v,T=state.T,
            cell_volume=state.volume,
            face_luminosity=state.face_luminosity if state.face_luminosity is not None else np.array([]),
            reference_mass_fractions=state.reference_mass_fractions if state.reference_mass_fractions is not None else state.dm/state.dm.sum(),
            ledger=json.dumps(state.ledger),steps=state.steps,next_dt=state.next_dt,config=json.dumps(config))
        f.flush(); os.fsync(f.fileno())
    temporary.replace(path)

def load_checkpoint(path):
    from .config import resolve
    with np.load(path,allow_pickle=False) as f:
        if int(f["schema_version"]) not in (1,2,3,4): raise ValueError("Unsupported checkpoint schema")
        state=State(float(f["time"]),float(f["bh_mass"]),f["dm"].copy(),f["r"].copy(),f["v"].copy(),f["T"].copy(),
                    json.loads(str(f["ledger"])),int(f["steps"]),float(f["next_dt"]))
        config=resolve(json.loads(str(f["config"])))
        if int(f['schema_version']) >= 2: state.cell_volume=f['cell_volume'].copy()
        if int(f['schema_version']) >= 3 and len(f['face_luminosity']):state.face_luminosity=f['face_luminosity'].copy()
        state.reference_mass_fractions=f['reference_mass_fractions'].copy() if int(f['schema_version'])>=4 else state.dm/state.dm.sum()
    if not (len(state.r)==len(state.dm)+1==len(state.v) and len(state.T)==len(state.dm)):
        raise ValueError("Malformed checkpoint dimensions")
    if np.any(np.diff(state.r)<=0) or np.any(state.dm<=0) or np.any(state.T<=0):
        raise ValueError("Invalid checkpoint physical state")
    geometric_volume=4*PI/3*np.diff(state.r)*(state.r[1:]**2+state.r[1:]*state.r[:-1]+state.r[:-1]**2)
    if state.volume.shape != state.dm.shape or np.any(state.volume<=0) or not np.allclose(state.volume,geometric_volume,rtol=1e-8,atol=0):
        raise ValueError('Checkpoint volume and geometry disagree')
    if state.face_luminosity is not None and (state.face_luminosity.shape!=state.r.shape or not np.all(np.isfinite(state.face_luminosity))):
        raise ValueError('Malformed checkpoint face luminosities')
    f=state.reference_mass_fractions
    if f.shape!=state.dm.shape or np.any(f<=0) or not np.isclose(f.sum(),1.,rtol=1e-12,atol=0):
        raise ValueError('Malformed checkpoint reference mass fractions')
    return state,config

def history_row(state,eos,opacity,config):
    e=rates(state,eos,config["engine"])
    L,d=luminosities(state,eos,opacity,config,e["luminosity"])
    radius=float(state.r[-1]); total=float(state.enclosed[-1])
    supply=value(config['feeding']['rate_msun_per_year'],state.time/YEAR)
    wind=wind_rate(state,eos,opacity,config['wind'],max(0.,float(L[-1])))
    er=state.energies(eos)
    transport=supply+wind*YEAR/MSUN+e['rest_mass_rate']*YEAR/MSUN
    row=dict(step=state.steps,time_years=state.time/YEAR,bh_mass_msun=state.bh_mass/MSUN,
        envelope_mass_msun=float(state.dm.sum()/MSUN),radius_cm=radius,
        luminosity_erg_s=float(L[-1]),engine_luminosity_erg_s=e["luminosity"],
        effective_temperature_diagnostic_K=float((max(L[-1],0)/(4*PI*radius**2*SIGMA))**.25),
        rest_accretion_msun_per_year=e["rest_mass_rate"]*YEAR/MSUN,
        bh_growth_msun_per_year=e["bh_growth_rate"]*YEAR/MSUN,
        supply_msun_per_year=supply,wind_msun_per_year=wind*YEAR/MSUN,
        photospheric_radius_cm=float('nan'),
        kelvin_helmholtz_time_years=abs(er['gravity'])/max(abs(float(L[-1])),1.)/YEAR,
        replacement_time_years=float(state.dm.sum()/MSUN)/transport if transport else float('inf'),
        dynamical_time_years=(radius**3/(G*total))**.5/YEAR,
        eddington_ratio=float(L[-1]/eddington(total,d["opacity"][-1])),
        max_hydrostatic_residual=float(max(abs(hydrostatic_residual(state,eos,config["boundary"]["pressure_cgs"])))),
        max_convective_saturation=float(max(d["saturation"],default=0.)),
        envelope_entropy_erg_K=float(np.dot(state.dm,eos.evaluate(state.rho,state.T)['s'])),
        **residuals(state,eos))
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
