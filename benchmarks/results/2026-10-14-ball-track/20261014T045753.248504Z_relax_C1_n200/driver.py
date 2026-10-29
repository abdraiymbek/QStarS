#!/usr/bin/env python3
"""Gate 6: registered replication of the Ball et al. 2011 fiducial quasi-star track.

Plan: benchmarks/ball2011_track_plan.json, registered by sha256 in
benchmarks/ball2011_track_plan.registration.json. Subcommands (the spec
Stage 7 order):

  relax     relax the registered seeds at one resolution (plan initial_model):
            convergence, seed uniqueness of row 0, the row-0 checkpoint (with
            its sha256 and the screens of row 0 in relax.json); relax.json is
            always written (an unexpected error fails that seed only)
  smoke     plan stop_rules.smoke_gate: an N = 200 track to M_BH = 100;
            smoke.json is always written
  track     one registered run (C0-C7, S1, S2, S4, V1-V3) from a row-0
            checkpoint to its termination, with the frozen-entropy fold
            classification (fold.py) when the termination trigger fires; with
            --restart, a restart segment (plan restart) from a checkpoint
  campaign  the registered runs in sequence, the restart segments R1-R3 of C0,
            then compare
  compare   comparison.json (C0 against Table 1, qualitative checks, endpoint,
            screen intervals) and campaign.json (convergence, conservation,
            restarts, verdict per plan verdict_rules, rejection causes per run,
            and the code of every run directory: comparisons between runs hold
            only between runs made with one code, source_consistency)

smoke, track and campaign start only from the designated row-0 checkpoint of a
relax run that ended converged_unique (row0_record: plan
initial_model.uniqueness), relaxed with the same configuration under the
registered plan and the same registered amendments, and, for registered
output, itself a registered run directory. A start row that violates an
always-enforced screen ends the run as a model_domain_exit of row 0 (a
nonfinite start state: numerical_failure). track.json is written whatever
happens to the run or to its fold classification. Every run directory
records the code it was made with (code_closure.json, finish()).

Every subcommand first recomputes the sha256 of the plan and of both target
CSVs and refuses to run unless they equal the registration record; every
amendment applied (--amendment) is verified against its own registration
record, must not be a draft, must cite a preserved failing run directory
whose record shows the failure, and may change only numerical-method
settings (AMENDABLE; the relaxation mesh cycles, the grey seed surface and
the photosphere bracket scan are off unless amended). --fallback F1 switches on the plan's first fallback
(stop_rules.fallbacks_in_order: pseudo-transient continuation inside a step,
quasistatic.ptc) for smoke, track and campaign; it is recorded in run.json
and in the configuration, and the relaxation never uses it. Every
run configuration is generated from the plan and resolved by qstars.config
(numerics.mechanics quasi_static); every subcommand also refuses unless
benchmarks/ball2011_track.yaml, the C0 configuration in the qstars schema,
equals that mapping (its MESA paths aside). Every run
directory, including failed and exited runs, is created under
benchmarks/results/<UTC date>-ball-track/ and keeps the plan registration,
the amendments, the driver source, the configuration and the source snapshot.
MESA EOS and the ABI-7 STARS opacity switch are required for every run.
"""
import os
for _name in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[_name]='1'
import argparse,copy,csv,hashlib,json,math,sys
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PLAN=ROOT/'benchmarks/ball2011_track_plan.json'
REGISTRATION=ROOT/'benchmarks/ball2011_track_plan.registration.json'
# C0 in the qstars configuration schema (Stage 6); the driver generates every run, C0 included, from
# the plan and refuses unless this file equals its C0 mapping (check_track_config).
TRACK_CONFIG=ROOT/'benchmarks/ball2011_track.yaml'
YEAR=365.25*86400
# Numerical-method settings an amendment may change (plan registration.amendment_rule:
# tolerances, primary physics, microphysics and targets are never amended).
AMENDABLE={'quasistatic':{'base','base_mass_fraction','photosphere_bracket'},'stepcontrol':{'metric'},
           'relaxation':{'mesh_cycles','mesh_tolerance','cycle_dt_years'},
           # Seed construction (plan initial_model.seed, not an immutable section): a seed is only the
           # initial guess of a relaxation whose result is tested for seed independence (uniqueness).
           'replication':{'seed_surface','seed_radii_rsun'}}
# Plan stop_rules.fallbacks_in_order that need no amendment: F1 (pseudo-transient continuation inside a
# step). F2 is introduced by amendment (quasistatic.base); F3 is the registered secant fallback of land().
FALLBACKS={'F1':('quasistatic','ptc',True)}
FALLBACK_TEXTS={'F1':'F1 pseudo-transient continuation inside a step'}
# The registered relaxation (plan initial_model.relaxation and converged_when, checked against
# RELAXATION_TEXTS). max_steps and max_retries are the driver's limits (the plan sets none). The STARS-Q
# mesh cycles of quasistatic.relax are not in the plan: off unless a registered amendment sets them.
RELAXATION=dict(dt_years=.01,dt_max_years=1e8,converged_dt_years=1e6,tolerance=1e-10,luminosity_tolerance=1e-6,
                max_steps=400,max_retries=20,mesh_cycles=0,mesh_tolerance=1e-5,cycle_dt_years=1e6)
RELAXATION_TEXTS={'relaxation':'M_BH frozen; dt starts at 0.01 yr and doubles up to 1e8 yr; injected L_BH dt is booked as relaxation_injected_energy',
                  'converged_when':'dt >= 1e6 yr, and two consecutive accepted steps have max |dln T| and |dln r| <= 1e-10 and |L*/L_BH - 1| <= 1e-6'}
# Plan initial_model.failure: a relaxation blocker is numerical unless a converging path reaches
# Gamma_ph >= 1. Evidence the driver requires for that label (its thresholds are driver choices,
# recorded in relax.json): every rejection of the final consecutive streak left the domain through the
# photosphere (DomainError message), and the photospheric Eddington factor of the last accepted steps
# never decreased and ended at or above the threshold.
GAMMA_PH_EXIT=dict(domain_message='Gamma_ph >= 1',last_accepted_steps=5,eddington_factor=.99)
# Registered must_cross items of the restart segments (plan restart.must_cross).
MUST_CROSS=('base_slice','split_merge','rejected_step')
CLAMP_LIMIT_MSUN=.01  # plan replication_mode: opacity_clamp mass at every compared row
# Directories an amendment may cite as its preserved failing run, and the record that shows the failure.
FAILURE_RECORDS=('smoke.json','relax.json','termination.json')
# Registered runs: the exact plan text (checked, so a different plan cannot silently
# reuse this mapping) and the configuration change it states.
RUNS={
    'C0':('convergence_campaign','production: N = 400, delta_target = 0.005, tolerance 1e-12, face boundary state, split/merge 2.0/0.5',{}),
    'C1':('convergence_campaign','N = 200',dict(zones=200)),
    'C2':('convergence_campaign','N = 800',dict(zones=800)),
    'C3':('convergence_campaign','delta_target = 0.01',dict(delta_target=.01)),
    'C4':('convergence_campaign','delta_target = 0.0025',dict(delta_target=.0025)),
    'C5':('convergence_campaign','tolerance 1e-10',dict(tolerance=1e-10)),
    'C6':('convergence_campaign','boundary_state = cell',dict(boundary_state='cell')),
    'C7':('convergence_campaign','split/merge 1.5/0.67',dict(split=1.5,merge=.67)),
    'S1':('sensitivity','MESA default opacity handle (OPAL Type 2 + Ferguson 2005 low T)',dict(tables='default')),
    'S2':('sensitivity','mlt_thesis_lossless: thesis eqs 2.6-2.9, (l/Hp)^2, delta = 1, alpha = 2, no cap',dict(convection='mlt_thesis_lossless')),
    'S4':('sensitivity','strict opacity domain, run only if C0 clamped any cell',dict(tables='stars2004_switch_strict')),
    'V1':('secondary','rho ~ r^(-1/2) cavity, M_cav = (8 pi/5) rho_b r0^3 (thesis 3.33, paper eq 46)',dict(profile='r-1/2')),
    'V2':('secondary','M_cav = 0 (thesis 3.3.2)',dict(profile='none')),
    'V3':('secondary','M* = 1e5 Msun',dict(total_mass_msun=1e5))}
# Plan comparison observables: Table 1 column (printed unit) and history column (QStarS unit).
OBSERVABLES={'t':('time_myr',1.,'time_myr'),'Mdot_BH':('bh_growth_1e-4_msun_year',1e-4,'mdot_bh_msun_per_year'),
             'M_cav':('cavity_mass_msun',1.,'cavity_mass_msun'),'L*':('luminosity_1e8_lsun',1e8,'luminosity_lsun'),
             'rho0':('rho_inner_cgs',1.,'base_density_cgs'),'T0':('T_inner_1e5_K',1e5,'base_temperature_K'),
             'Teff':('Teff_1e3_K',1e3,'effective_temperature_K'),'r_B':('r_inner_100_rsun',100.,'inner_radius_rsun'),
             'R*':('r_outer_1e4_rsun',1e4,'radius_rsun')}
PROFILE_FRACTIONS=(np.arange(20)+.5)/20  # 20 fixed exterior-mass fractions, equal mass weights
SEGMENT_STEPS=100  # plan restart.segment
RESTART_POINTS=(('R1','restart_bh20msun.npz'),('R2','restart_before_landing_500msun.npz'),('R3','restart_bh900msun.npz'))


# Plan texts whose numbers the comparison hard-codes (checked, like RUNS).
MAPPED_TEXTS={('endpoint','criteria','M_BH_end'):'|ln(M_end / 1194)| <= 0.10',('endpoint','criteria','t_end'):'|ln(t_end / 4.23 Myr)| <= 0.10',
              ('qualitative','Q1'):'every accepted step after t = 0 lies in [1.6e-4, 4.0e-4] Msun/yr',
              ('qualitative','Q2'):'outward density increase > 1e-3 relative present at rows 1-4 and absent at row 0',
              ('qualitative','Q3'):'integral of (3 Gamma1 - 4) p dV > 0 at rows 0-4 and at the last accepted model',
              ('qualitative','Q4'):'mass above the outermost convective boundary < 100 Msun at rows 2-4',
              ('qualitative','Q5'):'mass-weighted Gamma1 over convective cells in [1.33, 1.36] at rows 1-4',
              ('endpoint','classification_protocol',2):'the fold must lie within 2% of the terminal M_BH and agree between N = 400 and N = 800 to 2%',
              ('replication_mode','report_only_screens',4):'opacity_clamp (clamped mass must be <= 0.01 Msun at every compared row)',
              ('restart','must_cross'):['a base slice','a split/merge event','a rejected step','a landing (restart point 2)'],
              ('convergence_campaign','thresholds','observed_order'):'reported; expected >= 0.8 and monotone'}


class Refused(SystemExit):
    """The driver refuses to run (registration or amendment mismatch)."""


def sha256(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def verify_registration(plan=PLAN,registration=REGISTRATION,root=ROOT):
    """Refuse unless the plan and both target CSVs match their registered sha256."""
    plan,registration,root=Path(plan),Path(registration),Path(root)
    record=json.loads(registration.read_text());digest=sha256(plan)
    if digest!=record['plan_sha256']:
        raise Refused(f'Refusing to run: {plan} sha256 {digest} differs from the registered {record["plan_sha256"]}')
    for name,expected in record['target_files_sha256'].items():
        if sha256(root/name)!=expected:raise Refused(f'Refusing to run: {name} differs from its registered sha256')
    return dict(plan=str(plan),plan_sha256=digest,registration=str(registration),registration_sha256=sha256(registration),
                registered_utc=record['registered_utc'],target_files_sha256=record['target_files_sha256'])

def verify_amendment(path,plan_sha256,root=ROOT):
    """An amendment JSON and its <stem>.registration.json (amendment_sha256, plan_sha256).

    Its 'settings' mapping may name only AMENDABLE numerical-method keys, its
    status must not be a draft, and preserved_failing_run (a path, or a list)
    must name run directories under <root>/benchmarks/results/<date>-ball-track/
    whose record shows the failure (plan registration.amendment_rule).
    Returns (record, settings).
    """
    path=Path(path);registration=path.with_name(path.stem+'.registration.json')
    if not registration.exists():raise Refused(f'Refusing to run: amendment {path} has no registration record {registration.name}')
    record=json.loads(registration.read_text());digest=sha256(path)
    if digest!=record.get('amendment_sha256'):raise Refused(f'Refusing to run: amendment {path} differs from its registered sha256')
    if record.get('plan_sha256')!=plan_sha256:raise Refused(f'Refusing to run: amendment {path} was registered against another plan')
    amendment=json.loads(path.read_text());settings=amendment.get('settings')
    status=amendment.get('status')
    if not isinstance(status,str) or status.strip().upper().startswith('DRAFT'):
        raise Refused(f'Refusing to run: amendment {path} is a draft or has no status ({status!r})')
    runs=amendment.get('preserved_failing_run');runs=[runs] if isinstance(runs,str) else runs
    if not runs or not isinstance(runs,list):raise Refused(f'Refusing to run: amendment {path} cites no preserved failing run')
    for run in runs:
        failure=preserved_failure(run,root)
        if failure is not True:raise Refused(f'Refusing to run: amendment {path} preserved_failing_run {run!r}: {failure}')
    if not isinstance(settings,dict) or not settings:raise Refused(f'Refusing to run: amendment {path} has no machine-readable settings')
    for section,values in settings.items():
        if section not in AMENDABLE or not isinstance(values,dict) or set(values)-AMENDABLE[section]:
            raise Refused(f'Refusing to run: amendment {path} changes {section} {sorted(values) if isinstance(values,dict) else values}, '
                          'which is not an amendable numerical-method setting')
    return dict(path=str(path),sha256=digest,registration=str(registration),registration_sha256=sha256(registration),
                registered_utc=record.get('registered_utc'),preserved_failing_run=runs,settings=settings),settings

def preserved_failure(run,root=ROOT):
    """True when `run` is a run directory under <root>/benchmarks/results/<date>-ball-track/ whose
    smoke.json, relax.json or termination.json (first found) records a failure; else the reason.

    A failure is a smoke gate not passed, a relax run that ended relaxation_blocker or not_unique, or a
    track that did not complete. An interrupted or requested stop (status control_stop) is no failure:
    such a directory shows no defect and cannot back an amendment.
    """
    if not isinstance(run,str) or not run:return 'not a path'
    results=(Path(root)/'benchmarks/results').resolve();d=(Path(root)/run).resolve()
    if not (d.is_dir() and d.parent.parent==results and d.parent.name.endswith('-ball-track')):
        return 'not a run directory under benchmarks/results/<date>-ball-track/'
    name=next((f for f in FAILURE_RECORDS if (d/f).is_file()),None)
    if name is None:return f'no {", ".join(FAILURE_RECORDS)} in the directory'
    data=json.loads((d/name).read_text())
    status=data.get('status')
    failed={'smoke.json':data.get('passed') is not True and status!='control_stop',
            'relax.json':status in ('relaxation_blocker','not_unique'),
            'termination.json':status not in ('completed','control_stop')}[name]
    return True if failed else f'{name} records no failure'

def check_mapped_texts(plan):
    """Refuse when a plan text the comparison maps to numbers differs from MAPPED_TEXTS."""
    checks={q['id']:q['test'] for q in plan['qualitative_checks']}
    for key,text in MAPPED_TEXTS.items():
        if key[0]=='qualitative':found=checks.get(key[1])
        else:
            found=plan
            for k in key:found=found[k] if isinstance(found,(dict,list)) else None
        if found!=text:raise Refused(f'Refusing to compare: plan text {key} is {found!r}, the driver maps {text!r}')

def run_text(plan,run_id):
    if run_id not in RUNS:raise ValueError(f'Unknown registered run {run_id}; choose from {sorted(RUNS)}')
    kind,text,_=RUNS[run_id]
    if kind=='convergence_campaign':found=plan['convergence_campaign']['runs'].get(run_id)
    elif kind=='sensitivity':found=next((x['change'] for x in plan['sensitivity_runs_reported_only'] if x.get('id')==run_id),None)
    else:found=next((x['change'] for x in plan['secondary_published_checks_after_primary_verdict'] if x.get('id')==run_id),None)
    if found!=text:raise Refused(f'Refusing to run: plan text of {run_id} is {found!r}, the driver maps {text!r}')
    return kind,text

def build_config(plan,run_id,library='',mesa_dir='',amendments=()):
    """The registered configuration of one run (plan model, numerics and the run's change), resolved and
    cross-validated by qstars.config (numerics.mechanics quasi_static); C0 equals benchmarks/ball2011_track.yaml."""
    from qstars.config import resolve
    kind,text=run_text(plan,run_id);change=RUNS[run_id][2]
    model,numerics=plan['model'],plan['numerics'];controller=numerics['controller']
    zones=change.get('zones',400)
    for key,expected in RELAXATION_TEXTS.items():
        if plan['initial_model'][key]!=expected:raise Refused(f'Refusing to run: plan initial_model.{key} differs from the driver mapping')
    law=model['accretion_law'];stop=plan['stop_rules']
    raw=dict(name=f'ball2011_{run_id}_n{zones}',composition=model['composition'],eos='mesa',
             mesa=dict(library_path=str(library),mesa_dir=str(mesa_dir)),opacity=dict(model='mesa',tables=change.get('tables','stars2004_switch')),
             engine=dict(model=law['engine_model'],coupling=law['eta'],efficiency=law['epsilon']),
             transport=dict(convection=change.get('convection',model['convection']['closure']),mixing_length_alpha=model['convection']['alpha']),
             boundary=dict(inner='bondi_cavity',thermal='stars_photosphere',pressure_cgs=0.),
             numerics=dict(zones=zones,mechanics='quasi_static',conservation_tolerance=plan['conservation']['energy_cumulative_over_initial_internal_energy']),
             validity=dict(min_temperature_K=3000.,max_temperature_K=1e8,max_mach=.2,max_compactness=1e-3,require_optically_thick=True),
             evolution=dict(max_steps=stop['max_accepted_steps'],max_wall_seconds=stop['per_track_wall_hours']*3600.),
             cavity=dict(active=True,profile=change.get('profile','r-3/2'),b=model['inner_radius']['b'],boundary_state=change.get('boundary_state','face')),
             quasistatic=dict(tolerance=change.get('tolerance',numerics['solver']['tolerance']),
                              polish_iterations=numerics['solver']['polish_iterations'],roundoff_residual=1e-10,roundoff_step=1e-13,
                              max_passes=3,split=change.get('split',2.),merge=change.get('merge',.5),base='lagrangian',
                              energy_per_step=plan['conservation']['energy_cumulative_over_initial_internal_energy']/stop['max_accepted_steps'],
                              mass_per_step=plan['conservation']['mass_per_step_relative']),
             stepcontrol=dict(delta_target=change.get('delta_target',controller['delta_target']),
                              dt_initial_years=controller['dt_initial_years'],dt_max_years=controller['dt_max_years'],max_retries=12,metric='lagrangian'),
             relaxation=dict(RELAXATION),
             replication=dict(plan=plan['plan_id'],run=run_id,kind=kind,plan_text=text,policy=plan['replication_mode']['mode'],
                              total_mass_msun=change.get('total_mass_msun',model['total_mass_msun']),
                              initial_bh_mass_msun=model['initial_bh_mass_msun'],
                              seed_radii_rsun=plan['initial_model']['seed_radii_rsun']+plan['initial_model']['fallback_seed_radii_rsun'],
                              landings_msun=numerics['landing']['bh_masses_msun'],end_bh_mass_msun=0.,amendments=[]))
    for record,settings in amendments:
        for section,values in settings.items():raw[section].update(values)
        raw['replication']['amendments'].append(record['path'])
    return resolve(raw)

def _flat(d,prefix=''):
    out={}
    for k,v in d.items():
        if isinstance(v,dict):out.update(_flat(v,prefix+k+'.'))
        else:out[prefix+k]=v
    return out

def check_track_config(plan,path=None):
    """Refuse unless benchmarks/ball2011_track.yaml (its MESA paths aside) is the C0 configuration build_config
    maps from the plan; returns its path and sha256 (recorded in every run directory)."""
    import yaml
    from qstars.config import load
    path=Path(path or TRACK_CONFIG)
    try:c=load(path)
    except (OSError,ValueError,TypeError,KeyError,yaml.YAMLError) as exc:
        raise Refused(f'Refusing to run: {path.name} is not a valid qstars configuration ({exc})') from exc
    c['mesa']=dict(library_path='library',mesa_dir='mesa_dir')  # the driver substitutes --library and --mesa-dir
    a,b=_flat(c),_flat(build_config(plan,'C0','library','mesa_dir'))
    differ=sorted(k for k in set(a)|set(b) if a.get(k,'<missing>')!=b.get(k,'<missing>'))
    if differ:raise Refused(f'Refusing to run: {path.name} differs from the registered C0 configuration at {differ}')
    return dict(path=str(path),sha256=sha256(path))

def physics(config):
    """MESA 24.08.1 EOS and the run's opacity (opacity.tables: ABI 7 STARS switch, clamped per the STARS rule unless strict)."""
    from qstars.eos import Composition
    from qstars.mesa_microphysics import MesaMicrophysics
    from qstars.quasistatic import opacity_tables
    eos=MesaMicrophysics(Composition(**config['composition']),config['mesa']['library_path'],config['mesa']['mesa_dir'])
    if eos.abi_version<7:raise Refused('Refusing to run: the registered STARS opacity switch needs a MESA bridge ABI 7 build')
    return eos,opacity_tables(config,eos,eos.opacity)

def run_directory(parent,kind,run_id,config,registration,amendments,extra=None):
    from qstars import io
    path=io.create_directory(parent,f"{kind}_{run_id}_n{config['numerics']['zones']}")
    io.provenance(path,config)
    io.json_write(path/'plan_registration.json',registration)
    io.json_write(path/'amendments.json',[a for a,_ in amendments])
    (path/'driver.py').write_bytes(Path(__file__).read_bytes())
    library=Path(config['mesa']['library_path'])
    io.json_write(path/'run.json',dict(kind=kind,run=run_id,replication=config['replication'],
        track_config=dict(path=str(TRACK_CONFIG.relative_to(ROOT)),sha256=sha256(TRACK_CONFIG)),
        threads={k:os.environ[k] for k in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS')},
        library=str(library),library_sha256=sha256(library) if library.is_file() else None,**(extra or {})))
    print(path,flush=True)
    return path

def signature(config):
    """The configuration a row-0 model depends on: everything but the step controller, labels, the MESA
    paths and fallback F1 (quasistatic.ptc), which the relaxation never uses (bondi_step mode 'relax')."""
    c=copy.deepcopy(config)
    for key in ('stepcontrol','name'):c.pop(key,None)
    for key in ('run','kind','plan_text'):c.get('replication',{}).pop(key,None)
    c.get('mesa',{}).pop('library_path',None);c.get('mesa',{}).pop('mesa_dir',None)
    c.get('quasistatic',{}).pop('ptc',None)
    if c.get('quasistatic')=={}:c.pop('quasistatic')  # a configuration that held nothing but the switch
    return json.dumps(c,sort_keys=True)

def apply_fallbacks(config,fallbacks=()):
    """Switch on the registered fallbacks named on the command line (FALLBACKS); returns the configuration."""
    for name in fallbacks:
        if name not in FALLBACKS:raise Refused(f'Refusing to run: unknown fallback {name}; choose from {sorted(FALLBACKS)}')
        section,key,value=FALLBACKS[name]
        config.setdefault(section,{})[key]=value
    return config

def fallbacks_of(config):
    """The fallbacks a configuration has switched on (for run.json and restart checks)."""
    return [name for name,(section,key,value) in FALLBACKS.items() if (config.get(section) or {}).get(key)==value]

def row0_record(path,provenance=None):
    """(relax.json, None) when the relax run next to a row-0 checkpoint vouches for it, else (None, reason).

    Plan initial_model.uniqueness: a row-0 model is used only when the relax
    run that made it ended converged_unique, names this file as its row-0
    checkpoint and recorded its sha256, and its directory is a relax run
    directory (run.json of kind relax). provenance (smoke, track and
    campaign always pass it; row0_provenance): that run was made under the
    registered plan (plan_registration.json: the plan and registration
    sha256 of this run) with the same registered amendments (amendments.json:
    the same amendment sha256 in the same order), and, for registered output
    (results_root set), it is a run directory under
    <results_root>/<date>-ball-track/. A hand-edited record can still pass
    these checks; every run records the relax directory and the sha256 of
    its relax.json (run.json row0_relax) so the product can be audited.
    """
    path=Path(path);d=path.parent;record=d/'relax.json'
    if not record.is_file():return None,'no relax.json next to it (not the product of a relax run)'
    data=json.loads(record.read_text())
    if data.get('status')!='converged_unique':
        return None,f"its relax run ended {data.get('status')!r}, not 'converged_unique' (plan initial_model.uniqueness)"
    if data.get('row0_checkpoint')!=path.name:return None,f"its relax run names {data.get('row0_checkpoint')!r} as the row-0 checkpoint"
    if data.get('row0_sha256')!=sha256(path):return None,'its sha256 differs from the one its relax run recorded'
    run=_read_json(d/'run.json')
    if not isinstance(run,dict) or run.get('kind')!='relax':return None,'no run.json of kind relax next to it (not a relax run directory)'
    if provenance is not None:
        stored=_read_json(d/'plan_registration.json') or {};expected=provenance['registration']
        if not isinstance(stored,dict) or stored.get('plan_sha256')!=expected['plan_sha256'] or \
                stored.get('registration_sha256')!=expected['registration_sha256']:
            return None,'its relax run was not made under the registered plan (plan_registration.json differs)'
        stored=_read_json(d/'amendments.json')
        applied=[r['sha256'] for r,_ in provenance['amendments']]
        if not isinstance(stored,list) or [a.get('sha256') if isinstance(a,dict) else None for a in stored]!=applied:
            return None,'its relax run applied other registered amendments (amendments.json sha256 differ from this run\'s)'
        if provenance.get('results_root') is not None:
            results=Path(provenance['results_root']).resolve();rd=d.resolve()
            if not (rd.parent.parent==results and rd.parent.name.endswith('-ball-track')):
                return None,f'its relax run is not a registered run directory under {results}/<date>-ball-track/'
    return data,None

def _read_json(path):
    path=Path(path)
    if not path.is_file():return None
    try:return json.loads(path.read_text())
    except ValueError:return None

def row0_provenance(parent,registration,amendments):
    """What row0_record checks for a run written under `parent`: the plan registration and amendments of
    this run, and results_root benchmarks/results when the output is the registered one (parent under it)."""
    results=(ROOT/'benchmarks/results').resolve()
    registered=Path(parent).resolve().is_relative_to(results)
    return dict(registration=registration,amendments=list(amendments),results_root=results if registered else None)

def row0_for(config,candidates,reasons=None,provenance=None):
    """The first row-0 checkpoint (t = 0) relaxed with this run's physics and numerics by a relax run
    that ended converged_unique (row0_record with provenance); None if none fits. reasons (a list)
    receives why each candidate was refused."""
    return row0_pick(config,candidates,reasons,provenance)[0]

def row0_pick(config,candidates,reasons=None,provenance=None):
    """(state, path, relax.json record) of the first candidate row0_for accepts, else (None, None, None)."""
    from qstars.quasistatic import load_checkpoint
    for path in candidates:
        if path is None:continue
        state,stored=load_checkpoint(path)
        if not (state.time==0 and state.steps==0 and len(state.dm)==config['numerics']['zones'] and signature(stored)==signature(config)):
            why='not a t = 0 model with this run\'s zone count and row-0 configuration'
        else:
            record,why=row0_record(path,provenance)
            if record is not None:return state,path,record
        if reasons is not None:reasons.append(f'{path}: {why}')
    return None,None,None

def row0_relax(path):
    """run.json entry of the relax run that vouched for a row-0 checkpoint (audit trail)."""
    d=Path(path).parent
    return dict(directory=str(d),relax_json_sha256=sha256(d/'relax.json'),checkpoint_sha256=sha256(path))

def track_options(plan,config):
    """Landings, restart points and stop rules of quasistatic.run (plan numerics.landing, restart, stop_rules)."""
    MSUN=1.98847e33;stop=plan['stop_rules']
    landings=[m*MSUN for m in plan['numerics']['landing']['bh_masses_msun']]
    if stop['dt_floor']!='dt < 10 yr for 500 consecutive accepted steps before M_BH = 1000 is a numerical blocker':
        raise Refused('Refusing to run: the plan dt floor differs from the driver mapping')
    if plan['restart']['segment']!='restart from the checkpoint and run 100 accepted steps, or to termination':
        raise Refused('Refusing to run: the plan restart segment differs from the driver mapping')
    return dict(landings=landings,policy=plan['replication_mode']['mode'],max_wall_seconds=stop['per_track_wall_hours']*3600.,
                max_accepted_steps=stop['max_accepted_steps'],dt_floor=dict(dt=10*YEAR,steps=500,before_bh_mass=1000*MSUN))

# ---------------------------------------------------------------- relax and smoke

def relax_seeds(plan,config,eos,opacity,path,seeds,zones,log=print):
    """Plan initial_model: relax each seed (config['relaxation']); the first two converged seeds must agree to 1e-6 on row 0.

    Without two converged seeds the blocker is numerical unless no seed
    converged and a relaxation path meets relaxation_domain_exit (Gamma_ph
    >= 1 reached on a converging path); seed construction failures never
    count as that evidence, and neither does a relaxation that ends in an
    error other than a step rejection: that seed has failed (no domain-exit
    evidence) and the next seed is tried. An unexpected error of one seed
    (any exception, in its seed, its relaxation or the post-processing of a
    converged seed: checkpoints, profile, observables, screens) is that
    seed's failure, recorded with its traceback; such a seed never counts as
    converged. relax.json is always written: an interrupt (KeyboardInterrupt,
    SystemExit) leaves status control_stop, which is no result and no
    failure (preserved_failure), and is raised again.
    """
    import traceback
    from qstars import ledger,quasistatic,io
    from qstars.constants import MSUN,RSUN
    from qstars.mass_transport import DomainError
    rep=config['replication'];results=[];converged=[]
    out=dict(zones=zones,relaxation=config['relaxation'],seeds=results,status='control_stop',
             classification='interrupted before the relaxation ended: no result, no failure')
    def failed(entry,stage,exc,status='failed'):
        entry.update(status=status,reason=f'{type(exc).__name__}: {exc}',error_stage=stage,traceback=traceback.format_exc())
        if status=='seed_failed':entry['seed_domain_error']=isinstance(exc,DomainError) and stage=='seed'
        else:entry['domain_exit_evidence']=dict(domain_exit=False,domain_trials=[])
    try:
        for radius in seeds:
            entry=dict(seed_radius_rsun=radius)
            try:
                state,info=quasistatic.seed(config,eos,opacity,zones,rep['initial_bh_mass_msun']*MSUN,rep['total_mass_msun']*MSUN,radius*RSUN)
                entry['seed']={k:v for k,v in info.items() if not isinstance(v,np.ndarray)}
            except (ValueError,RuntimeError,DomainError) as exc:
                entry.update(status='seed_failed',reason=f'{type(exc).__name__}: {exc}',seed_domain_error=isinstance(exc,DomainError))
                results.append(entry);log(f'seed {radius} Rsun failed: {exc}');continue
            except Exception as exc:  # unexpected: that seed fails, the record is kept
                failed(entry,'seed',exc,'seed_failed');entry['seed_domain_error']=False
                results.append(entry);log(f'seed {radius} Rsun failed: {type(exc).__name__}: {exc}');continue
            try:
                ledger.initialize(state,eos,cavity=True)
                with (path/f'relax_{radius:g}rsun.jsonl').open('w') as events:
                    state,record=quasistatic.relax(state,eos,opacity,config,config['relaxation'],
                                                   log=lambda kind,data:events.write(json.dumps(data)+'\n'))
            except (ValueError,RuntimeError,ArithmeticError,DomainError) as exc:  # not a step rejection: relax retries those
                entry.update(status='failed',reason=f'{type(exc).__name__}: {exc}',
                             domain_exit_evidence=dict(domain_exit=False,domain_trials=[]))
                results.append(entry);log(f'seed {radius} Rsun: relaxation error: {exc}');continue
            except Exception as exc:
                failed(entry,'relaxation',exc)
                results.append(entry);log(f'seed {radius} Rsun: relaxation error: {type(exc).__name__}: {exc}');continue
            entry.update(status=record['status'],reason=record['reason'],accepted_steps=len(record['steps']),
                         rejections=len(record['rejections']),mesh_cycles=record['cycles'],
                         final_rejection_streak=record.get('final_rejection_streak'),
                         final_rejections_all_domain_exits=record.get('final_rejections_all_domain_exits'))
            try:
                if record['status']!='converged':entry['domain_exit_evidence']=relaxation_domain_exit(record)
                quasistatic.save_checkpoint(path/f'relaxed_{radius:g}rsun.npz',state,config)
                if record['status']=='converged':
                    row0=quasistatic.zero_time(state,eos)
                    quasistatic.save_checkpoint(path/f'row0_{radius:g}rsun.npz',row0,config)
                    quasistatic.save_profile(path/f'profile_row0_{radius:g}rsun.npz',row0,eos,opacity,config)
                    done=dict(row0=row0_observables(row0,eos,opacity,config),row0_checkpoint=f'row0_{radius:g}rsun.npz')
                    done.update(row0_sha256=sha256(path/done['row0_checkpoint']),row0_zones=len(row0.dm),
                                row0_screens=row0_screens(row0,eos,opacity,config,plan['replication_mode']['mode']))
                    entry.update(done);converged.append(entry)
            except Exception as exc:  # post-processing of a relaxed seed: that seed fails, never converged
                relaxed=entry['status']
                for key in ('row0','row0_checkpoint','row0_sha256','row0_zones','row0_screens'):entry.pop(key,None)
                failed(entry,'post-processing of a '+str(relaxed)+' relaxation',exc)
                entry['relaxation_status']=relaxed
            results.append(entry);log(f"seed {radius} Rsun: {entry['status']}")
            if len(converged)==2:break
        out.pop('classification')
        if len(converged)>=2:
            a,b=converged[0]['row0'],converged[1]['row0']
            diffs={k:abs(math.log(a[k]/b[k])) for k in a}
            out['uniqueness']=dict(seeds=[converged[0]['seed_radius_rsun'],converged[1]['seed_radius_rsun']],tolerance=1e-6,
                                   max_abs_dln=diffs,unique=all(v<=1e-6 for v in diffs.values()))
            out['row0_checkpoint']=converged[0]['row0_checkpoint'];out['row0_sha256']=converged[0]['row0_sha256']
            # Reported here, enforced by the track: a row-0 model that violates an always-enforced screen ends
            # every run started from it as a model_domain_exit of row 0 (quasistatic.run).
            out['row0_enforced_screen_errors']=converged[0]['row0_screens']['enforced_screen_errors']
            out['status']='converged_unique' if out['uniqueness']['unique'] else 'not_unique'
        else:
            evidence=[e['domain_exit_evidence'] for e in results if 'domain_exit_evidence' in e]
            domain=not converged and any(x['domain_exit'] for x in evidence)
            out['classification']=('model_domain_exit of row 0 (a converging path reaches Gamma_ph >= 1)' if domain
                                   else 'numerical relaxation blocker')
            out['domain_exit_rule']=dict(GAMMA_PH_EXIT,converged_seeds=len(converged),
                                         note='a row-0 model exists when any seed converged; thresholds are driver choices')
            out['domain_trials']={str(e['seed_radius_rsun']):e['domain_exit_evidence']['domain_trials'] for e in results if 'domain_exit_evidence' in e}
            out['status']='relaxation_blocker'
    except Exception as exc:  # the assembly itself failed: a numerical blocker, never converged_unique
        out.update(status='relaxation_blocker',classification='numerical relaxation blocker',
                   driver_error=f'{type(exc).__name__}: {exc}',traceback=traceback.format_exc())
        for key in ('uniqueness','row0_checkpoint','row0_sha256','row0_enforced_screen_errors'):out.pop(key,None)
    finally:
        io.json_write(path/'relax.json',_clean(out))
    return out

def relaxation_domain_exit(record,rule=GAMMA_PH_EXIT):
    """Evidence that a failed relaxation reached Gamma_ph >= 1 on a converging path (see GAMMA_PH_EXIT).

    The final consecutive rejections (record['final_rejection_streak']) must
    all carry the photospheric domain message, and the photosphere Eddington
    factor of the last `last_accepted_steps` accepted steps must be
    nondecreasing and end at or above `eddington_factor`.
    """
    streak=int(record.get('final_rejection_streak') or 0)
    final=record['rejections'][-streak:] if streak else []
    photosphere=bool(final) and all(rule['domain_message'] in (r.get('domain_message') or '') for r in final)
    k=rule['last_accepted_steps'];factors=[s.get('photosphere_eddington_factor') for s in record['steps'][-k:]]
    rising=len(factors)==k and all(f is not None for f in factors) and all(b>=a for a,b in zip(factors,factors[1:])) \
        and factors[-1]>=rule['eddington_factor']
    return dict(final_rejections=len(final),final_rejections_photosphere=photosphere,
                domain_trials=[int(r.get('domain_trials') or 0) for r in final],last_eddington_factors=factors,
                eddington_factor_rising_to_boundary=bool(rising),domain_exit=bool(photosphere and rising))

def row0_screens(state,eos,opacity,config,policy):
    """Screens of a relaxed row-0 model under the plan policy (reported in relax.json; the track enforces them)."""
    from qstars import validity
    try:
        errors,diag,_=validity.check(state,eos,opacity,config,policy=policy,face_luminosity=state.face_luminosity)
    except (RuntimeError,ValueError,ArithmeticError,KeyError) as exc:  # StepRejected (geometry, conservation) is a RuntimeError
        return dict(policy=str(policy),enforced_screen_errors=[f'screens could not be evaluated: {type(exc).__name__}: {exc}'],report_only=[])
    return dict(policy=diag.get('policy'),enforced_screen_errors=list(errors),report_only=list(diag.get('screen_warnings',[])))

def row0_observables(state,eos,opacity,config):
    """Row-0 observables of the plan comparison (relaxed model, t = 0), in history units."""
    from qstars.quasistatic import history_row
    row=history_row(state,eos,opacity,config,event='start')
    return {name:float(row[col]) for name,(_,_,col) in OBSERVABLES.items() if name!='t'}

def smoke_summary(path,report,plan):
    """Plan stop_rules.smoke_gate: M_BH = 100 reached, median dt >= 50 yr, <= 20% rejected attempts."""
    rows=read_history(path) if (Path(path)/'history.csv').is_file() else []  # a run that died before its start row
    dts=[float(r['step_dt']) for r in rows if r['event'] in ('step','landing')]
    attempts=report.get('attempts',0);rejected=report.get('rejections',0)/attempts if attempts else None
    reached=report['status']=='completed'
    median=float(np.median(dts)) if dts else None
    out=dict(rule=plan['stop_rules']['smoke_gate'],reached_bh_mass_100=reached,status=report['status'],reason=report['reason'],
             accepted_steps=len(dts),attempts=attempts,rejections=report.get('rejections',0),median_dt_years=median,
             history_written=bool(rows),
             rejected_fraction=rejected,passed=bool(reached and median is not None and median>=50 and rejected is not None and rejected<=.2))
    return out

# ---------------------------------------------------------------- history and profiles

def read_history(path):
    with open(Path(path)/'history.csv',newline='') as f:return list(csv.DictReader(f))

def number(x):
    return None if x in ('',None) else float(x)

def landed_rows(rows):
    """Rows 1-4: the landing rows, keyed by target (Msun)."""
    return {float(r['landed_target_msun']):r for r in rows if r['event']=='landing'}

def profile_quantities(path):
    """Mass-coordinate profiles of a saved profile for the campaign norms and checks Q4, Q5."""
    with np.load(path) as f:
        p={k:f[k].copy() for k in f.files}
    dm=p['dm_g'];menv=float(np.sum(dm))
    exterior=(menv-(np.cumsum(dm)-.5*dm))/menv
    L=p['luminosity_edges_erg_s'];conv=np.r_[False,p['convective_interior_edges'],False]
    cells=conv[:-1]|conv[1:]
    rc=(.5*(p['r_edges_cm'][:-1]**3+p['r_edges_cm'][1:]**3))**(1/3)
    outer=np.flatnonzero(conv)
    return dict(exterior=exterior,ln_rho=np.log(p['density_cgs']),ln_T=np.log(p['temperature_K']),ln_r=np.log(rc),
                L_over_LBH=.5*(L[:-1]+L[1:])/L[0],dm=dm,gamma1=p['gamma1'],convective_cells=cells,
                mass_above_outer_convective_boundary_msun=float(np.sum(dm[outer[-1]:]))/1.98847e33 if len(outer) else menv/1.98847e33,
                convective_gamma1=float(np.sum(dm[cells]*p['gamma1'][cells])/np.sum(dm[cells])) if np.any(cells) else None)

def profile_norm(a,b,key):
    """Mass-weighted L1 difference at PROFILE_FRACTIONS of exterior mass (equal mass weights)."""
    fa=np.interp(PROFILE_FRACTIONS,a['exterior'][::-1],a[key][::-1]);fb=np.interp(PROFILE_FRACTIONS,b['exterior'][::-1],b[key][::-1])
    return float(np.mean(abs(fa-fb)))

# ---------------------------------------------------------------- comparison against Table 1

def rounding(text):
    """Half a unit in the last printed digit of a CSV string, relative to its value (inf for 0)."""
    text=text.strip();mantissa,_,exponent=text.lower().partition('e')
    decimals=len(mantissa.split('.')[1]) if '.' in mantissa else 0
    value=float(text)
    return math.inf if value==0 else .5*10**(-decimals)*10**(int(exponent) if exponent else 0)/abs(value)

def table_rows(path):
    with open(path,newline='') as f:return list(csv.DictReader(f))

def spread(plan,name):return plan['comparison']['columns'][name]['realization_spread']

def cells(plan,targets,rows):
    """The 44 plan cells: d = |ln(Q_qstars/Q_ball)| against tau_model + realization spread + rounding.

    rows: {0: row-0 history row, 100: landing row, ...} (history.csv dicts).
    """
    out=[];columns=plan['comparison']['columns'];strict=plan['comparison']['allowance_components']['numerical']
    report_only=set(_report_only(plan));enforceable=_enforceable(plan)
    for index,target in enumerate(targets[:5]):
        mass=float(target['bh_mass_msun']);row=rows.get(0 if index==0 else mass)
        names=plan['comparison']['cells']['row_0' if index==0 else 'rows_1_to_4']
        for name in names:
            cell=dict(row=index,bh_mass_msun=mass,observable=name)
            if name.startswith('M_cav <'):
                limit=float(name.split('<')[1].split()[0])
                q=None if row is None else number(row['cavity_mass_msun'])
                cell.update(kind='upper_limit',limit_msun=limit,qstars=q,printed=target['cavity_mass_msun'],
                            passed=None if q is None else q<limit)
            else:
                column,factor,history=OBSERVABLES[name]
                ball=float(target[column])*factor;q=None if row is None else number(row[history])
                rnd=rounding(target[column]);eff=columns[name]['tau_model']+spread(plan,name)+rnd
                d=None if q is None or q<=0 else abs(math.log(q/ball))
                cell.update(kind='ratio',ball=ball,printed=target[column],qstars=q,deviation=d,effective_tolerance=eff,
                            tau_model=columns[name]['tau_model'],realization_spread=spread(plan,name),rounding=rnd,
                            strict_tolerance=strict+spread(plan,name)+rnd,
                            passed=None if d is None else d<=eff,strict_passed=None if d is None else d<=strict+spread(plan,name)+rnd)
            if row is not None:
                fired=[k for k in report_only if row.get('screen_'+k+'_violated')=='True']
                # Always-enforced screens (plan replication_mode.always_enforced). Accepted steps cannot
                # violate them (the step is rejected); the start row (row 0) can, and the run then ends
                # there (a model_domain_exit of row 0, or a numerical_failure for a nonfinite state). Such a
                # row is not a model inside the domain: its cells are not evaluable (passed None), so the
                # verdict can be neither PASS nor, through these cells, FAIL (plan verdict_rules: FAIL needs
                # a cell outside tolerance; an exit the plan does not define leaves the verdict INCONCLUSIVE).
                enforced=[k for k in enforceable if row.get('screen_'+k+'_violated')=='True']
                cell['qstars_production_valid']=not fired and not enforced;cell['report_only_screens_fired']=fired
                cell['enforced_screens_fired']=enforced
                if enforced:
                    cell['not_evaluable']='always-enforced screen violated at this row: '+', '.join(enforced)
                    cell.update(passed=None,strict_passed=None) if cell['kind']=='ratio' else cell.update(passed=None)
                cell['opacity_clamped_mass_msun']=number(row.get('opacity_clamped_mass_msun'))
            else:cell['missing']=True
            out.append(cell)
    return out

def clamp_check(targets,rows,history=()):
    """Plan replication_mode opacity_clamp: clamped mass <= 0.01 Msun at every compared row (row 0, landings).

    valid is False when a compared row exceeds the limit (the comparison is
    then invalid and S4 is required), None when a compared row does not
    report it, else True. S4 (strict opacity domain) runs whenever any
    accepted row of the run clamped a cell.
    """
    values={}
    for index,target in enumerate(targets[:5]):
        key=0 if index==0 else float(target['bh_mass_msun']);row=rows.get(key)
        values[key]=None if row is None else number(row.get('opacity_clamped_mass_msun'))
    offending=[k for k,v in values.items() if v is not None and v>CLAMP_LIMIT_MSUN]
    unreported=[k for k,v in values.items() if v is None]
    clamped=any((number(r.get('opacity_clamped_mass_msun')) or 0.)>0 for r in history)
    return dict(limit_msun=CLAMP_LIMIT_MSUN,values=values,offending_rows=offending,unreported_rows=unreported,
                valid=False if offending else None if unreported else True,any_clamped_row=clamped,
                S4_required=clamped,follow_up='comparison invalid: run S4 (strict opacity domain)' if offending else None)

def _report_only(plan):
    from qstars.validity import POLICIES
    return POLICIES[plan['replication_mode']['mode']]['report_only']

def _enforceable(plan):
    """Screens the plan policy enforces wherever they apply (a history row flags a screen as violated only
    when it applies): every 'screen' of the registry that is neither report-only nor replaced."""
    from qstars.validity import POLICIES,SCREENS
    rules=POLICIES[plan['replication_mode']['mode']]
    return [k for k,(kind,_) in SCREENS.items() if kind=='screen' and k not in rules['report_only'] and k not in rules['replaced']]

def all3(values):
    """Three-valued AND: False if any check failed, None if any could not be evaluated, else True."""
    values=[None if v is None else bool(v) for v in values]
    return False if any(v is False for v in values) else None if any(v is None for v in values) else True

def qualitative(plan,rows,history,profiles):
    """Plan qualitative_checks Q1-Q8 (required ones decide; the rest are reported).

    Required checks are three-valued: passed is False when an evaluated
    epoch fails, else None when a needed row or profile is missing (listed
    in 'missing'), else True. A profile without convective cells fails Q5.
    """
    accepted=[r for r in history if r['event'] in ('step','landing')]
    mdot=[number(r['mdot_bh_msun_per_year']) for r in accepted];known=[x for x in mdot if x is not None]
    epoch=lambda k:rows.get(k)
    def evaluate(keys,value,test):
        values={k:value(k) for k in keys if value(k) is not None}
        return values,[k for k in keys if k not in values],all3(test(k,values[k]) if k in values else None for k in keys)
    out={}
    out['Q1']=dict(required=True,minimum=min(known,default=None),maximum=max(known,default=None),missing_steps=len(mdot)-len(known),
                   passed=all3(None if x is None else 1.6e-4<=x<=4e-4 for x in mdot) if mdot else None)
    row_value=lambda column:lambda k:number(epoch(k)[column]) if epoch(k) else None
    inv,missing,passed=evaluate((0,100.,200.,500.,1000.),row_value('max_density_inversion'),lambda k,v:v<=1e-3 if k==0 else v>1e-3)
    out['Q2']=dict(required=True,values=inv,missing=missing,passed=passed)
    stab,missing,passed=evaluate((0,100.,200.,500.,1000.),row_value('stability_ratio'),lambda k,v:v>0)
    last=number(history[-1]['stability_ratio']) if history else None
    out['Q3']=dict(required=True,values=stab,last=last,missing=missing+(['last'] if last is None else []),
                   passed=all3([passed,None if last is None else last>0]))
    present=lambda key:lambda k:profiles[k][key] if k in profiles else None
    above,missing,passed=evaluate((200.,500.,1000.),present('mass_above_outer_convective_boundary_msun'),lambda k,v:v<100)
    out['Q4']=dict(required=True,values=above,missing=missing,passed=passed)
    g1={k:profiles[k]['convective_gamma1'] for k in (100.,200.,500.,1000.) if k in profiles}
    out['Q5']=dict(required=True,values=g1,missing=[k for k in (100.,200.,500.,1000.) if k not in g1],
                   no_convective_cells=[k for k,v in g1.items() if v is None],
                   passed=all3((g1[k] is not None and 1.33<=g1[k]<=1.36) if k in g1 else None for k in (100.,200.,500.,1000.)))
    t=[number(r['time_years']) for r in history];R=[number(r['radius_rsun']) for r in history]
    rate=[abs(R[i+1]-R[i])/(t[i+1]-t[i]) for i in range(len(t)-1) if t[i+1]>t[i]]
    out['Q6']=dict(required=False,max_abs_dRdt_rsun_per_year=max(rate,default=None))
    out['Q7']=dict(required=False,values={k:abs(number(epoch(k)['luminosity_over_law'])-1) for k in (0,100.,200.,500.,1000.) if epoch(k)})
    if accepted:
        m=[number(r['bh_mass_msun']) for r in accepted]
        late=[i for i,x in enumerate(m) if x>=200]
        pick=lambda key,f:m[late[f(range(len(late)),key=lambda j:number(accepted[late[j]][key]))]] if late else None
        out['Q8']=dict(required=False,early_mdot_by_bh_mass=[(m[i],mdot[i]) for i in range(len(m)) if m[i]<=150][::max(1,len(m)//200)],
                       late_mdot_max_bh_mass=pick('mdot_bh_msun_per_year',max),late_teff_min_bh_mass=pick('effective_temperature_K',min),
                       late_radius_max_bh_mass=pick('radius_rsun',max))
    else:out['Q8']=dict(required=False)
    return out

FOLD_AGREEMENT=.02  # plan endpoint.classification_protocol[2]: N = 400 vs N = 800, |ln(M_fold ratio)|

def endpoint(plan,track,track_n800=None,n800_missing=None):
    """Plan endpoint outcome from the C0 track summary (track.json) and C2's (N = 800).

    The fold is confirmed only when C0 and C2 each confirm one (fold.py, within
    2% of their terminal M_BH) and the two fold masses agree to 2%; a C0 fold
    without that agreement is a numerical_failure (never relabelled).
    n800_missing: why no C2 summary is given (compare: C2 made with another code).
    """
    ball=plan['endpoint']['ball']
    if track is None:return dict(outcome='INCONCLUSIVE',reason='no C0 track summary')
    if track.get('final_bh_mass_msun') is None or track.get('final_time_years') is None:  # the run raised (track.json of a failed run)
        return dict(outcome='INCONCLUSIVE',label=track.get('label'),reason=f"termination classified as {track.get('label')}: {track.get('reason')}")
    m,t=track['final_bh_mass_msun'],track['final_time_years']/1e6
    dm,dt=abs(math.log(m/ball['M_BH_msun'])),abs(math.log(t/ball['t_myr'])) if t>0 else math.inf
    folded=lambda x:x is not None and x.get('label')=='model_domain_exit' and bool((x.get('fold') or {}).get('fold_confirmed'))
    fold400,fold800=folded(track),folded(track_n800)
    masses=[(x.get('fold') or {}).get('bh_mass') if x else None for x in (track,track_n800)]
    agreement=abs(math.log(masses[0]/masses[1])) if fold400 and fold800 and all(masses) else None
    confirmed=agreement is not None and agreement<=FOLD_AGREEMENT
    label='model_domain_exit' if confirmed else 'numerical_failure' if fold400 else track.get('label')
    out=dict(label=label,c0_label=track.get('label'),final_bh_mass_msun=m,final_time_myr=t,dln_bh_mass=dm,dln_time=dt,
             fold_confirmed_n400=fold400,fold_confirmed_n800=fold800,fold_bh_mass_g_n400=masses[0],fold_bh_mass_g_n800=masses[1],
             dln_fold_n400_n800=agreement,fold_agreement_tolerance=FOLD_AGREEMENT,fold_confirmed=confirmed,
             M_BH_end_passed=dm<=.1,t_end_passed=dt<=.1)
    if confirmed:out['outcome']='PASS' if dm<=.1 and dt<=.1 else 'FAIL'
    elif m>1.1*ball['M_BH_msun']:out['outcome']='FAIL';out['reason']='track continues past 1.1 x 1194 Msun with no confirmed fold'
    else:
        out['outcome']='INCONCLUSIVE'
        if fold400:
            why=(n800_missing or 'no N = 800 (C2) track') if track_n800 is None else 'C2 confirms no fold' if not fold800 else \
                f'C0 and C2 folds differ by {agreement:.3g} in ln M_BH (> {FOLD_AGREEMENT:g})'
            out['reason']=f'termination classified as numerical_failure: C0 fold not confirmed at N = 800 ({why})'
        else:out['reason']=f"termination classified as {track.get('label')}"
    return out

# ---------------------------------------------------------------- campaign checks

def compared_values(plan,rows):
    """{(row, observable): value} for every compared Table 1 observable present in rows."""
    out={}
    for key,row in rows.items():
        names=plan['comparison']['cells']['row_0' if key==0 else 'rows_1_to_4']
        for name in names:
            if name in OBSERVABLES and row is not None:
                v=number(row[OBSERVABLES[name][2]])
                if v is not None and v>0:out[(key,name)]=v
    return out

def pair(plan,a,b,profiles_a,profiles_b,threshold):
    va,vb=compared_values(plan,a),compared_values(plan,b)
    common=sorted(set(va)&set(vb),key=str)
    d={f'{k[0]}:{k[1]}':abs(math.log(va[k]/vb[k])) for k in common}
    norms={}
    for mass in sorted(set(profiles_a)&set(profiles_b)):
        for key in ('ln_rho','ln_T','L_over_LBH','ln_r'):norms[f'{mass}:{key}']=profile_norm(profiles_a[mass],profiles_b[mass],key)
    missing=sorted(f'{k[0]}:{k[1]}' for k in set(va)^set(vb))
    worst=max(list(d.values())+list(norms.values()),default=None)
    return dict(threshold=threshold,observables=d,profile_norms=norms,missing=missing,max=worst,
                passed=bool(d) and not missing and worst is not None and worst<=threshold)

def richardson(plan,coarse,base,fine):
    """C0 error per cell from three runs refined by 2, against min(0.01, tau_model/3).

    Monotone differences (e1 e2 > 0, |e1| > |e2|): Richardson, |e2| 2^p/(2^p - 1)
    with p = log2(|e1|/|e2|) > 0. Otherwise (the plan only reports order and
    monotonicity: 'reported; expected >= 0.8 and monotone') the spread of the
    three solutions, |e1| + |e2|, bounds the error; 'estimator' says which.
    """
    v1,v0,v2=(compared_values(plan,x) for x in (coarse,base,fine))
    out={}
    for k in sorted(set(v1)&set(v0)&set(v2),key=str):
        e1,e2=math.log(v1[k]/v0[k]),math.log(v0[k]/v2[k])
        tau=plan['comparison']['columns'][k[1]]['tau_model'];limit=min(.01,tau/3)
        monotone=e1*e2>0 and abs(e1)>abs(e2)
        p=math.log2(abs(e1/e2)) if e1!=0 and e2!=0 else None
        if monotone:error,estimator=abs(e2)*2**p/(2**p-1),'richardson'
        else:error,estimator=abs(e1)+abs(e2),'spread'
        out[f'{k[0]}:{k[1]}']=dict(order=p,error=error,estimator=estimator,limit=limit,monotone=monotone,
                                   expected_order_met=monotone and p>=.8,passed=error<=limit)
    return out

def conservation(plan,history):
    """Plan conservation thresholds at every accepted step of one run.

    reconstruction_event: the step's largest |EOS-inversion defect| of one
    event booked to remap_energy_defect (step_max_event_defect: a split/merge
    group or the fixed-mass base regrid) over the initial internal energy.
    A value a row does not report counts as a violation (missing).
    """
    c=plan['conservation'];bad=dict(mass_per_step=0,mass_cumulative=0,energy_cumulative=0,engine=0,geometry=0,
                                    reconstruction_event=0,inflow_flagged=0)
    worst=dict.fromkeys(bad,0.);missing=dict.fromkeys(bad,0);steps=0
    magnitude=lambda key:None if number(r.get(key)) is None else abs(number(r[key]))
    for r in history:
        if r['event'] not in ('step','landing'):continue
        steps+=1
        scale=number(r.get('ledger_initial_energy_scale'));event=magnitude('step_max_event_defect')
        engine=magnitude('engine_conversion_residual_erg')
        checks=dict(mass_per_step=(magnitude('step_mass_residual_change'),c['mass_per_step_relative']),
                    mass_cumulative=(magnitude('mass_residual_relative'),c['mass_cumulative_relative']),
                    energy_cumulative=(magnitude('energy_residual_relative'),c['energy_cumulative_over_initial_internal_energy']),
                    engine=(None if engine is None else engine/max(number(r['ledger_central_generated']),1.),c['engine_identity_relative']),
                    geometry=(number(r['screen_volume_geometry']) or 0.,c['geometry']),
                    reconstruction_event=(None if event is None or not scale else event/scale,c['reconstruction_defect_per_event_relative']),
                    inflow_flagged=(magnitude('step_inflow_heating_over_bh_energy') or 0.,c['inflow_numerical_heating_over_L_BH']['flag_if_above']))
        for k,(v,limit) in checks.items():
            if v is None:missing[k]+=1;bad[k]+=1;continue
            worst[k]=max(worst[k],v)
            if v>limit:bad[k]+=1
    required=('mass_per_step','mass_cumulative','energy_cumulative','engine','geometry','reconstruction_event')
    return dict(accepted_steps=steps,violations=bad,missing=missing,worst=worst,passed=steps>0 and all(bad[k]==0 for k in required))

def restart_check(c0,segment,label):
    """Plan restart: the segment's rows equal C0's bitwise and its end checkpoint equals C0's segment_end."""
    from qstars.quasistatic import load_checkpoint
    a={r['step']:r for r in read_history(c0)};b=read_history(segment)
    start=int(b[0]['step']);rows_equal=all(r==a.get(r['step']) for r in b[1:])
    end=int(b[-1]['step']);reference=Path(c0)/f'segment_end_{end:07d}.npz'
    if not reference.exists() and end==int(read_history(c0)[-1]['step']):reference=Path(c0)/'checkpoint.npz'  # ended by termination
    checkpoint_equal=False
    if reference.exists():
        x,_=load_checkpoint(reference);y,_=load_checkpoint(Path(segment)/'checkpoint.npz')
        checkpoint_equal=all(np.array_equal(getattr(x,f),getattr(y,f)) for f in ('r','T','dm','v','cell_volume','radius_compensation','face_luminosity')) \
            and all(getattr(x,f)==getattr(y,f) for f in ('bh_mass','cavity_mass','photosphere_density','time','steps')) \
            and json.dumps(x.aux,sort_keys=True)==json.dumps(y.aux,sort_keys=True) and json.dumps(x.ledger,sort_keys=True)==json.dumps(y.ledger,sort_keys=True)
    seg=b[1:]
    def rejected(path,steps=None):
        events=Path(path)/'events.jsonl'
        if not events.exists():return []
        out=[json.loads(x) for x in open(events) if json.loads(x)['event']=='rejected']
        return [{k:e[k] for k in ('step','dt_years','message')} for e in out if steps is None or e['step'] in steps]
    replayed=rejected(segment)  # the rejections this segment itself replayed
    term=Path(segment)/'termination.json'
    stopped=not term.exists() or json.loads(term.read_text()).get('status')=='control_stop'
    attempted=range(start,end if stopped else end+1)  # a control stop makes no attempt at its last step
    crossed=dict(base_slice=any(number(r['step_slice_mass']) for r in seg),
                 split_merge=any(int(number(r['step_splits']) or 0)+int(number(r['step_merges']) or 0) for r in seg),
                 rejected_step=bool(replayed),landing=any(r['event']=='landing' for r in seg))
    return dict(label=label,start_step=start,end_step=end,rows=len(seg),rows_bitwise_equal=rows_equal,
                end_checkpoint=str(reference),end_checkpoint_bitwise_equal=checkpoint_equal,crossed=crossed,
                rejections_replayed=len(replayed),rejections_equal_c0=replayed==rejected(c0,attempted),
                passed=bool(seg) and rows_equal and checkpoint_equal)

def must_cross(restarts):
    """Plan restart.must_cross over the segments together; the landing within R2's segment."""
    crossed={k:any(r['crossed'][k] for r in restarts.values()) for k in MUST_CROSS}
    landing=bool(restarts.get('R2',{}).get('crossed',{}).get('landing'))
    unmet=[k for k in MUST_CROSS if not crossed[k]]+([] if landing else ['landing (restart point 2)'])
    return dict(crossed=crossed,R2_landing=landing,unmet=unmet,passed=not unmet)

def verdict(cells_out,end,quality,convergence_ok,conservation_ok,restarts_ok,plan_ok=True,clamp_ok=True):
    """Plan verdict_rules, in order: FAIL, INCONCLUSIVE, PASS.

    A required qualitative check fails the verdict only when it is False; one
    that could not be evaluated (None) leaves it INCONCLUSIVE. clamp_ok is
    clamp_check()['valid']: False (a compared row above the opacity-clamp
    limit) invalidates the comparison, so the verdict is INCONCLUSIVE whatever
    else holds; None (not reported) blocks PASS.
    """
    ratio=[c for c in cells_out if c.get('passed') is not None]
    all_cells=len(ratio)==len(cells_out) and all(c['passed'] for c in ratio)
    failed_cells=[c for c in cells_out if c.get('passed') is False]
    required=[q['passed'] for q in quality.values() if q.get('required')]
    if clamp_ok is False:return 'INCONCLUSIVE'
    if end['outcome']=='FAIL' or any(v is False for v in required) or (failed_cells and convergence_ok):
        return 'FAIL'
    if not (convergence_ok and conservation_ok and restarts_ok and plan_ok and clamp_ok) or end['outcome']=='INCONCLUSIVE' \
            or not all_cells or any(v is None for v in required):
        return 'INCONCLUSIVE'
    return 'PASS'

def compare(plan,runs,restarts,out_dir,registration):
    """comparison.json and campaign.json from finished run directories (runs: id -> dir)."""
    from qstars import io
    check_mapped_texts(plan)
    targets=table_rows(ROOT/plan['sources']['targets']['primary'])
    def load(run_id):
        if run_id not in runs:return None
        path=Path(runs[run_id]);history=read_history(path)
        rows={0:history[0] if history and history[0]['event']=='start' else None,**landed_rows(history)}
        profiles={m:profile_quantities(path/f'profile_landing_{m:g}msun.npz') for m in landed_rows(history)
                  if (path/f'profile_landing_{m:g}msun.npz').exists()}
        track=json.loads((path/'track.json').read_text()) if (path/'track.json').exists() else None
        return dict(path=str(path),history=history,rows=rows,profiles=profiles,track=track)
    data={k:load(k) for k in runs}
    c0=data.get('C0')
    comparison=dict(plan_sha256=registration['plan_sha256'],runs={k:v['path'] for k,v in data.items() if v})
    # One code: the comparisons between runs (convergence pairs, Richardson, the N = 400 vs 800 fold
    # agreement, restarts) hold only between runs made with one established code (source_consistency);
    # C0 and its row-0 relax run must share it for any verdict but INCONCLUSIVE.
    code=source_consistency({**runs,**{f'restart {k}':v for k,v in restarts.items()}})
    if c0:
        cell_list=cells(plan,targets,c0['rows'])
        quality=qualitative(plan,c0['rows'],c0['history'],c0['profiles'])
        c2_same='C2' in runs and same_code(code,'C0','C2')
        end=endpoint(plan,c0['track'],(data.get('C2') or {}).get('track') if c2_same else None,
                     None if c2_same or 'C2' not in runs else 'C2 was not made with the code of C0 (source_consistency)')
        clamp=clamp_check(targets,c0['rows'],c0['history'])
        warnings=sorted({w for r in c0['history'] for w in r['screen_warnings'].split(';') if w})
    else:cell_list,quality,end,clamp,warnings=[],{},endpoint(plan,None),dict(valid=None),[]
    comparison.update(cells=cell_list,cell_count=len(cell_list),cells_passed=sum(1 for c in cell_list if c.get('passed')),
                      failed_cells=[c for c in cell_list if c.get('passed') is False],qualitative_checks=quality,endpoint=end,
                      opacity_clamp=clamp,
                      report_only_screens_fired=warnings,
                      label=plan['replication_mode']['label'] if warnings else 'replica within the QStarS production screens')
    t=plan['convergence_campaign']['thresholds'];campaign=dict(plan_sha256=registration['plan_sha256'])
    def check(name,a,b,threshold):
        if not (data.get(a) and data.get(b)):return dict(threshold=threshold,passed=False,missing=[x for x in (a,b) if not data.get(x)])
        return pair(plan,data[a]['rows'],data[b]['rows'],data[a]['profiles'],data[b]['profiles'],threshold)
    campaign['C0_vs_C2']=check('C0_vs_C2','C0','C2',t['C0_vs_C2_and_C0_vs_C4'])
    campaign['C0_vs_C4']=check('C0_vs_C4','C0','C4',t['C0_vs_C2_and_C0_vs_C4'])
    campaign['C0_vs_C5']=check('C0_vs_C5','C0','C5',t['C0_vs_C5'])
    campaign['C0_vs_C7']=check('C0_vs_C7','C0','C7',t['C0_vs_C7'])
    campaign['C6']=dict(reported=check('C0_vs_C6','C0','C6',math.inf),
                        note='plan: reported; convergence of C6 with N needs C6 at more than one N, which the plan does not register')
    for name,(a,b,c) in (('richardson_N',('C1','C0','C2')),('richardson_dt',('C3','C0','C4'))):
        campaign[name]=richardson(plan,data[a]['rows'],data[b]['rows'],data[c]['rows']) if all(data.get(x) for x in (a,b,c)) else None
    ends=[data[x]['track']['final_bh_mass_msun'] if data.get(x) and data[x]['track'] else None for x in ('C0','C2')]
    campaign['endpoint_N400_vs_N800']=dict(threshold=t['endpoint_M_BH_between_N400_and_N800'],values=ends,
        dln=abs(math.log(ends[0]/ends[1])) if all(ends) else None,
        passed=all(ends) and abs(math.log(ends[0]/ends[1]))<=t['endpoint_M_BH_between_N400_and_N800'])
    campaign['escalation_N1600']=bool(campaign['C0_vs_C2'].get('max') is not None and campaign['C0_vs_C2']['max']>t['C0_vs_C2_and_C0_vs_C4'])
    campaign['conservation']={k:conservation(plan,v['history']) for k,v in data.items() if v and k.startswith('C')}
    campaign['restarts']={label:restart_check(runs['C0'],path,label) for label,path in restarts.items()} if 'C0' in runs else {}
    campaign['restart_must_cross']=must_cross(campaign['restarts'])
    rich=[x for name in ('richardson_N','richardson_dt') for x in (campaign[name] or {None:dict(passed=False)}).values()]
    convergence_ok=all(campaign[k]['passed'] for k in ('C0_vs_C2','C0_vs_C4','C0_vs_C5','C0_vs_C7','endpoint_N400_vs_N800')) \
        and all(x['passed'] for x in rich) and code['one_code']
    conservation_ok=all(f'C{i}' in campaign['conservation'] and campaign['conservation'][f'C{i}']['passed'] for i in range(8))
    restarts_ok=len(campaign['restarts'])==3 and all(r['passed'] for r in campaign['restarts'].values()) \
        and campaign['restart_must_cross']['passed'] and code['one_code']
    campaign.update(convergence_met=convergence_ok,conservation_met=conservation_ok,restarts_bitwise=restarts_ok,
                    opacity_clamp_valid=clamp['valid'],one_code=code['one_code'],c0_one_code=code['c0_one_code'])
    campaign['verdict']=verdict(cell_list,end,quality,convergence_ok,conservation_ok,restarts_ok,clamp_ok=clamp['valid']) \
        if c0 and code['c0_one_code'] else 'INCONCLUSIVE'
    campaign['label']=comparison['label']
    campaign['source_consistency']=code
    campaign['rejection_causes']={k:rejection_causes(v) for k,v in {**runs,**{f'restart {k}':v for k,v in restarts.items()}}.items()}
    comparison['verdict']=campaign['verdict']
    io.json_write(Path(out_dir)/'comparison.json',_clean(comparison))
    io.json_write(Path(out_dir)/'campaign.json',_clean(campaign))
    return comparison,campaign

def _clean(x):
    from qstars.quasistatic import _json
    return _json(x)

# The Ball import closure: every qstars module a driver process imports (relax, smoke, track with
# the fold classification, compare). Its digest (code_closure, from the sources as imported) is the
# code of a run; code_closure.json also lists any qstars module a process loaded outside it, and such a
# run has no established code (source_consistency).
BALL_MODULES=('__init__','cavity','config','constants','controls','engine','eos','fold','io','ledger','mass_transport',
              'mesa_microphysics','mesh','mlt','nonlinear','numerics','photosphere','quasistatic','state','stepcontrol',
              'structure','transport','validity','winds')

def import_ball_modules():
    """Import the whole Ball closure up front (main), so every driver process holds the same modules."""
    import importlib
    for module in BALL_MODULES:
        if module!='__init__':importlib.import_module(f'qstars.{module}')

def code_closure():
    """sha256 of the Ball import closure as this process imported it (qstars.io keeps the package
    sources of the process), and the qstars modules loaded outside BALL_MODULES."""
    from qstars import io
    sources=io._PROCESS_SOURCES;digest=hashlib.sha256();missing=[]
    for module in BALL_MODULES:
        name='qstars/__init__.py' if module=='__init__' else f'qstars/{module}.py'
        if name not in sources:missing.append(name)
        digest.update(name.encode()+b'\0'+sources.get(name,b'')+b'\0')
    loaded=sorted(m.split('.',1)[1] for m in sys.modules if m.startswith('qstars.'))
    return dict(modules=list(BALL_MODULES),sha256=digest.hexdigest(),missing=missing,
                loaded_outside_closure=[m for m in loaded if m not in BALL_MODULES])

def finish(path):
    """io.finish_provenance and code_closure.json (the code this process ran); returns source_snapshot_verified."""
    from qstars import io
    verified=io.finish_provenance(path)
    io.json_write(Path(path)/'code_closure.json',code_closure())
    return verified

def source_consistency(runs):
    """The code each run directory was made with, and whether one code made them all (enforced by compare).

    Per directory (tracks, restart segments, and the relax run of each
    track's row-0 model): the digest of the Ball import closure
    (code_closure.json), the sha256 of the archived driver.py and of the
    MESA library (run.json), source_snapshot_verified (no package file
    changed while the run was in progress), and the qstars package digest
    (every file under src/qstars, reported only: other workflows edit
    modules outside the closure). A run's code is established when all four
    are recorded, the snapshot is verified and the process loaded no qstars
    module outside the closure. one_code: every listed directory has
    established code and they all share it. c0_one_code: the same for C0
    and the relax run of its row-0 model.
    """
    out={}
    def entry(path):
        path=Path(path);meta=_read_json(path/'metadata.json') or {};run=_read_json(path/'run.json') or {}
        closure=_read_json(path/'code_closure.json') or {};driver=path/'driver.py'
        e=dict(path=str(path),closure_sha256=closure.get('sha256'),loaded_outside_closure=closure.get('loaded_outside_closure'),
               driver_sha256=sha256(driver) if driver.is_file() else None,library_sha256=run.get('library_sha256'),
               source_snapshot_verified=meta.get('source_snapshot_verified'),package_sha256=meta.get('process_package_sha256'))
        e['established']=bool(e['closure_sha256'] and e['driver_sha256'] and e['library_sha256'] and e['source_snapshot_verified'] is True
                              and e['loaded_outside_closure']==[] and not closure.get('missing'))
        e['code']=[e['closure_sha256'],e['driver_sha256'],e['library_sha256']] if e['established'] else None
        return e
    for label,path in runs.items():
        out[label]=entry(path)
        row0=(_read_json(Path(path)/'run.json') or {}).get('row0')
        if row0:out[label+' row0 relax']=entry(Path(row0).parent)
    def one(labels):
        codes=[out[k]['code'] for k in labels if k in out]
        return bool(codes) and all(c is not None for c in codes) and all(c==codes[0] for c in codes)
    packages={e['package_sha256'] for e in out.values()};drivers={e['driver_sha256'] for e in out.values()}
    return dict(runs=out,one_code=one(list(out)),c0_one_code=one(['C0','C0 row0 relax']) if 'C0' in out else False,
                one_package=len(packages)==1 and None not in packages,one_driver=len(drivers)==1 and None not in drivers,
                all_snapshots_verified=bool(out) and all(e['source_snapshot_verified'] is True for e in out.values()))

def same_code(consistency,a,b):
    """Runs a and b (with their row-0 relax runs) were made with one established code."""
    labels=[k for k in (a,b,a+' row0 relax',b+' row0 relax') if k in consistency['runs']]
    codes=[consistency['runs'][k]['code'] for k in labels]
    return bool(codes) and all(c is not None for c in codes) and all(c==codes[0] for c in codes)

# Rejection causes (events.jsonl), reported per run: the plan's reject rule names delta_actual > 2
# delta_target, a Newton failure and more than 3 predictor passes; the step also rejects on the per-step
# conservation budgets, enforced screens and errors before the solve. A Newton failure whose final
# residual is at or below the roundoff-exit residual (1e-10) is counted as at the tolerance floor.
def rejection_causes(path):
    import re
    events=Path(path)/'events.jsonl'
    out=dict(total=0,causes={},controller_metric={},newton_failures_at_tolerance_floor=0,newton_failure_residuals=[])
    if not events.is_file():return out
    for line in events.read_text().splitlines():
        e=json.loads(line)
        if e.get('event')!='rejected':continue
        out['total']+=1;message=e.get('message') or ''
        if e.get('kind')=='screen':cause='enforced_screen'
        elif message.startswith('delta_actual'):
            cause='controller_metric';metric=message.rsplit('(',1)[-1].rstrip(')') if message.endswith(')') else 'unknown'
            out['controller_metric'][metric]=out['controller_metric'].get(metric,0)+1
        elif message.startswith('Quasi-static solve failed'):
            cause='newton_failure_domain' if e.get('domain_exit') else 'newton_failure'
            found=re.search(r'residual ([0-9.eE+-]+|inf|nan)\)',message)
            residual=float(found.group(1)) if found else None
            out['newton_failure_residuals'].append(residual)
            if residual is not None and residual<=1e-10:out['newton_failures_at_tolerance_floor']+=1
        elif message.startswith('Per-step conservation budget'):cause='conservation_budget'
        elif message.startswith('Sink correction'):cause='predictor_passes'
        elif message.startswith('Landing'):cause='landing'
        elif re.match(r'P\d ',message):cause='domain_exit_before_solve' if e.get('domain_exit') else 'error_before_solve'
        else:cause='other'
        out['causes'][cause]=out['causes'].get(cause,0)+1
    return out

# ---------------------------------------------------------------- subcommands

def common(a):
    registration=verify_registration(a.plan,a.registration)
    plan=json.loads(Path(a.plan).read_text())
    amendments=[verify_amendment(p,registration['plan_sha256']) for p in a.amendment]
    library=Path(a.library).resolve()
    mesa_dir=a.mesa_dir or (json.loads((library.parent/'build.json').read_text())['mesa_dir'] if library.is_file() else '')
    parent=Path(a.output) if a.output else ROOT/'benchmarks/results'/(datetime.now(timezone.utc).strftime('%Y-%m-%d')+'-ball-track')
    check_track_config(plan)
    return registration,plan,amendments,library,mesa_dir,parent

ROW0_EXIT_LABEL='model_domain_exit of row 0 (enforced screen)'
RESTART_ROW_EXIT_LABEL='model_domain_exit at the restart row (enforced screen)'

def failed_report(exc):
    """The record of a run that raised instead of returning its termination report (never lost, never relabelled)."""
    import traceback
    return dict(status='numerical_failure',reason=f'driver error: {type(exc).__name__}: {exc}',accepted_steps_this_run=0,attempts=0,
                rejections=0,final_step=None,final_time_years=None,final_bh_mass_msun=None,landed_msun=[],last_rejection=None,
                traceback=traceback.format_exc(),v1_validated=False)

def track(a,plan,registration,amendments,library,mesa_dir,parent,run_id=None,row0=None,restart=None,segment=None,fallbacks=()):
    """One registered track (or restart segment) into its own run directory; returns (path, summary).

    track.json is always written: an error raised by the run or by the fold
    classification is recorded as a numerical_failure of that stage (the run
    directory, its termination.json and events.jsonl are kept), never lost
    and never relabelled. A start row that violates an always-enforced
    screen is a model_domain_exit of row 0 and is not passed to the fold
    protocol (no termination trigger fired).
    """
    from qstars import fold,io,quasistatic
    MSUN=1.98847e33
    if restart is not None:
        state,config=quasistatic.load_checkpoint(restart);run_id=config['replication']['run']
        if config['replication']['amendments']!=[r['path'] for r,_ in amendments]:
            raise Refused('Refusing to run: the restart checkpoint was made with other amendments')
        if fallbacks_of(config)!=sorted(fallbacks):
            raise Refused(f'Refusing to run: the restart checkpoint was made with fallbacks {fallbacks_of(config)}, not {sorted(fallbacks)}')
        kind='restart'
    else:
        config=apply_fallbacks(build_config(plan,run_id,library,mesa_dir,amendments),fallbacks)
        candidates=row0 if isinstance(row0,(list,tuple)) else [row0]
        reasons=[];state,row0,_=row0_pick(config,candidates,reasons,row0_provenance(parent,registration,amendments))
        if state is None:
            raise Refused(f'Refusing to run: no t = 0 row-0 model relaxed with the configuration of {run_id} among {candidates} '
                          f'(relax --run {run_id}): '+'; '.join(reasons))
        kind='track'
    eos,opacity=physics(config)
    path=run_directory(parent,kind,run_id,config,registration,amendments,
                       dict(row0=str(row0) if row0 else None,row0_relax=row0_relax(row0) if row0 else None,
                            restart=str(restart) if restart else None,restart_sha256=sha256(restart) if restart else None,
                            segment=segment,fallbacks=fallbacks_of(config)))
    options=track_options(plan,config)
    summary=dict(run=run_id,kind=kind)
    try:
        try:
            if restart is not None:
                state_out,report=quasistatic.run(state,eos,opacity,config,path,stop_after_steps=segment,**options)
            else:
                state_out,report=quasistatic.run(state,eos,opacity,config,path,restart_masses=[20*MSUN,900*MSUN],
                                                 restart_landing=500*MSUN,segment_steps=SEGMENT_STEPS,**options)
        except Exception as exc:  # quasistatic.run records its own failures; this one escaped it
            state_out,report=None,failed_report(exc)
            io.json_write(path/'driver_error.json',_clean(report))
        summary.update(status=report['status'],reason=report['reason'],final_bh_mass_msun=report['final_bh_mass_msun'],
                       final_time_years=report['final_time_years'],final_step=report['final_step'],landed_msun=report['landed_msun'],
                       attempts=report['attempts'],rejections=report['rejections'],fallbacks=fallbacks_of(config))
        if report.get('start_row_screen_errors'):
            # No step was attempted. A nonfinite start state is numerical (quasistatic.run status
            # numerical_failure); any other enforced screen is an exit of the start row: of row 0 for a
            # track (an applicability exit the plan does not define; endpoint INCONCLUSIVE), of the
            # restart row for a restart segment.
            label=('numerical_failure' if report['status']=='numerical_failure' else
                   RESTART_ROW_EXIT_LABEL if restart is not None else ROW0_EXIT_LABEL)
            summary.update(label=label,label_reason=report['reason'],start_row_screen_errors=report['start_row_screen_errors'])
        elif report['status'] in ('termination_trigger','model_domain_exit') and restart is None and state_out is not None:
            # Both are the plan's termination trigger (dt < t_dyn or 12 consecutive retries): the fold protocol
            # labels the end; retries exhausted by enforced screens are recorded next to that label.
            try:
                fold_out=fold.classify(state_out,eos,opacity,config)
                io.json_write(path/'fold.json',_clean(fold_out))
                summary.update(label=fold_out['label'],label_reason=fold_out['reason'],fold=dict(fold_confirmed=fold_out['fold_confirmed'],
                               **{k:fold_out['fold'][k] for k in ('bh_mass','relative_offset','sqrt_exponent','det_sign_change','within_window')}
                               if fold_out['fold'] else {}))
            except Exception as exc:  # a failed classification is a numerical failure of the protocol, never a fold
                import traceback
                error=dict(error=f'{type(exc).__name__}: {exc}',traceback=traceback.format_exc())
                io.json_write(path/'fold_error.json',error)
                summary.update(label='numerical_failure',label_reason=f"fold classification failed: {error['error']}",
                               fold=dict(fold_confirmed=False),fold_error=error['error'])
            if report['status']=='model_domain_exit':
                summary['enforced_screen_exit']=dict(reason=report['reason'],last_rejection=report['last_rejection'])
        else:
            summary['label']={'model_domain_exit':'model_domain_exit (enforced screen)','stop_rule':'numerical_blocker (stop rule)',
                              'completed':'completed','control_stop':'control_stop'}.get(report['status'],'numerical_failure')
    finally:
        # Written whatever happened above (KeyboardInterrupt included), so the directory always says how far it got.
        summary.setdefault('status','control_stop');summary.setdefault('label','control_stop (track summary incomplete)')
        try:summary['rejection_causes']=rejection_causes(path)
        except Exception as exc:summary['rejection_causes']=dict(error=f'{type(exc).__name__}: {exc}')
        io.json_write(path/'track.json',_clean(summary))
        finish(path)
    return path,summary

def smoke(a,plan,registration,amendments,library,mesa_dir,parent,row0,fallbacks=()):
    """Plan stop_rules.smoke_gate: the N = 200 run (C1) to M_BH = 100; returns (path, smoke.json content).

    smoke.json is always written, also when the run raised or never wrote
    its start row (history.csv): such a gate has not passed.
    """
    from qstars import io,quasistatic
    config=apply_fallbacks(build_config(plan,'C1',library,mesa_dir,amendments),fallbacks)
    if 'N = 200' not in plan['stop_rules']['smoke_gate']:raise Refused('Refusing to run: the smoke gate is not an N = 200 run')
    reasons=[];state=row0_for(config,[row0],reasons,row0_provenance(parent,registration,amendments))
    if state is None:
        raise Refused(f'Refusing to run: {row0} is not a t = 0 row-0 model of C1 (relax --run C1): '+'; '.join(reasons))
    eos,opacity=physics(config)
    path=run_directory(parent,'smoke','C1',config,registration,amendments,dict(row0=str(row0),row0_relax=row0_relax(row0),
                                                                           fallbacks=fallbacks_of(config)))
    options=track_options(plan,config)
    report=dict(status='control_stop',reason='smoke run interrupted before its report')
    try:
        try:
            _,report=quasistatic.run(state,eos,opacity,config,path,end_bh_mass=options['landings'][0],
                                     **{**options,'landings':options['landings'][:1]})
        except Exception as exc:
            report=failed_report(exc);io.json_write(path/'driver_error.json',_clean(report))
    finally:
        out=dict(smoke_summary(path,report,plan),fallbacks=fallbacks_of(config),
                 start_row_screen_errors=report.get('start_row_screen_errors') or [])
        try:out['rejection_causes']=rejection_causes(path)
        except Exception as exc:out['rejection_causes']=dict(error=f'{type(exc).__name__}: {exc}')
        io.json_write(path/'smoke.json',_clean(out));finish(path)
    return path,out

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--plan',default=str(PLAN));p.add_argument('--registration',default=str(REGISTRATION))
    p.add_argument('--amendment',action='append',default=[],help='registered amendment JSON to apply (repeatable)')
    p.add_argument('--fallback',action='append',default=[],choices=sorted(FALLBACKS),
                   help='plan stop_rules.fallbacks_in_order switched on for smoke, track and campaign (F1: pseudo-transient '
                        'continuation inside a step; recorded in run.json and the configuration)')
    p.add_argument('--library',default=os.environ.get('QSTARS_MESA_LIBRARY',str(ROOT/'backends/mesa/build_ball_abi7/libqstars_mesa.dylib')))
    p.add_argument('--mesa-dir',default=os.environ.get('MESA_DIR'),help='default: mesa_dir of the build manifest')
    p.add_argument('--output',default=None,help='parent directory; default benchmarks/results/<UTC date>-ball-track')
    sub=p.add_subparsers(dest='command',required=True)
    r=sub.add_parser('relax',help='relax the registered seeds (plan initial_model)')
    r.add_argument('--run',default='C0',help='registered run whose configuration is relaxed (C1, C0, C2: N = 200, 400, 800)')
    s=sub.add_parser('smoke',help='plan smoke gate (N = 200 to M_BH = 100)')
    s.add_argument('--row0',required=True,help='the row-0 checkpoint of a converged_unique relax run at N = 200 (relax --run C1)')
    t=sub.add_parser('track',help='one registered run, or a restart segment')
    t.add_argument('--run',default='C0');t.add_argument('--row0',help='row-0 checkpoint with the run\'s zone count')
    t.add_argument('--restart',help='checkpoint of a restart point');t.add_argument('--segment',type=int,default=SEGMENT_STEPS)
    c=sub.add_parser('campaign',help='the registered runs, the C0 restart segments, then compare')
    c.add_argument('--runs',nargs='+',default=[f'C{i}' for i in range(8)])
    c.add_argument('--row0',action='append',default=[],help='RUN=PATH row-0 checkpoint (relax --run RUN); a run without one '
                   'uses any listed row-0 model relaxed with its configuration (repeatable)')
    m=sub.add_parser('compare',help='comparison.json and campaign.json from finished runs')
    m.add_argument('--run',action='append',default=[],help='ID=DIR of a finished track (repeatable)')
    m.add_argument('--restart',action='append',default=[],help='LABEL=DIR of a restart segment of C0 (repeatable)')
    a=p.parse_args(argv)
    import_ball_modules()
    registration,plan,amendments,library,mesa_dir,parent=common(a)
    from qstars import io
    fallbacks=sorted(set(a.fallback))
    for name in fallbacks:
        if FALLBACK_TEXTS[name] not in plan['stop_rules']['fallbacks_in_order']:
            raise Refused(f'Refusing to run: fallback {name} is not in the plan stop_rules.fallbacks_in_order')
    if fallbacks and a.command in ('relax','compare'):
        p.error('--fallback applies to smoke, track and campaign: the relaxation is the registered one (bondi_step mode relax never uses F1)')
    if a.command=='relax':
        config=build_config(plan,a.run,library,mesa_dir,amendments)
        eos,opacity=physics(config)
        path=run_directory(parent,'relax',a.run,config,registration,amendments)
        try:out=relax_seeds(plan,config,eos,opacity,path,config['replication']['seed_radii_rsun'],config['numerics']['zones'])
        finally:finish(path)
        print(json.dumps({k:out[k] for k in out if k!='seeds'},default=str,indent=1))
        return 0 if out['status']=='converged_unique' else 1
    if a.command=='smoke':
        path,out=smoke(a,plan,registration,amendments,library,mesa_dir,parent,a.row0,fallbacks=fallbacks)
        print(json.dumps(out,indent=1,default=str));return 0 if out['passed'] else 1
    if a.command=='track':
        if not (a.row0 or a.restart):p.error('track needs --row0 or --restart')
        path,summary=track(a,plan,registration,amendments,library,mesa_dir,parent,run_id=a.run,row0=a.row0,
                           restart=a.restart,segment=a.segment if a.restart else None,fallbacks=fallbacks)
        print(json.dumps(summary,indent=1,default=str));return 0
    if a.command=='campaign':
        rows0=dict(x.split('=',1) for x in a.row0)
        reasons=[];provenance=row0_provenance(parent,registration,amendments)
        missing=[r for r in a.runs if row0_for(apply_fallbacks(build_config(plan,r,library,mesa_dir,amendments),fallbacks),
                                               [rows0.get(r)]+list(rows0.values()),reasons,provenance) is None]
        if missing:
            raise Refused(f'Refusing to run: no row-0 model relaxed with the configuration of {missing} (relax --run ...): '
                          +'; '.join(sorted(set(reasons))))
        index=dict(started_utc=datetime.now(timezone.utc).isoformat(),runs={},restarts={},fallbacks=fallbacks)
        runs={}
        for run_id in a.runs:
            path,summary=track(a,plan,registration,amendments,library,mesa_dir,parent,run_id=run_id,
                               row0=[rows0.get(run_id)]+list(rows0.values()),fallbacks=fallbacks)
            runs[run_id]=str(path);index['runs'][run_id]=dict(path=str(path),**summary)
        restarts={}
        if 'C0' in runs:
            for label,ckpt in RESTART_POINTS:
                source=Path(runs['C0'])/ckpt
                if source.exists():
                    path,summary=track(a,plan,registration,amendments,library,mesa_dir,parent,restart=source,segment=SEGMENT_STEPS,
                                       fallbacks=fallbacks)
                    restarts[label]=str(path);index['restarts'][label]=dict(path=str(path),**summary)
                else:index['restarts'][label]=dict(missing=str(source))
        out=io.create_directory(parent,'campaign');io.json_write(out/'campaign_runs.json',_clean(index))
        io.json_write(out/'plan_registration.json',registration)
        compare(plan,runs,restarts,out,registration);print(out);return 0
    if a.command=='compare':
        runs=dict(x.split('=',1) for x in a.run);restarts=dict(x.split('=',1) for x in a.restart)
        out=io.create_directory(parent,'compare');io.json_write(out/'plan_registration.json',registration)
        io.json_write(out/'inputs.json',dict(runs=runs,restarts=restarts,amendments=[r for r,_ in amendments]))
        comparison,campaign=compare(plan,runs,restarts,out,registration)
        print(out);print(json.dumps(dict(verdict=campaign['verdict'],cells_passed=comparison['cells_passed'],
                                         cell_count=comparison['cell_count'],endpoint=comparison['endpoint']['outcome']),indent=1))
        return 0

if __name__=='__main__':raise SystemExit(main())
