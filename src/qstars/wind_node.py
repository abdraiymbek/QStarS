"""Contract for a CMFGEN-derived self-selected wind node (evolution coupling).

A node is one converged steady CMFGEN wind at fixed base data. Duplicated
fields must satisfy exact identities; nothing is fitted or inferred. Gate
flags are the producer's own claims and are reported as such. No Gate 4
evidence verifier exists, so every node is reported unvalidated.
"""
import json,math,re
from pathlib import Path
import numpy as np
from .constants import A,G,PI,SIGMA
from .outer_interface import InterfaceState,compare_faces
from .mass_transport import DomainError

SCHEMA,SET_SCHEMA,VERSION='qstars_wind_node','qstars_wind_node_set',1
# Positive finite CGS scalars. BASE selects/guards a node; OUTPUT enters evolution or diagnostics.
BASE={'M_total_g':'total mass inside the base (point-mass gravity), g',
      'R_base_cm':'CMFGEN inner boundary = envelope outer face, cm',
      'L_base_erg_s':'lab radiative luminosity through the base, 4 pi R^2 F, erg/s',
      'rho_base':'base mass density, g/cm^3','T_base':'base gas temperature, K'}
STATE={'P_gas_base':'base gas pressure, dyn/cm^2','P_rad_base':'base lab radiation pressure P_rr (K moment), dyn/cm^2'}
OUTPUT={'mdot_g_s':'self-selected mass-loss eigenvalue, g/s','v_inf_cm_s':'terminal speed, cm/s',
        'L_emergent_erg_s':'escaping lab luminosity at the outer boundary, erg/s',
        'R_phot_cm':'photospheric radius (Rosseland tau=2/3), cm',
        'Teff_emergent_K':'(L_emergent/(4 pi sigma R_phot^2))^(1/4), K'}
FACE=('radius','enclosed_mass','density','velocity','gas_pressure','gas_internal_energy','radiation_energy',
      'radiation_flux','radiation_pressure_rr','potential','energy_gauge','gas_temperature')
# Gate claims the node producer must state: one flag per docs/nlte_backend.md release-boundary
# gate (mapping in docs/wind_node.md). Contract proposal pending Gate 4 and user approval.
# Envelope admissibility and the evolution-EOS energy zero are checked by the run, not flagged.
VALIDATION=('atomic_data_validated','radiation_force_validated','material_energy_validated',
            'nonlinear_converged','electron_charge_radiation_audits','radius_converged','frequency_converged',
            'atom_size_converged','export_invariant','sonic_regular','photosphere_resolved','escapes',
            'independent_ledger','mach_one_mlt_guard')
PROVENANCE=('cmfgen_executable_sha256','vadat_sha256','atomic_data_sha256','solution_sha256')
RTOL=1e-10  # duplicated/derived fields come from one double-precision export
ID=r'[A-Za-z0-9_.:-]{1,100}'
# Evolution EOS energy zeros. MESA's EOS blend has no single documented zero
# (FreeEOS uses ground-state H2); it stays undeclared until a bridge exists.
EVOLUTION_GAUGES={'ideal_gas_radiation':'fully_ionized_H_He_at_rest'}

def _number(x):
    return not isinstance(x,bool) and isinstance(x,(int,float)) and math.isfinite(x)

def _close(a,b):
    return abs(a-b)<=RTOL*max(abs(a),abs(b))

def validate(node):
    if not isinstance(node,dict):raise ValueError('Wind node must be a JSON object')
    allowed={'schema','schema_version','node_id','energy_gauge','base_face','validation','provenance',*BASE,*STATE,*OUTPUT}
    missing,unknown=allowed-set(node),set(node)-allowed-{'notes'}
    if missing or unknown:raise ValueError(f'Wind node keys: missing {sorted(missing)}, unknown {sorted(unknown)}')
    if node['schema']!=SCHEMA or node['schema_version']!=VERSION:raise ValueError('Unsupported wind node schema')
    if not isinstance(node['node_id'],str) or not re.fullmatch(ID,node['node_id']):raise ValueError('Invalid wind node_id')
    if not isinstance(node.get('notes',''),str):raise ValueError('Wind node notes must be a string')
    for k in (*BASE,*STATE,*OUTPUT):
        if not _number(node[k]) or node[k]<=0:raise ValueError(f'Wind node {k} must be positive and finite')
    if not isinstance(node['energy_gauge'],str) or not node['energy_gauge']:raise ValueError('Wind node needs an explicit energy gauge')
    f=node['base_face']
    if not isinstance(f,dict) or set(f)!=set(FACE):raise ValueError(f'base_face needs exactly {list(FACE)}')
    if any(not _number(f[k]) for k in FACE if k!='energy_gauge'):raise ValueError('base_face values must be finite numbers')
    face=InterfaceState(**f);face.validate()
    if face.energy_gauge!=node['energy_gauge']:raise ValueError('base_face and node energy gauges differ')
    R,M,rho=node['R_base_cm'],node['M_total_g'],node['rho_base']
    identities={'radius':(face.radius,R),'enclosed_mass':(face.enclosed_mass,M),'density':(face.density,rho),
        'gas_pressure':(face.gas_pressure,node['P_gas_base']),'radiation_pressure_rr':(face.radiation_pressure_rr,node['P_rad_base']),
        'gas_temperature':(face.gas_temperature,node['T_base']),'potential':(face.potential,-G*M/R),
        'L_base_erg_s':(4*PI*R**2*face.radiation_flux,node['L_base_erg_s']),
        'mdot_g_s':(4*PI*R**2*rho*face.velocity,node['mdot_g_s']),
        'Teff_emergent_K':((node['L_emergent_erg_s']/(4*PI*SIGMA*node['R_phot_cm']**2))**.25,node['Teff_emergent_K'])}
    bad=[k for k,(a,b) in identities.items() if not _close(a,b)]
    if bad:raise ValueError(f'Wind node identities fail at relative {RTOL:g}: {bad}')
    if node['R_phot_cm']<R:raise ValueError('Photosphere lies inside the wind base')
    v=node['validation']
    if not isinstance(v,dict) or set(v)!=set(VALIDATION) or any(not isinstance(x,bool) for x in v.values()):
        raise ValueError(f'validation needs boolean {list(VALIDATION)}')
    p=node['provenance']
    if not isinstance(p,dict) or not isinstance(p.get('source_run'),str) or not p['source_run']:
        raise ValueError('provenance needs a source_run description')
    if any(not isinstance(p.get(k),str) or not re.fullmatch(r'[0-9a-f]{64}',p[k]) for k in PROVENANCE):
        raise ValueError(f'provenance needs sha256 hex digests {list(PROVENANCE)}')
    return node

def nodes_of(data):
    """One node, or an ordered sequence {schema: qstars_wind_node_set, nodes: [...]}."""
    if isinstance(data,dict) and data.get('schema')==SET_SCHEMA:
        if set(data)-{'schema','schema_version','nodes','notes'} or data.get('schema_version')!=VERSION:
            raise ValueError('Unsupported wind node set')
        nodes=data.get('nodes')
        if not isinstance(nodes,list) or not nodes:raise ValueError('A wind node set needs an ordered node list')
    else:nodes=[data]
    for n in nodes:validate(n)
    if len({n['node_id'] for n in nodes})!=len(nodes):raise ValueError('Duplicate wind node_id')
    if len({n['energy_gauge'] for n in nodes})!=1:raise ValueError('Wind nodes use different energy gauges')
    P=np.log([[n[k] for k in BASE] for n in nodes])
    if np.any(np.all(np.diff(P,axis=0)==0,axis=1)):raise ValueError('Consecutive wind nodes share base data')
    return tuple(nodes)

_CACHE={}
def load(path):
    p=Path(path).resolve();st=p.stat();key=(str(p),st.st_mtime_ns,st.st_size)
    if key not in _CACHE:_CACHE[key]=nodes_of(json.loads(p.read_text()))
    return _CACHE[key]

def interface_check(nodes,eos,eos_name):
    """Thin-face fluxes: node base face against the evolution EOS at the same state.

    Same geometry, potential, mdot and lab flux; evolution gas P/u and LTE
    radiation E=aT^4, P=E/3. compare_faces refuses a gauge mismatch, which is
    a rejection. A flux defect is recorded only; it does not change the ledger.
    """
    gauge=EVOLUTION_GAUGES.get(eos_name,'')
    out=[]
    for n in nodes:
        b=InterfaceState(**n['base_face']);T=b.gas_temperature;E=A*T**4
        q=eos.evaluate(b.density,T)
        e=InterfaceState(radius=b.radius,enclosed_mass=b.enclosed_mass,density=b.density,velocity=b.velocity,
            gas_pressure=float(q['pg']),gas_internal_energy=float(q['u'])-E/b.density,radiation_energy=E,
            radiation_flux=b.radiation_flux,radiation_pressure_rr=E/3,potential=b.potential,energy_gauge=gauge,gas_temperature=T)
        try:r=compare_faces(e,b)
        except ValueError as exc:
            raise DomainError(f"Wind node {n['node_id']} base face rejected (evolution EOS gauge {gauge or 'undeclared'!r}, "
                              f"node gauge {b.energy_gauge!r}): {exc}") from exc
        out.append(dict(node_id=n['node_id'],evolution_gauge=gauge,node_gauge=b.energy_gauge,
                        conservation_passed=r['conservation_passed'],
                        relative_defects={k:v['relative'] for k,v in r['defects'].items()},
                        node_flags_all_true=all(n['validation'].values())))  # claims, not verified evidence
    return out

def json_schema():
    """Draft 2020-12 form of the field contract; identities are enforced by validate()."""
    num={'type':'number','exclusiveMinimum':0}
    sha={'type':'string','pattern':'^[0-9a-f]{64}$'}
    node={'type':'object','additionalProperties':False,
          'required':['schema','schema_version','node_id','energy_gauge',*BASE,*STATE,*OUTPUT,'base_face','validation','provenance'],
          'properties':{'schema':{'const':SCHEMA},'schema_version':{'const':VERSION},'node_id':{'type':'string','pattern':f'^{ID}$'},
              'energy_gauge':{'type':'string','minLength':1},'notes':{'type':'string'},
              **{k:dict(num,description=d) for k,d in {**BASE,**STATE,**OUTPUT}.items()},
              'base_face':{'type':'object','additionalProperties':False,'required':list(FACE),
                           'properties':{k:{'type':'string','minLength':1} if k=='energy_gauge' else {'type':'number'} for k in FACE}},
              'validation':{'type':'object','additionalProperties':False,'required':list(VALIDATION),
                            'description':'Producer claims; a proposal pending Gate 4 and user approval. Not verified evidence.',
                            'properties':{k:{'type':'boolean'} for k in VALIDATION}},
              'provenance':{'type':'object','required':['source_run',*PROVENANCE],
                            'properties':{'source_run':{'type':'string','minLength':1},**{k:sha for k in PROVENANCE}}}}}
    return {'$schema':'https://json-schema.org/draft/2020-12/schema','$id':f'{SCHEMA}_v{VERSION}','$defs':{'node':node},
            'oneOf':[{'$ref':'#/$defs/node'},
                     {'type':'object','additionalProperties':False,'required':['schema','schema_version','nodes'],
                      'properties':{'schema':{'const':SET_SCHEMA},'schema_version':{'const':VERSION},'notes':{'type':'string'},
                                    'nodes':{'type':'array','minItems':1,'items':{'$ref':'#/$defs/node'}}}}],
            'description':f'Cross-field identities (relative {RTOL:g}) are enforced by qstars.wind_node.validate.'}
