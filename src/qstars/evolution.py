"""Transactional evolution: rejected candidates never mutate accepted state."""
import json,time,resource,copy,hashlib
from .constants import YEAR
from .eos import Composition,IdealGasRadiation
from .opacity import ElectronScattering,OpacityTable
from .structure import initialize
from . import ledger
from .engine import rates
from .winds import rate as wind_rate
from .mass_transport import exchange,DomainError
from .transport import luminosities
from .numerics import thermal_mechanical_step,StepRejected
from .validity import check,assumptions
from .controls import next_event
from . import io

def physics(config):
    eos=IdealGasRadiation(Composition(**config["composition"]))
    op=ElectronScattering(eos.composition) if config["opacity"]["model"]=="electron_scattering" else OpacityTable(config["opacity"]["table_path"])
    return eos,op

def attempt_step(old,dt,eos,opacity,config):
    engine=rates(old,eos,config["engine"])
    L,_=luminosities(old,eos,opacity,config,engine["luminosity"])
    wind=wind_rate(old,eos,opacity,config["wind"],max(0,float(L[-1])))
    transported,stream=exchange(old,dt,eos,config,engine,wind)
    state,diag=thermal_mechanical_step(transported,dt,eos,opacity,config,engine["luminosity"])
    if stream["wind_cost"]>max(0.,dt*diag["luminosity"])*(1+1e-10):
        raise DomainError("Outgoing luminosity cannot provide the requested escaping wind energy")
    state.ledger["radiated_energy"]-=stream["wind_cost"]
    state.ledger["wind_energy"]+=stream["wind_cost"]
    state.ledger["wind_mechanical_energy"]+=stream["wind_cost"]
    errors,validity,conservation=check(state,eos,opacity,config)
    if errors: raise DomainError("; ".join(errors))
    return state,{**diag,**stream,**validity,**conservation}

def run(config,output_parent="results",restart=None,stop_after_steps=None):
    config=copy.deepcopy(config)
    path=io.create_directory(output_parent,config["name"])
    if config['opacity']['model']=='table':
        from pathlib import Path
        source=Path(config['opacity']['table_path'])
        data=source.read_bytes()
        table=(path/'opacity_table.npz').resolve();table.write_bytes(data)
        config['opacity']['table_path']=str(table)
        io.json_write(path/'external_inputs.json',{'opacity_source_path':str(source),
                      'opacity_sha256':hashlib.sha256(data).hexdigest()})
    io.provenance(path,config)
    io.json_write(path/"validity_warnings.json",assumptions(config))
    start=time.perf_counter(); accepted=0; retries=0; state=None
    status="numerical_failure"; reason="initialization did not complete"
    def log(event,**data):
        with (path/"events.jsonl").open("a") as f:
            f.write(json.dumps(dict(event=event,**data),default=float,allow_nan=False)+"\n")
    try:
        eos,opacity=physics(config)
        if restart is not None: state=restart.copy()
        else:
            state,error=initialize(config,eos); ledger.initialize(state,eos)
            log("initialized",hydrostatic_residual=error)
        errors,_,_=check(state,eos,opacity,config)
        if errors: raise DomainError("; ".join(errors))
        io.append_history(path/"history.csv",io.history_row(state,eos,opacity,config))
        io.save_profile(path/f"profile_{state.steps:07d}.npz",state,eos,opacity,config)
        end=config["evolution"]["end_time_years"]*YEAR
        c=config["numerics"]
        while state.time<end*(1-1e-14):
            if config['evolution']['max_wall_seconds'] and time.perf_counter()-start>=config['evolution']['max_wall_seconds']:
                status,reason='control_stop','configured wall time budget reached';break
            if state.steps>=config["evolution"]["max_steps"]:
                status,reason="control_stop","maximum accepted steps reached"; break
            if stop_after_steps is not None and accepted>=stop_after_steps:
                status,reason="control_stop","requested checkpoint stop"; break
            event=min(next_event(config[s][k],state.time/YEAR) for s,k in
                [("feeding","rate_msun_per_year"),("feeding","temperature_K"),("engine","coupling"),
                 ("engine","rest_mass_rate_msun_per_year"),("wind","rate_msun_per_year"),("boundary","luminosity_erg_s")])
            dt=min(state.next_dt,c["max_dt_years"]*YEAR,end-state.time,event*YEAR-state.time)
            last_exception=None
            for retry in range(c["max_retries"]+1):
                try:
                    candidate,diag=attempt_step(state,dt,eos,opacity,config)
                    break
                except (StepRejected,DomainError,ValueError,FloatingPointError) as exc:
                    retries+=1; last_exception=exc
                    log("rejected",time_seconds=state.time,dt_seconds=dt,error=str(exc))
                    dt*=.5
                    if dt<c["min_dt_years"]*YEAR: raise exc
            else: raise last_exception
            state=candidate; accepted+=1
            state.next_dt=min(dt*1.25,c["max_dt_years"]*YEAR)
            if c["refine_every_steps"] and state.steps % c["refine_every_steps"]==0 and len(state.dm)*2<=c["max_zones"]:
                import numpy as np
                from .mesh import remap
                before=state.energies(eos)["total"]
                remapped=remap(state,np.repeat(state.dm/2,2),eos)
                errors,_,_=check(remapped,eos,opacity,config)
                if errors: log("remap_rejected",errors=errors)
                else:
                    state=remapped
                    log("remapped",zones=len(state.dm),energy_change=state.energies(eos)["total"]-before)
            log("accepted",step=state.steps,dt_seconds=dt,**diag)
            io.append_history(path/"history.csv",io.history_row(state,eos,opacity,config))
            if state.steps%config["output"]["profile_every_steps"]==0:
                io.save_profile(path/f"profile_{state.steps:07d}.npz",state,eos,opacity,config)
            if state.steps%config["output"]["checkpoint_every_steps"]==0:
                io.save_checkpoint(path/"checkpoint.npz",state,config)
        else: status,reason="completed_requested_interval","requested end time reached within configured domain"
    except DomainError as exc:
        status,reason="model_domain_exit",str(exc); log(status,error=reason)
    except KeyboardInterrupt:
        status,reason='control_stop','execution interrupted; last accepted state preserved'
        log(status,reason=reason)
    except (RuntimeError,ValueError,ArithmeticError) as exc:
        status,reason="numerical_failure",str(exc); log(status,error=reason)
    finally:
        if state is not None and state.ledger:
            io.save_checkpoint(path/"checkpoint.npz",state,config)
            try: io.save_profile(path/f"profile_{state.steps:07d}.npz",state,eos,opacity,config)
            except (ValueError,RuntimeError,ArithmeticError) as exc:
                log('final_profile_unavailable',error=str(exc))
        # macOS ru_maxrss is bytes; Linux reports KiB.
        import sys
        memory=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform=="darwin" else 1024)
        report=dict(status=status,reason=reason,accepted_steps_this_run=accepted,retries=retries,
                    runtime_seconds=time.perf_counter()-start,process_peak_rss_MiB=memory,
                    final_time_years=state.time/YEAR if state is not None else 0.,v1_validated=False)
        if state is not None and state.ledger: report["conservation"]=ledger.residuals(state,eos)
        io.json_write(path/"termination.json",report)
    return path,report
