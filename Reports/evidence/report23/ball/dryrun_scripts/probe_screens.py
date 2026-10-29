"""DRYRUN: enforced / report-only screens of the dry-run relaxed row 0 (no table comparison)."""
from dry import *
from qstars import validity
import sys
for d in sys.argv[1:]:
    state,config=quasistatic.load_checkpoint(d+'/DRYRUN_row0.npz')
    eos,opacity=vbt.physics(config)
    errors,diag,_=validity.check(state,eos,opacity,config,policy='published_assumption',face_luminosity=state.face_luminosity)
    print('DRYRUN',d,'enforced errors:',errors,'warnings:',diag['screen_warnings'])
    print('   ',{k:(v['value'],v['violated'],v['action']) for k,v in diag['screens'].items() if v['action'] in ('enforced','report_only')})
    print('    clamped mass Msun',diag['opacity_clamped_mass_g']/MSUN,'photosphere clamped',diag['opacity_photosphere_clamped'],'Gamma_ph',diag['photosphere_eddington_factor'])
    s=quasistatic.structure_summary(state,eos,opacity,config)
    print('    structure',{k:v for k,v in s.items() if k in ('max_density_inversion','stability_ratio','convective_mass_fraction')})
