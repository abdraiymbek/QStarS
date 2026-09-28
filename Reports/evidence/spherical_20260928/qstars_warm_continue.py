exec(open('/tmp/qstars_wind_probe.py').read().split('rows=[]')[0].replace('spherical_wind_probe','warm_matching_continuation'))
import json
x=json.load(open('results/20260928T172021.293589Z_spherical_wind_validation/report8_base_attempts.json'))[1]
b=min(x['trials'],key=lambda t:sum(v*v for v in t['residual']));guess=[b['sonic_radius'],b['sonic_density'],b['sonic_temperature']]
t=time.perf_counter()
try:
 w,s=match_base(FreeEOSOpacity(eos,mesa),**x['inputs'],sonic_guess=guess,rtol=1e-9,max_nfev=400)
 f,rs=w.integrate(inner_ratio=x['inputs']['radius']/w.rs,outer_ratio=4,rtol=1e-9,sonic_offset=1e-7)
 s.update(status='matched_local_component',wind=rs);np.savez_compressed(out/'wind.npz',**f)
except Exception as e:s=dict(status='matching_failure',reason=str(e),trials=getattr(e,'trials',[]),traceback=traceback.format_exc())
s.update(seconds=time.perf_counter()-t,guess=guess,inputs=x['inputs']);io.json_write(out/'match.json',s)
io.json_write(out/'provenance_check.json',{'source_snapshot_verified':io.finish_provenance(out)})
print({k:v for k,v in s.items() if k not in ['trials','traceback']},flush=True)
