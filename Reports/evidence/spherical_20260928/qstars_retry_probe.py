exec(open('/tmp/qstars_wind_probe.py').read().split('rows=[]')[0].replace('spherical_wind_probe','wind_retry_probe'))
class ExplainMicro(FreeEOSOpacity):
 def __call__(self,rho,T):
  try:return super().__call__(rho,T)
  except ValueError as e:raise ValueError(f'{e}; rho={rho}, T={T}') from e
micro=ExplainMicro(eos,mesa);w=DiffusionWind(micro,100*MSUN,50*RSUN,9.85e-9,2.5e5)
rows=[]
for h in [.02,.01,.005]:
 try:
  f,s=w.integrate(inner_ratio=.1,outer_ratio=100,endpoint='paper_photosphere',rtol=1e-9,sonic_offset=1e-7,max_log_step=h)
  np.savez_compressed(out/f'step{h}.npz',**f)
 except Exception as e:s={'reason':str(e)}
 rows.append(s);print(h,s.get('reason',s.get('terminations')),flush=True)
io.json_write(out/'results.json',rows);io.json_write(out/'provenance_check.json',{'source_snapshot_verified':io.finish_provenance(out)})
