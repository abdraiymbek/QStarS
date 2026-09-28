exec(open('/tmp/qstars_wind_probe.py').read().split('rows=[]')[0].replace('spherical_wind_probe','cool_sonic_probe'))
rows=[]
micro=FreeEOSOpacity(eos,mesa)
for T in [1e4,1.2e4,1.6e4,2e4]:
 for rho in [1e-11,3e-11,1e-10]:
  name=f'T{T}_rho{rho}';t=time.perf_counter()
  try:
   w=DiffusionWind(micro,1e4*MSUN,1.152e15,rho,T)
   f,s=w.integrate(inner_ratio=.99,outer_ratio=2,endpoint='diffusion_limit',samples=301,rtol=1e-7,max_log_step=.002)
   np.savez_compressed(out/(name+'.npz'),**f)
  except Exception as e:s={'status':'numerical_or_domain_failure','reason':str(e),'traceback':traceback.format_exc()}
  s.update(name=name,seconds=time.perf_counter()-t);rows.append(s);io.json_write(out/'cases.json',rows)
  print(name,{k:s[k] for k in s if k in ['reason','terminations','sonic_slope','endpoint_temperature','endpoint_radius','shell_mass_fraction','seconds']},flush=True)
io.json_write(out/'provenance_check.json',{'source_snapshot_verified':io.finish_provenance(out)})
