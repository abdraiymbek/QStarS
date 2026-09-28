exec(open('/tmp/qstars_wind_probe.py').read().split('rows=[]')[0].replace('spherical_wind_probe','report8_base_matching_probe'))
micro=FreeEOSOpacity(eos,mesa);rows=[]
for teff,g in [(6000,1),(8000,3)]:
 f=np.load(f'results/20260928T044002.683023Z_atmosphere_reference/T{teff}_g{float(g)}_baseline.npz')
 rp=np.sqrt(G*1e4*MSUN/g);ip=np.argmin(abs(f['tau']-2/3));rb=rp-(f['inward_depth'][-1]-f['inward_depth'][ip]);rhob=f['density'][-1];Tb=f['temperature'][-1];L=4*np.pi*rp**2*mesa.radiation_constant*C/4*teff**4
 inputs=dict(mass=1e4*MSUN,radius=rb,density=rhob,temperature=Tb,luminosity=L)
 print(teff,inputs,flush=True)
 for guess in [[rb*1.001,rhob,10000],[rb*1.003,rhob*.1,12000]]:
  t=time.perf_counter()
  try:
   w,s=match_base(micro,**inputs,sonic_guess=guess,rtol=1e-8,max_nfev=60)
   ff,ss=w.integrate(inner_ratio=rb/w.rs,outer_ratio=2,endpoint='diffusion_limit',rtol=1e-8)
   s.update(wind=ss);np.savez_compressed(out/f'T{teff}_guess{len(rows)}.npz',**ff)
  except Exception as e:s={'status':'numerical_or_domain_failure','reason':str(e),'traceback':traceback.format_exc()}
  s.update(inputs=inputs,guess=guess,seconds=time.perf_counter()-t);rows.append(s);io.json_write(out/'cases.json',rows);print(teff,{k:v for k,v in s.items() if k in ['max_log_residual','reason','wind','seconds']},flush=True)
io.json_write(out/'provenance_check.json',{'source_snapshot_verified':io.finish_provenance(out)})
