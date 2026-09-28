#!/usr/bin/env python3
"""Figure 2 vector samples; fixed PDF/tick calibration, never a model fit."""
from pathlib import Path
import csv,hashlib,json
import numpy as np
from pypdf import PdfReader
from pypdf.generic import ContentStream
source=Path('sources/cache/nakauchi2017.pdf')
expected='05cba5d4373c94269d13a37245aad00996e05c97ef3711ad56e6b0e19682537a'
if hashlib.sha256(source.read_bytes()).hexdigest()!=expected:raise ValueError('Recalibrate changed PDF')
page=PdfReader(source).pages[3];paths=[];current=[]
for val,op in ContentStream(page.get_contents(),page.pdf).operations:
 if op==b'm':current=[list(map(float,val))]
 elif op==b'l':current.append(list(map(float,val)))
 elif op==b'S':
  if len(current)>10:paths.append(np.array(current))
  current=[]
if len(paths)!=473:raise ValueError('Unexpected path inventory')
# Visually checked Figure 2 on page 4. Units are PDF stream coordinates.
# x: log10(r/Rsun)=1,2,3 at 3699.32,4263.76,4828.19.
# y: T=1e4,1e5,1e6 at 5071.69,5318.63,5565.57;
# v=1e6,1e7 at 6617.51,6987.92; log rho=-12,-10 at 5738.52,5950.19.
blocks={'temperature':[(0,34),(34,69),(69,96),(96,120)],
        'velocity':[(120,154),(189,229),(263,291),(318,343)],
        'gas_sound_speed':[(154,189),(229,263),(291,318),(343,367)],
        'density':[(367,397),(397,431),(431,454),(454,473)]}
calibration={'temperature':(5071.69,246.94,4.),'velocity':(6617.51,370.41,6.),
             'gas_sound_speed':(6617.51,370.41,6.),'density':(5738.52,105.835,-12.)}
rows=[]
for key,ranges in blocks.items():
 for (start,end),density in zip(ranges,[9.71e-9,9.65e-9,9.85e-9,10.5e-9]):
  p=np.concatenate(paths[start:end]);p=p[np.r_[0,np.arange(20,len(p)-1,20),len(p)-1]]
  y0,scale,log0=calibration[key]
  for x,y in p:
   rows.append(dict(field=key,sonic_density=density,radius_rsun=10**(1+(x-3699.32)/564.435),
                    value=10**(log0+(y-y0)/scale),log_y_uncertainty=(9.79918/2+1)/scale*np.log(10),
                    log_r_uncertainty=(9.79918/2+1)/564.435*np.log(10),x_pdf=x,y_pdf=y))
p=Path('benchmarks/nakauchi2017_figure2.csv')
with p.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
p.with_suffix('.json').write_text(json.dumps(dict(source='https://arxiv.org/abs/1611.09601',sha256=expected,
 page=4,figure=2,calibration=calibration,blocks=blocks,rows=len(rows),
 uncertainty='Half line width plus one stream unit in both log coordinates; not original model data',
 no_coefficients_fitted=True),indent=2)+'\n')
print(p,len(rows))
