#!/usr/bin/env python3
"""Extract the plotted vector paths of FR16 Figure 4, not source model data.

Axis calibration and path IDs are explicit and checked against this exact PDF.
Uncertainty includes half the visible stroke width and the published base-Mach
range (the latter is propagated separately by validate_wind.py). No curve is
used to fit any physical coefficient.
"""
import csv,hashlib,json
from pathlib import Path
import numpy as np
from pypdf import PdfReader
from pypdf.generic import ContentStream

source=Path('sources/cache/fiacconi2016.pdf')
page=PdfReader(source).pages[6]
paths=[];current=[];dash=[];width=0.
for values,op in ContentStream(page.get_contents(),page.pdf).operations:
    if op==b'd':dash=values
    elif op==b'w':width=float(values[0])
    elif op==b'm':current=[list(map(float,values))]
    elif op==b'l':current.append(list(map(float,values)))
    elif op==b'S':
        if len(current)>10:paths.append(dict(points=np.asarray(current),width=width,dash=str(dash)))
        current=[]
if len(paths)!=22:raise ValueError('PDF vector paths changed; calibrate again before extracting')
# x at r/R=1, y at value=1, y units per decade, visible y limits.
axes=[(935.441,6655.65,693.14,6139.45,7348.79),
      (2447.11,6819.72,229.925,6139.45,7348.79),
      (3958.77,6774.91,635.46,6139.45,7348.79),
      (935.441,5679.88,441.305,4930.12,6139.45),
      (2447.11,5461.27,171.51,4930.12,6139.45),
      (3958.77,5763.25,416.565,4930.12,6139.45)]
groups=[[(0,'velocity_ratio'),(1,'sound_speed_ratio'),(2,'mach')],
        [(3,'optical_depth_proxy'),(4,'u'),(5,'density_ratio')],
        [(6,'kinetic_luminosity'),(7,'advected_luminosity'),(8,'luminosity')],
        [(9,'velocity_ratio'),(11,'sound_speed_ratio'),(13,'mach')],
        [(15,'optical_depth_proxy'),(16,'u'),(17,'density_ratio')],
        [(19,'kinetic_luminosity'),(20,'advected_luminosity'),(21,'luminosity')]]
rows=[]
for panel,group in enumerate(groups):
    x0,y0,yscale,ymin,ymax=axes[panel]
    for index,field in group:
        curve=paths[index];half=curve['width']/2
        for x,y in curve['points']:
            log_r=(x-x0)/291.06125
            if not -1e-6<=log_r<=4.00001 or not ymin+half<y<ymax-half:continue
            rows.append(dict(alpha=.4 if panel<3 else 7.64,field=field,
              radius_ratio=10**max(0,log_r),value=10**((y-y0)/yscale),
              log_uncertainty=(half+1)/yscale*np.log(10),path_id=index,x_pdf=x,y_pdf=y))
out=Path('benchmarks/fiacconi2016_figure4.csv')
with out.open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
meta=dict(source_url='https://arxiv.org/pdf/1509.07511',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
 page_one_based=7,figure=4,method='Vector path vertex digitization with explicit logarithmic tick calibration',
 extraction_scope='Plotted polyline vertices inside visible axes; not original simulation samples',
 axes=axes,x_units_per_decade=291.06125,rows=len(rows),
 uncertainty='Half stroke width plus one PDF coordinate unit; published base Mach 0.8 +/- 0.05 propagated separately',
 parameters=dict(energy=5.,beta=100.,base_mach=.8,base_mach_range=[.75,.85]),
 no_physical_parameters_fitted=True)
out.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
print(out,len(rows))
