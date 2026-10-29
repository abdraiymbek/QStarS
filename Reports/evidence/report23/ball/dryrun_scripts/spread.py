"""DRYRUN: seed-to-seed spread of the row-0 observables (max |dln| only; values are not printed or compared with any table)."""
import json,sys,math
a,b=(json.load(open(d+'/DRYRUN_summary.json')) for d in sys.argv[1:3])
if 'row0' not in a or 'row0' not in b:print('DRYRUN spread: not both converged',a['status'],b['status']);sys.exit()
d={k:abs(math.log(a['row0'][k]/b['row0'][k])) for k in a['row0']}
print('DRYRUN spread',sys.argv[1],'vs',sys.argv[2],'zones_final',a.get('zones_final'),b.get('zones_final'),'max |dln| = %.3e'%max(d.values()),{k:'%.2e'%v for k,v in d.items()})
