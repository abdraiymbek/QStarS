"""Isolated-process performance measurement and saved-output figures."""
import json,subprocess,sys,time
from pathlib import Path

def main():
    reports=[]
    for file in ['closed_relaxation','hot_accreting','two_stream_pulse','bra_luminosity','dotan_launch','pressure_confined']:
        start=time.perf_counter()
        proc=subprocess.run([sys.executable,'-m','qstars','run',f'configs/{file}.yaml'],capture_output=True,text=True)
        stdout=proc.stdout.strip().splitlines()
        if not stdout:raise RuntimeError(proc.stderr)
        path=Path(stdout[0]);report=json.loads((path/'termination.json').read_text())
        report.update(config=f'configs/{file}.yaml',path=str(path),process_wall_seconds=time.perf_counter()-start,exit_code=proc.returncode)
        reports.append(report)
        print(file,report['status'],f"{report['runtime_seconds']:.3f}s",flush=True)
        if (path/'history.csv').exists():
            subprocess.run([sys.executable,'-m','qstars','plot',str(path)],check=True)
    out=Path('Reports/evidence');out.mkdir(parents=True,exist_ok=True)
    stamp=time.strftime('%Y%m%dT%H%M%S')
    (out/f'performance_{stamp}.json').write_text(json.dumps(reports,indent=2)+'\n')
    return reports

if __name__=='__main__':main()
