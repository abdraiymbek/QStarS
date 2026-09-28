"""Fetch primary-source PDFs for local audit; do not redistribute their text."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib, json, requests
from pypdf import PdfReader
PAPERS = {'bra2008':'0711.4078','ball2011':'1102.5098','dotan2011':'1107.3562',
'ball_thesis':'1207.5972','fiacconi2016':'1509.07511','coughlin2024':'2405.00084',
'quest2025':'2507.12651','santarelli2025':'2510.11772','hassan2025':'2510.18301',
'cantiello2025':'2512.17997','roman2026':'2603.21714','begelman2026':'2507.09085',
'atmosphere2026':'2605.12141'}
def fetch(item):
    key, aid = item
    p = Path('sources/cache') / f'{key}.pdf'
    try:
        if not p.exists():
            r = requests.get(f'https://arxiv.org/pdf/{aid}',timeout=60); r.raise_for_status(); p.write_bytes(r.content)
        reader = PdfReader(p)
        p.with_suffix('.txt').write_text('\n'.join(f'\nPAGE {i+1}\n'+page.extract_text() for i,page in enumerate(reader.pages)))
        return dict(id=key, url=f'https://arxiv.org/abs/{aid}', sha256=hashlib.sha256(p.read_bytes()).hexdigest(),pages=len(reader.pages),status='downloaded')
    except Exception as e: return dict(id=key,status='failed',error=str(e))
if __name__ == '__main__':
    out=list(ThreadPoolExecutor(4).map(fetch,PAPERS.items()))
    Path('sources/acquisition.json').write_text(json.dumps(out,indent=2)+'\n')
    for row in out: print(row)
