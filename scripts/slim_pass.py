#!/usr/bin/env python3
"""Provenance-safe disk compaction of completed CMFGEN pass directories.

A pass directory starts as a full copy of its source (continuation.json "source"; copytree
preserves size and mtime). After the run, this removes
  (a) files identical by name, size and mtime to the source copy (restorable from the source);
  (b) bulky native outputs that no QStarS step reads from a pass: JH_AT_CURRENT_TIME,
      STEQ_VALS, J_COMP, EDDFACTOR, GreyDiagnostics, MOM_J_ERRORS, QSTARS_STEQ_TRACE.bin.
Never removed: QSTARS_RESID_TRACE.bin, QSTARS_BA_TRACE.bin, QSTARS_RJ.bin (unless --drop-rj),
QSTARS_JV_RAD.bin, QSTARS_DIRS.bin, QSTARS_STEP.bin, pairs.npz, run_manifest.json, VADAT, OUTGEN,
logs and every file changed by the run. Each removal is listed with size (and, for outputs,
sha256) in QSTARS_SLIM.json, so the directory remains an auditable record.
"""
import argparse, hashlib, json, os
from pathlib import Path

BULKY = ('JH_AT_CURRENT_TIME', 'STEQ_VALS', 'J_COMP', 'EDDFACTOR', 'GreyDiagnostics', 'MOM_J_ERRORS',
         'QSTARS_STEQ_TRACE.bin')
KEEP = {'QSTARS_RESID_TRACE.bin', 'QSTARS_BA_TRACE.bin', 'QSTARS_RJ.bin', 'QSTARS_JV_RAD.bin', 'QSTARS_DIRS.bin',
        'QSTARS_STEP.bin', 'pairs.npz', 'run_manifest.json', 'run_result.json', 'VADAT', 'OUTGEN', 'native.log',
        'continuation.json', 'closure_control.json', 'QSTARS_SLIM.json', 'IN_ITS', 'IT_SPECIFIER', 'MODEL'}


def sha256(p, chunk=1 << 22):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(chunk), b''):
            h.update(b)
    return h.hexdigest()


def slim(d, drop_rj=False, dry=False):
    d = Path(d)
    if (d / 'QSTARS_SLIM.json').exists() or not (d / 'run_manifest.json').exists():
        return 0
    cont = d / 'continuation.json'
    src = Path(json.loads(cont.read_text())['source']) if cont.exists() else None
    removed, freed = [], 0
    for f in sorted(d.iterdir()):
        if not f.is_file() or f.is_symlink():
            continue
        name = f.name
        if name in KEEP and not (drop_rj and name == 'QSTARS_RJ.bin'):
            continue
        st = f.stat()
        why = None
        if name in BULKY or (drop_rj and name == 'QSTARS_RJ.bin'):
            why = 'bulky output'
        elif src is not None and (src / name).is_file():
            s2 = (src / name).stat()
            if s2.st_size == st.st_size and int(s2.st_mtime) == int(st.st_mtime):
                why = 'identical to source'
        if why is None:
            continue
        rec = dict(name=name, bytes=st.st_size, reason=why)
        if why == 'bulky output':
            rec['sha256'] = sha256(f)
        removed.append(rec)
        freed += st.st_size
        if not dry:
            f.unlink()
    if not dry:
        (d / 'QSTARS_SLIM.json').write_text(json.dumps(dict(source=str(src) if src else None, removed=removed,
                                                             freed_bytes=freed), indent=0) + '\n')
    return freed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('roots', nargs='+', type=Path)
    ap.add_argument('--drop-rj', action='store_true', help='also remove stored J (only for passes never replayed again)')
    ap.add_argument('--dry', action='store_true')
    a = ap.parse_args()
    total = 0
    for r in a.roots:
        for dirpath, dirnames, filenames in os.walk(r):
            if 'run_manifest.json' in filenames:
                total += slim(dirpath, a.drop_rj, a.dry)
    print(f'freed {total / 1e9:.2f} GB' + (' (dry run)' if a.dry else ''))


if __name__ == '__main__':
    main()
