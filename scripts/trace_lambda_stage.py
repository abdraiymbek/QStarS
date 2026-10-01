#!/usr/bin/env python3
"""Read-only analysis of the QSTARS step-scale traces from a fixed-T lambda run.

Parses a run directory's QSTARS_STEQ_TRACE.bin (per-call records: int32 call
counter, int32 NT, int32 ND, then POPS(NT,ND) and STEQ(NT,ND) as Fortran
column-major float64) and QSTARS_SCALE_TRACE (one formatted line per traced
depth), and the OUTGEN text after the manifest-attested inherited prefix.
Writes <run>/lambda_trajectory.json and prints a per-iteration summary.

Population index map (1-based): H0 1, HI 2-31, HII 32, HeI 33-101,
He2 102-131, HeIII 132, ne 133, T 134.
"""
import argparse, json, re, struct, sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
for _p in (str(_ROOT / 'src'), str(_ROOT / 'backends/cmfgen')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from run_component import current_native_log

TRACE_DEPTHS = list(range(14, 27))                      # 1-based depths to tabulate
POP_INDEX = dict(HII=32, HeI_ground=33, HeII_ground=102, HeIII=132, ne=133)
IT_MARK = re.compile(r'Current great iteration count is\s+(\d+)')
LUM = re.compile(r'Luminosity of star \(d=1,ND\)\(iteration\s+(\d+)\) is:\s+(\S+)\s+(\S+)')
MIN_SCALE = re.compile(r'minimum value of scale for Major species is:\s+(\S+)')
MAX_CHG = re.compile(r'Maximum % (increase|decrease) at depth\s+(\d+) is\s+(\S+)\s+\(([^)]+)\)(.*)')


def read_steq_bin(path):
    recs = []
    data = path.read_bytes()
    off = 0
    while off + 12 <= len(data):
        cnt, nt, nd = struct.unpack_from('<3i', data, off)
        n = nt * nd
        need = 12 + 16 * n
        if off + need > len(data):
            break  # partial trailing record (e.g. run in progress)
        pops = np.frombuffer(data, '<f8', n, off + 12).reshape(nt, nd, order='F')
        steq = np.frombuffer(data, '<f8', n, off + 12 + 8 * n).reshape(nt, nd, order='F')
        recs.append(dict(call=cnt, nt=nt, nd=nd, pops=pops, steq=steq))
        off += need
    return recs, off


def last_segment(items, key):
    """Drop records inherited from a continuation source: the SAVE-d call
    counter restarts at 1 in a new process, so a decrease marks the seam."""
    cut = 0
    for i in range(1, len(items)):
        if key(items[i]) < key(items[i - 1]):
            cut = i
    return items[cut:]


def read_scale_trace(path):
    rows = []
    for line in path.read_text().splitlines():
        f = line.split()
        if len(f) < 12:
            continue
        try:
            rows.append(dict(call=int(f[0]), lambda_it=f[1] == 'T', depth=int(f[2]),
                             scale=float(f[3]), t1=float(f[4]), jdec=int(f[5]),
                             t2=float(f[6]), jinc=int(f[7]), pops_jinc=float(f[8]),
                             pops_jdec=float(f[9]), ne=float(f[10]), steq_t=float(f[11])))
        except ValueError:
            continue
    return rows


def iteration_blocks(text):
    marks = list(IT_MARK.finditer(text))
    return [(int(m.group(1)), text[m.start():marks[k + 1].start() if k + 1 < len(marks) else len(text)])
            for k, m in enumerate(marks)]


def analyze(directory):
    directory = Path(directory)
    recs, used = read_steq_bin(directory / 'QSTARS_STEQ_TRACE.bin')
    recs = last_segment(recs, lambda r: r['call'])
    bin_size = (directory / 'QSTARS_STEQ_TRACE.bin').stat().st_size
    scales = last_segment(read_scale_trace(directory / 'QSTARS_SCALE_TRACE'),
                          lambda r: r['call'])
    by_call = {}
    for r in scales:
        by_call.setdefault(r['call'], []).append(r)
    blocks = iteration_blocks(current_native_log(directory))
    iterations = []
    for i, rec in enumerate(recs):
        pops, steq, nt = rec['pops'], rec['steq'], rec['nt']
        amag = np.abs(steq)
        j, d = np.unravel_index(np.argmax(amag), amag.shape)
        ne = pops[nt - 2, :]                                    # row 133 = electron density
        major = pops >= 1.0e-10 * ne[None, :]                   # FIDDLE major-species cut
        amaj = np.where(major, amag, np.nan)
        jm, dm = np.unravel_index(np.nanargmax(amaj), amaj.shape)
        dt_row = np.abs(steq[nt - 1, :])                        # row 134 = T correction
        dt = int(np.argmax(dt_row))
        rows = by_call.get(rec['call'], [])
        smin = min(rows, key=lambda r: r['scale']) if rows else None
        oit, blk = blocks[i] if i < len(blocks) else (None, '')
        lm = LUM.search(blk)
        sm = MIN_SCALE.search(blk)
        chg = [dict(kind=m.group(1), depth=int(m.group(2)), value=float(m.group(3)),
                    mode=m.group(4).strip(), tail=m.group(5).strip())
               for m in MAX_CHG.finditer(blk)]
        depth_tab = {}
        for dp in TRACE_DEPTHS:
            if dp <= rec['nd']:
                depth_tab[str(dp)] = {name: float(pops[idx - 1, dp - 1])
                                      for name, idx in POP_INDEX.items()}
        iterations.append(dict(
            call=rec['call'], outgen_iteration=oit,
            max_abs_steq_all=dict(value=float(amag[j, d]), index=int(j) + 1, depth=int(d) + 1),
            max_abs_steq_major=dict(value=float(amaj[jm, dm]), index=int(jm) + 1, depth=int(dm) + 1),
            max_abs_steq_t=dict(value=float(dt_row[dt]), depth=dt + 1),
            max_changes=chg,
            min_scale=(dict(value=smin['scale'], depth=smin['depth'],
                            jinc=smin['jinc'], jdec=smin['jdec']) if smin else None),
            outgen_min_scale=(float(sm.group(1)) if sm else None),
            luminosity=([float(lm.group(2)), float(lm.group(3))] if lm else None),
            poor_flux_accuracy='poor flux accuracy' in blk,
            populations=depth_tab))
    result = dict(directory=str(directory), trace_records=len(recs),
                  trace_bin_bytes=bin_size, trace_bin_bytes_parsed=used,
                  expected_record_bytes=(12 + 16 * recs[0]['nt'] * recs[0]['nd'] if recs else None),
                  scale_trace_rows=len(scales), outgen_iterations=[b[0] for b in blocks],
                  iterations=iterations)
    out = directory / 'lambda_trajectory.json'
    out.write_text(json.dumps(result, indent=1) + '\n')
    hdr = (f"{'it':>3} {'call':>4} {'max|STEQ| all':>14} {'@idx':>4} {'@d':>3} {'max|STEQ| maj':>14} {'@idx':>4} {'@d':>3} "
           f"{'max|STEQ(134)|':>14} {'@d':>3} {'min SCALE':>10} {'@d':>3} {'L(1)':>11} {'L(ND)':>11} flux  maxchg")
    print(hdr)
    for it in iterations:
        a, m, s, t = it['max_abs_steq_all'], it['max_abs_steq_major'], it['min_scale'], it['max_abs_steq_t']
        L = it['luminosity'] or [float('nan'), float('nan')]
        ch = ';'.join(f"{c['kind']} {c['value']:.2E}@d{c['depth']}({c['mode']})" for c in it['max_changes'])
        print(f"{it['outgen_iteration'] or -1:>3} {it['call']:>4} {a['value']:>14.6E} {a['index']:>4} {a['depth']:>3} "
              f"{m['value']:>14.6E} {m['index']:>4} {m['depth']:>3} {t['value']:>14.6E} {t['depth']:>3} "
              f"{(s['value'] if s else float('nan')):>10.3E} {(s['depth'] if s else -1):>3} "
              f"{L[0]:>11.5E} {L[1]:>11.5E} {'POOR' if it['poor_flux_accuracy'] else 'ok'}  {ch}")
    print(f"wrote {out}")
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('directory', type=Path)
    a = p.parse_args()
    analyze(a.directory)
