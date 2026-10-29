#!/usr/bin/env python3
"""Keep a long cool_newton_krylov.py run inside the disk budget (device transfer, 2026-10-29).

    retain_runs.py OUT [--keep-iterates 3] [--idle-min 30] [--offload DEST] [--slim] [--budget-gb 100] [--dry-run]

Pass directories of OUT (a driver --out directory) that the driver no longer reads are packed with
archive_results.py (full archive byte-exact and verified, record, manifest; the raw directory is then removed).
Never touched: the last --keep-iterates iterates (itNN_*, cols_itNN_* by their iterate number), the bases of those
history.json entries and the next base (accepted trial of the last entry) with their cols_/_from_/_tp/_ali
directories, it00_base and start_base, --protect NAMEs, and every directory with a file modified in the last
--idle-min minutes. The driver reads newton_step.npy of the last --augment (2) iterates and the same-base
it*_from_<base> passes, both inside that set.

--offload DEST   move each new full archive to DEST (local path, e.g. an external SSD, or host:path for rsync over
                 ssh) and keep only its record and manifest here; rsync checks every transfer.
--slim           run slim_pass.py first (lossy: removes bulky native outputs no QStarS step reads, with their
                 sha256 in QSTARS_SLIM.json, and files identical to the pass's source copy); for budget pressure only.
Exit status 3 (and OUT/RETENTION_OVER_BUDGET) when the project (repo + ~/qsrd-runtime) exceeds --budget-gb after
retention: start no new iterate until space is made. Every action is appended to OUT/RETENTION.log."""
import argparse, json, re, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import archive_results as ar  # noqa: E402
import slim_pass  # noqa: E402

RUNTIME = Path.home() / 'qsrd-runtime'


def iterate_no(name):
    m = re.match(r'(?:cols_)?it(\d+)', name)
    return int(m.group(1)) if m else None


def protected(out, keep, extra):
    hist = json.loads((out / 'history.json').read_text()) if (out / 'history.json').exists() else []
    # the first base is the reference of scripts/fixed_metric.py (Report 24 measures progress against it)
    bases = {Path(h['base']).name for h in hist[:1] + hist[-keep:]}
    if hist and hist[-1].get('accepted') is not None:
        bases.add(f"it{int(hist[-1]['it']) + 1:02d}_trial_{hist[-1]['accepted']}")
    bases |= {'it00_base', 'start_base'} | set(extra)
    nums = [n for n in (iterate_no(d.name) for d in out.iterdir() if d.is_dir()) if n is not None]
    newest = max(nums, default=0)

    def keep_dir(name):
        n = iterate_no(name)
        if n is not None and n > newest - keep:
            return True
        return any(name == b or name == f'cols_{b}' or name.endswith(f'_from_{b}') or name in (f'{b}_tp', f'{b}_ali')
                   for b in bases)
    return keep_dir


def idle(d, minutes):
    cut = time.time() - 60 * minutes
    return all(p.lstat().st_mtime < cut for p in d.rglob('*'))


def project_gb():
    sizes = [int(subprocess.run(['du', '-sk', str(p)], capture_output=True, text=True).stdout.split()[0])
             for p in (ROOT, RUNTIME) if p.exists()]
    return sum(sizes) * 1024 / 1e9


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('out', type=Path)
    ap.add_argument('--keep-iterates', type=int, default=3)
    ap.add_argument('--idle-min', type=float, default=30)
    ap.add_argument('--protect', nargs='*', default=[])
    ap.add_argument('--offload')
    ap.add_argument('--slim', action='store_true')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--budget-gb', type=float, default=100)
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()
    out = a.out.resolve()
    log = None if a.dry_run else open(out / 'RETENTION.log', 'a')

    def say(msg):
        line = f"{time.strftime('%Y-%m-%d %H:%M')}: {msg}"
        print(line, flush=True)
        if log:
            log.write(line + '\n'); log.flush()
    keep_dir = protected(out, a.keep_iterates, a.protect)
    ix = ar.load_index()
    for d in sorted(p for p in out.iterdir() if p.is_dir() and not p.is_symlink()):
        unit = str(d.relative_to(ar.RESULTS))
        if keep_dir(d.name) or unit in ix['units'] or not idle(d, a.idle_min):
            continue
        if a.dry_run:
            say(f'would pack {unit}')
            continue
        if a.slim:
            freed = slim_pass.slim(d)
            if freed:
                say(f'slimmed {unit}: {freed / 1e9:.2f} GB (QSTARS_SLIM.json)')
        e = ar.pack(unit, a.threads, True, None)
        if a.offload:
            full = ar.paths(unit)[0]
            subprocess.run(['rsync', '-a', '--remove-source-files', str(full), a.offload.rstrip('/') + '/'], check=True)
            e['full_location'] = a.offload
        ar.save_unit(unit, e)
        say(f"packed {unit}: {e['raw_bytes'] / 1e9:.2f} GB -> full {e['full_bytes'] / 1e9:.2f} GB"
            f"{' (offloaded to ' + a.offload + ')' if a.offload else ''}, record {e['record_bytes'] / 1e6:.1f} MB")
    gb = project_gb()
    flag = out / 'RETENTION_OVER_BUDGET'
    if gb > a.budget_gb:
        say(f'OVER BUDGET: project {gb:.1f} GB > {a.budget_gb:.0f} GB')
        if not a.dry_run:
            flag.write_text(f'{gb:.1f} GB > {a.budget_gb:.0f} GB at {time.ctime()}\n')
        sys.exit(3)
    flag.unlink(missing_ok=True)
    say(f'project {gb:.1f} GB of {a.budget_gb:.0f} GB')


if __name__ == '__main__':
    main()
