"""Verify retained Report 16 evidence and freeze the final project interface.

Run from the repository root. This performs no native science computation.
External inputs and executables must remain available for local attestation.
"""
import gzip
import hashlib
import json
import re
import subprocess
import zipfile
from pathlib import Path

E = Path('Reports/evidence/moment_closure_20260930')
A = Path('benchmarks/results/2026-09-30-moment-closure')
cache = {}


def digest(path, compressed=False):
    path = Path(path)
    key = (path.resolve(), compressed)
    if key not in cache:
        h = hashlib.sha256()
        count = 0
        with (gzip.open(path, 'rb') if compressed else path.open('rb')) as f:
            for block in iter(lambda: f.read(8 * 1024**2), b''):
                h.update(block)
                count += len(block)
        cache[key] = (h.hexdigest(), count)
    return cache[key]


def check(path, expected, size=None, compressed=False):
    h, n = digest(path, compressed)
    if h != expected or (size is not None and n != size):
        raise RuntimeError(f'Evidence mismatch: {path}')
    return h


def read(path):
    return json.loads(Path(path).read_text())


def write(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + '\n')


manifest = read(A / 'archive_manifest.json')
for row in manifest['files']:
    check(row['original'], row['sha256'], row['bytes'])
    check(A / row['archived'], row['sha256'], row['bytes'], row['compressed'])
for row in manifest['omitted_local_files']:
    check(row['original'], row['sha256'], row['bytes'])
for row in manifest['external_atomic_links']:
    p = Path(row['original'])
    if not p.is_symlink() or str(p.resolve()) != row['target']:
        raise RuntimeError(f'External link mismatch: {p}')
    if row['sha256'] is not None:
        check(p, row['sha256'])
local_archive = dict(files_verified=len(manifest['files']),
                     omitted_local_files_verified=len(manifest['omitted_local_files']),
                     external_links_verified=len(manifest['external_atomic_links']),
                     failures=0)
write(E / 'local_archive_verification.json', local_archive)

builds = []
for number in (9, 10, 11, 12):
    root = Path(f'backends/cmfgen/build_{number}/cur_cmf')
    cap = read(E / f'build{number}.json')
    exe = root / 'exe/cmfgen_dev.exe'
    check(exe, cap['executable_sha256'])
    for name, h in cap['source_hashes'].items():
        check(root / name, h)
    if read(exe.with_suffix('.capabilities.json')) != cap:
        raise RuntimeError('Build capabilities changed')
    builds.append(dict(build=number, executable_sha256=cap['executable_sha256'],
                       recorded_sources_verified=len(cap['source_hashes'])))

for filename, before, after in [('seed_trace_manifest.json', 8, 9),
                               ('edge_trace_manifest.json', 9, 10)]:
    m = read(E / filename)
    source = 'plane/mom_jrel_v9.f'
    check(Path(f'backends/cmfgen/build_{before}/cur_cmf') / source, m['original'])
    check(Path(f'backends/cmfgen/build_{after}/cur_cmf') / source, m['traced'])
for filename, after in [('raw_nh_manifest.json', 11), ('node_moment_manifest.json', 12)]:
    m = read(E / filename)
    for name, row in m['sources'].items():
        check(Path('backends/cmfgen/build_9/cur_cmf') / name, row['original'])
        check(Path(f'backends/cmfgen/build_{after}/cur_cmf') / name, row['patched'])
check('backends/cmfgen/build_11/cur_cmf/plane/var_jrel_v5.f',
      read(E / 'raw_nh_manifest.json')['unchanged_VAR_sha256'])

campaign = read(E / 'campaign_summary.json')
for case in campaign['cases']:
    d = Path('results') / case['case']
    run = read(d / 'run_manifest.json')
    check(run['executable'], case['executable_sha256'])
    for role in ('runner', 'adapter'):
        h = case['source_code'][role]
        if run[f'{role}_sha256'] != h:
            raise RuntimeError('Run-code attribution mismatch')
        check(E / f'{h}.py', h)
    if read(d / 'run_result.json')['status'] != 'numerical_failure':
        raise RuntimeError('Failed case classification changed')

analyses = []
for folder in ('seed_analysis', 'seed_analysis_v2'):
    d = read(E / folder / 'analysis.json')
    check(E / folder / 'driver.py', d['driver_sha256'])
    check('results/cmfgen_closure_seed/QSTARS_MOM_SEED', d['source_sha256'])
    check('results/cmfgen_closure_node/QSTARS_NEGATIVE_J',
          d['native_node_candidate_failure']['sha256'])
    analyses.append(dict(path=str(E / folder / 'driver.py'), sha256=d['driver_sha256']))
d = read(E / 'edge_realizability.json')
check(E / 'edge_analysis_driver.py', d['driver_sha256'])
check('results/cmfgen_closure_edge/QSTARS_MOM_EDGE', d['input_sha256'])
analyses.append(dict(path=str(E / 'edge_analysis_driver.py'), sha256=d['driver_sha256']))
probe = read(E / 'native_response_replayed/validation.json')
check('scripts/check_cmfgen_moment_response.py', probe['driver_sha256'])
check(E / 'native_response_replayed/driver.py', probe['driver_sha256'])
check('backends/cmfgen/moment_response_probe.f90', probe['source_sha256'])
check(E / 'native_response_replayed/probe.f90', probe['source_sha256'])
check('results/cmfgen_moment_probe_node_replay2/probe', probe['binary_sha256'])
for name, h in probe['native_sources'].items():
    check(Path('backends/cmfgen/build_12/cur_cmf') / name, h)
if not probe['passed'] or any(v['max_relative'] >= probe['target_relative']
                              for v in probe['errors'].values()):
    raise RuntimeError('Native response result is not a bounded pass')
replay = read(E / 'node_source_replay.json')
check('backends/cmfgen/apply_node_moment.py', replay['driver_sha256'])
for name in ('mom_jrel_v9.f', 'var_jrel_v5.f'):
    if digest(Path('results/cmfgen_moment_source_replay/plane') / name) != digest(
            Path('backends/cmfgen/build_12/cur_cmf/plane') / name):
        raise RuntimeError('Source replay mismatch')

# Final interface source is distinct from historical run code preserved above.
paths = subprocess.check_output(['git', 'ls-files', '--cached', '--others',
                                  '--exclude-standard', '-z']).decode().split('\0')
paths = sorted({p for p in paths if p and not p.startswith(('Reports/', 'benchmarks/results/'))}
               | {'Reports/Report_16.md'})
rows = []
with zipfile.ZipFile(E / 'interface_source.zip', 'w', zipfile.ZIP_DEFLATED) as z:
    for name in paths:
        p = Path(name)
        if not p.is_file():
            raise RuntimeError(f'Missing final source: {p}')
        z.write(p, name)
        h, n = digest(p)
        rows.append(dict(path=name, sha256=h, bytes=n))
with zipfile.ZipFile(E / 'interface_source.zip') as z:
    for row in rows:
        b = z.read(row['path'])
        if hashlib.sha256(b).hexdigest() != row['sha256'] or len(b) != row['bytes']:
            raise RuntimeError('Final source archive mismatch')
scope = ('Final Report 16 project interface, experimental patch/analysis tools, plans, tests, '
         'documentation and source registry. External native/atomic/binary assets excluded and '
         'separately attested; this final snapshot is not the historical code used by every trial.')
write(E / 'interface_source_manifest.json', dict(scope=scope,
      baseline_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode().strip(),
      files=rows, files_verified=len(rows)))

log = (E / 'pytest_native.log').read_text()
if not re.search(r'123 passed in 53\.70s', log) or 'skipped' in log or 'failed' in log:
    raise RuntimeError('Regression evidence mismatch')
registry = read('src/qstars/data/sources.json')['prescriptions']
if len(registry) != 59:
    raise RuntimeError('Source registry count mismatch')
artifact_paths = sorted(p for p in E.rglob('*') if p.is_file()
                        and p.name not in ('release_verification.json', 'verify_release.log'))
artifact_paths += [A / 'archive_manifest.json', A / 'verification.json',
                   Path('Reports/evidence/completion_gates.json')]
artifacts = [dict(path=str(p), sha256=digest(p)[0]) for p in artifact_paths]
release = dict(schema_version=1, report='Reports/Report_16.md',
    baseline_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode().strip(),
    outcome='NUMERICAL BLOCKER', software_version='0.1.0', validated_closure=False,
    self_selected_wind_pass=False, v1_complete=False, physical_nonexistence_demonstrated=False,
    regression=dict(passed=123, failed=0, skipped=0, seconds=53.70,
                    mesa_bridge_abi=6, freeeos_enabled=True),
    cli=dict(rerun=False, retained_prior_report=15, prior_passed=29, prior_checks=31,
             prior_failed=2, scope='No Report 16 CLI/science validation claim.'),
    native_cases=len(campaign['cases']),
    sum_native_process_wall_seconds=campaign['sum_process_wall_seconds'],
    source_registry_entries=len(registry), archive=local_archive, builds=builds,
    analysis_driver_verification=analyses,
    source_snapshot=dict(files_verified=len(rows), sha256=digest(E / 'interface_source.zip')[0], scope=scope),
    native_response=dict(passed=True, eta_max_relative=probe['errors']['eta']['max_relative'],
                        chi_max_relative=probe['errors']['chi']['max_relative'],
                        scope=probe['scope']), artifacts=artifacts,
    artifact_hash_scope='Verification output log excluded because it is still open during this verification; science logs and immutable analysis outputs are hashed.',
    limitations=['The retained seed is decomposed exactly; the absolute earliest full-spectrum ancestor remains open.',
                 'Common-face reconstruction is rejected; interior cancellation is not a global energy certificate.',
                 'Native eta/chi response is fixed-factor only; full relaxed nonlocal response is not validated.',
                 'No jointly converged cool radiation/population/temperature/electron/energy solution or two-resolution gate.',
                 'No physical selected wind, convective handoff or evolutionary track; Mach-one MLT guard remains.',
                 'No universal numerical impossibility, physical nonexistence or Linux/cloud agreement claim.'])
write(E / 'release_verification.json', release)
print(json.dumps(dict(archive=local_archive, builds_verified=len(builds),
                      source_files_verified=len(rows), outcome=release['outcome']), indent=2))
