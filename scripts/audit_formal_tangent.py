#!/usr/bin/env python3
"""Build and audit the analytic tangent against native CMFGEN finite differences.

Run from an initialized MESA SDK shell. Native libraries must carry build-17
formal-ALI patches; neither this driver nor the component patches the solver.
The report distinguishes this component gate from rate-equation convergence.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    'thin': (1000, 10000, 0, 2, 0.03, False, 'DIFFUSION'),
    'front': (1000, 10000, 20, 2, 0.03, False, 'DIFFUSION'),
    'inserted_front': (1000, 100, 20, 2, 0.03, False, 'DIFFUSION'),
    'thick_inserted': (1000, 100, 20, 2, 0.03, True, 'DIFFUSION'),
    'thick_absorption': (1000, 100, 100, 1e6, 1e4, True, 'DIFFUSION'),
    'scattering_dominated': (1000, 100, 20, 0.001, 0.0001, True, 'DIFFUSION'),
    'zero_flux': (1000, 100, 20, 2, 0.03, False, 'ZERO_FLUX'),
}
DIRECTIONS = {1: 'blue_emissivity', 2: 'blue_absorption', 3: 'electron_scattering',
              4: 'mixed_profile', 5: 'diffusion_boundary', 6: 'absorption_emission'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse(text):
    out = {}
    for line in text.splitlines():
        if line.startswith('RESULT '):
            _, key, value = line.split()
            out[key] = {'T': True, 'F': False}.get(value)
            if out[key] is None:
                out[key] = float(value)
    required = ({f'error_moment_{i}' for i in range(1, 6)} |
                {f'krylov_vs_dense_{i}' for i in range(1, 6)}) | {
        'scattering_tangent_closure', 'primal_replay_error',
        'scattering_fixed_point_error', 'max_ray_points', 'ray_points_inserted'}
    if not required <= out.keys():
        raise ValueError(f'incomplete native probe: {required - out.keys()}')
    return out


def gate_epsilon(case, mode):
    # Fixed after the initial 20-point step-size study, then tested unchanged
    # on the refined 39-point fixture. No minimum-error selection at runtime.
    if case == 'thick_absorption' and mode in (4, 5):
        return 1e-8  # remain on the near-zero H face-interpolation branch
    if (case == 'thick_absorption' and mode == 3 or
            case == 'scattering_dominated' and mode in (1, 2)):
        return 1e-3  # weak response: smaller differences amplify roundoff
    return 1e-5


def audit(cmf_root, output, cases=None, epsilons=(1e-3, 1e-4, 1e-5), modes=None, nd=20):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    cmf_root = cmf_root.resolve()
    build = output / 'build'
    build.mkdir(exist_ok=True)
    compiler = shutil.which('gfortran')
    if compiler is None:
        raise RuntimeError('gfortran unavailable; initialize the MESA SDK first')
    sources = [ROOT / 'backends/cmfgen' / name for name in
               ('formal_tangent.f90', 'formal_tangent_probe.f90')]
    libs = [cmf_root / 'lib' / name for name in
            ('libsubs.a', 'libmod.a', 'libmy_lapack.a', 'libmy_blas.a')]
    caps = json.loads((cmf_root / 'exe/cmfgen_dev.capabilities.json').read_text())
    if not caps.get('qstars_formal_ali'):
        raise ValueError('native build must have the recorded formal-ALI capability')
    flags = ['-cpp', f'-DNDV={nd}', '-O2', '-g', '-fcheck=bounds', '-fopenmp', '-ffree-line-length-none',
             f'-I{cmf_root / "mod"}', f'-I{build}']
    log = []
    for command in ([compiler, *flags, '-c', str(sources[0])],
                    [compiler, *flags, str(sources[1]), str(build / 'formal_tangent.o'),
                     *map(str, libs), '-framework', 'Accelerate', '-o', str(build / 'probe')]):
        r = subprocess.run(command, cwd=build, capture_output=True, text=True, check=True)
        log.append(r.stdout + r.stderr)
    (build / 'compile.log').write_text(''.join(log))
    provenance_files = sources + libs + [Path(compiler), Path(__file__)]
    provenance_files += [cmf_root / name for name in (
        'plane/cmf_formal_rel_v4.f', 'plane/solve_cmf_formal_v2.f',
        'plane/rel_variables.f', 'plane/optdepth_v4.f',
        'subs/mon_int_funs_v2.f', 'subs/mon_interp.f')]
    result = dict(scope='Fixed-geometry formal-transfer tangent component; no population Jacobian or wind certification',
                  source_hashes={str(p): sha(p) for p in provenance_files},
                  native_executable_sha256=caps['executable_sha256'],
                  native_libraries_linked=True, finite_difference_products_in_tangent=0,
                  target_relative_error=1e-6, radial_points=nd, cases={}, production_certified=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    for name in cases or CASES:
        case = CASES[name]
        rows = result['cases'][name] = {}
        for mode in modes or DIRECTIONS:
            runs = rows[DIRECTIONS[mode]] = []
            gate_eps = gate_epsilon(name, mode)
            steps = list(epsilons)
            if gate_eps == 1e-8:
                steps += [1e-6, 1e-7, 1e-8, 1e-9]
            for eps in steps:
                d = output / name / DIRECTIONS[mode] / f'eps_{eps:.0e}'
                d.mkdir(parents=True, exist_ok=True)
                inputs = ' '.join(map(str, case[:5])) + f' {"T" if case[5] else "F"} {mode} {eps} {case[6]}\n'
                start = time.monotonic()
                r = subprocess.run([str(build / 'probe')], input=inputs, cwd=d, env=env,
                                   capture_output=True, text=True, timeout=180, check=False)
                (d / 'input.txt').write_text(inputs)
                (d / 'out.txt').write_text(r.stdout + r.stderr)
                if r.returncode != 0:
                    raise RuntimeError(f'native probe exited {r.returncode}; see {d / "out.txt"}')
                p = d / 'response.txt'
                with gzip.open(d / 'response.txt.gz', 'wb') as f:
                    f.write(p.read_bytes())
                p.unlink()
                metrics = parse(r.stdout)
                metrics.update(epsilon=eps, gate_evaluation=(eps == gate_eps), seconds=time.monotonic()-start,
                               response_sha256=sha(d / 'response.txt.gz'))
                runs.append(metrics)
                print(name, DIRECTIONS[mode], eps,
                      f'error={max(metrics[f"error_moment_{i}"] for i in range(1,6)):.3e}', flush=True)
        # Retain truncation/limiter/noise failures; gate only the declared
        # response-resolving perturbation, identically on both fixture grids.
        final = [next(r for r in runs if r['gate_evaluation']) for runs in rows.values()]
        result.setdefault('gates', {})[name] = all(
            max(row[f'error_moment_{i}'] for i in range(1, 6)) < 1e-6
            and max(row[f'krylov_vs_dense_{i}'] for i in range(1, 6)) < 1e-8
            and row['scattering_tangent_closure'] < 1e-10
            and row['primal_replay_error'] < 1e-10
            and row['scattering_fixed_point_error'] < 1e-12 for row in final)
    result['pass'] = all(result['gates'].values())
    (output / 'audit.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cmf-root', type=Path, default=ROOT / 'backends/cmfgen/build_17/cur_cmf')
    p.add_argument('--output', type=Path, default=ROOT / 'results/formal_tangent/audit')
    p.add_argument('--case', action='append', choices=CASES)
    p.add_argument('--nd', type=int, default=20)
    p.add_argument('--mode', action='append', type=int, choices=DIRECTIONS)
    a = p.parse_args()
    r = audit(a.cmf_root, a.output, a.case, modes=a.mode, nd=a.nd)
    print('Component gate:', 'PASS' if r['pass'] else 'FAIL')
    raise SystemExit(0 if r['pass'] else 2)
