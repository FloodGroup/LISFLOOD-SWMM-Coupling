"""Short synthetic regression checks; writes only to a new output directory."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time


def read_depth(path):
    lines = path.read_text().splitlines()
    header = dict(line.lower().split() for line in lines[:6])
    values = [float(v) for line in lines[6:] for v in line.split()]
    assert len(values) == int(header['ncols']) * int(header['nrows']), path
    values = [v for v in values if v != float(header['nodata_value'])]
    assert values and all(math.isfinite(v) and v >= 0 for v in values), path
    return values, float(header['cellsize']) ** 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cuda', action='store_true')
    parser.add_argument('--runtime-dir', type=Path)
    args = parser.parse_args()
    exe = args.exe.resolve(strict=True)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    template = Path(__file__).resolve().parents[1] / 'examples/02_coupled_synthetic'
    env = {k.upper(): v for k, v in os.environ.items()}
    env['OMP_NUM_THREADS'] = '4'
    if args.runtime_dir:
        env['PATH'] = str(args.runtime_dir.resolve()) + os.pathsep + env.get('PATH', '')
    results = []

    def run(case, gpu, duration):
        started = time.perf_counter()
        with (case / 'console.log').open('w') as log:
            process = subprocess.run([str(exe), 'case.par'], cwd=case, env=env,
                                     stdout=log, stderr=subprocess.STDOUT, timeout=180)
        console = (case / 'console.log').read_text(errors='replace')
        assert process.returncode == 0, (case.name, process.returncode, console[-3000:])
        if gpu:
            assert 'Using CUDA device:' in console, case
        depths = sorted((case / 'results').glob('*.wd'))
        assert depths, case
        for depth in depths:
            read_depth(depth)
        mass_file, = (case / 'results').glob('*.mass')
        mass_rows = [[float(v) for v in line.split()]
                     for line in mass_file.read_text().splitlines()[1:] if line.strip()]
        assert mass_rows and all(math.isfinite(v) for r in mass_rows for v in r), case
        assert abs(mass_rows[-1][0] - duration) < .001, case
        row = {'case': case.name, 'exit_code': process.returncode,
               'wall_seconds': time.perf_counter() - started}
        return console, row, read_depth(depths[-1])

    for gpu in ([False, True] if args.cuda else [False]):
        for solver in ('acceleration', 'fv1'):
            case = output / (('cuda_' if gpu else 'cpu_') + solver)
            case.mkdir()
            for name in ('dem.asc', 'network.inp', 'floodconfig.txt'):
                shutil.copy2(template / name, case / name)
            header = (template / 'start.asc').read_text().splitlines()[:6]
            (case / 'start.asc').write_text('\n'.join(header + [' '.join(['0.1'] * 10)] * 10) + '\n')
            (case / 'case.par').write_text(
                'DEMfile dem.asc\nstartfile start.asc\nresroot syn\ndirroot results\n'
                'sim_time 60\ninitial_tstep 0.2\nmassint 5\nsaveint 10\n'
                'output_precision 9\nelevoff\n' + solver + '\n' + ('cuda\n' if gpu else '') +
                'swmminp network.inp\nswmmrpt network.rpt\nswmmout network.out\n'
                'floodconfig floodconfig.txt\nflooderror flood_error.txt\ncoupleint step\n')
            console, row, (depth, area) = run(case, gpu, 60)
            row['final_surface_volume_m3'] = sum(depth) * area
            assert abs(row['final_surface_volume_m3'] - 1000) > .001, row
            if gpu:
                with (case / 'results/syn.water_balance.csv').open() as stream:
                    audit = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(stream)]
                assert audit and all(math.isfinite(v) for r in audit for v in r.values())
                assert abs(audit[-1]['time_2d_s'] - 60) < .001
                assert audit[-1]['inlet_to_1d_m3'] > 0
                assert audit[-1]['outfall_to_1d_m3'] + audit[-1]['outfall_to_2d_m3'] > 0
                row['max_balance_residual_m3'] = max(abs(r['combined_balance_residual_m3']) for r in audit)
                assert row['max_balance_residual_m3'] < .02, row
            results.append(row)
            print(json.dumps(row), flush=True)

        case = output / ('boundary_cuda_acc' if gpu else 'boundary_cpu_acc')
        case.mkdir()
        shutil.copy2(template / 'dem.asc', case / 'dem.asc')
        (case / 'case.par').write_text(
            'DEMfile dem.asc\nresroot pts\ndirroot results\nsim_time 2\n'
            'initial_tstep 0.1\nmassint 1\nsaveint 1\noutput_precision 9\n'
            'elevoff\nacceleration\nbcifile points.bci\nbdyfile points.bdy\n' + ('cuda\n' if gpu else ''))
        (case / 'points.bci').write_text(
            'P 45 45 QVAR inflow_a\nP 45 45 QVAR inflow_b\n'
            'P -1 45 QFIX 100\nP 65 65 QVAR missing\n')
        (case / 'points.bdy').write_text(
            'Synthetic boundary regression\ninflow_a\n2 seconds\n0.01 0\n0.01 10\n'
            'inflow_b\n2 seconds\n0.02 0\n0.02 10\n')
        console, row, (depth, area) = run(case, gpu, 2)
        assert 'outside the DEM' in console and 'were disabled' in console
        row['surface_volume_m3'] = sum(depth) * area
        assert abs(row['surface_volume_m3'] - .6) < .0001, row
        results.append(row)
        print(json.dumps(row), flush=True)

    report = {'passed': True, 'exe_sha256': hashlib.sha256(exe.read_bytes()).hexdigest(),
              'threads': 4, 'tests': results}
    (output / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
