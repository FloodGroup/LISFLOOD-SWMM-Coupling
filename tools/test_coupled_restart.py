"""Regression checks for fixed-interval FV1 coupling and pending-volume restart."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import time


def grid(path):
    values=[float(v) for line in path.read_text().splitlines()[6:] for v in line.split()]
    assert values and all(math.isfinite(v) and v>=0 for v in values),path
    return values


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--exe',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--cuda',action='store_true')
    ap.add_argument('--runtime-dir',type=Path)
    ap.add_argument('--rain',action='store_true',help='Add positive and negative NetCDF net rain')
    args=ap.parse_args()
    exe=args.exe.resolve(strict=True); out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    env={k.upper():v for k,v in os.environ.items()}; env['OMP_NUM_THREADS']='4'
    if args.runtime_dir: env['PATH']=str(args.runtime_dir.resolve())+os.pathsep+env.get('PATH','')
    report={'exe_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'net_rain':args.rain,'tests':[],'rejections':[]}

    def setup(case,solver,gpu,mode='2',resume=None,write=True):
        case.mkdir()
        header='ncols 10\nnrows 10\nxllcorner 0\nyllcorner 0\ncellsize 10\nNODATA_value -9999\n'
        (case/'dem.asc').write_text(header+(' '.join(['10']*10)+'\n')*10)
        (case/'start.asc').write_text(header+(' '.join(['0.1']*10)+'\n')*10)
        (case/'flood.config').write_text('node_num 16\nID Index C0 A0 Qout Cw W Qin C1 A1\n'+
            '\n'.join(f'J{i} {i-1} 0.6 1 0 0.5 1 1 0.6 1' for i in range(1,17))+'\n')
        network='''[OPTIONS]
FLOW_UNITS CMS
FLOW_ROUTING DYNWAVE
START_DATE 01/01/2020
START_TIME 00:00:00
REPORT_START_DATE 01/01/2020
REPORT_START_TIME 00:00:00
END_DATE 01/01/2020
END_TIME 00:00:30
REPORT_STEP 00:00:01
ROUTING_STEP 0.2
VARIABLE_STEP 0.75
MINIMUM_STEP 0.01
MAX_TRIALS 30
HEAD_TOLERANCE 0.0015
THREADS 4
ALLOW_PONDING NO
[JUNCTIONS]
J1 9 1 1 0 0
[OUTFALLS]
OF1 8.5 COUPLED NO
[CONDUITS]
C1 J1 OF1 20 0.013 0 0 0 0
[XSECTIONS]
C1 CIRCULAR 0.5 0 0 0 1
[INFLOWS]
J1 FLOW constant FLOW 1.0 1.0
[TIMESERIES]
constant 0:00 1.0
constant 0:01 1.0
[COORDINATES]
J1 45 45
OF1 55 45
[REPORT]
NODES ALL
LINKS ALL
'''
        network=network.replace('C1 J1 OF1 20 0.013 0 0 0 0',
            '\n'.join(f'C{i} J{i} OF{i} 20 0.013 0 0 0 0' for i in range(1,17)))
        network=network.replace('C1 CIRCULAR 0.5 0 0 0 1',
            '\n'.join(f'C{i} CIRCULAR 0.5 0 0 0 1' for i in range(1,17)))
        network=network.replace('J1 9 1 1 0 0','\n'.join(f'J{i} 9 1 1 0 0' for i in range(1,17)))
        network=network.replace('OF1 8.5 COUPLED NO','\n'.join(f'OF{i} 8.5 COUPLED NO' for i in range(1,17)))
        network=network.replace('J1 FLOW constant FLOW 1.0 1.0',
            '\n'.join(f'J{i} FLOW constant FLOW 1.0 1.0' for i in range(1,17)))
        network=network.replace('J1 45 45\nOF1 55 45',
            '\n'.join(f'J{i} 45 45' for i in range(1,17))+'\n'+
            '\n'.join(f'OF{i} 55 45' for i in range(1,17)))
        (case/'network.inp').write_text(network)
        par=('DEMfile dem.asc\nstartfile start.asc\nresroot check\ndirroot results\n'
             'sim_time 12\ninitial_tstep 0.2\nmassint 0.5\nsaveint 1\n'
             'output_precision 9\nelevoff\n'+solver+'\n'+('cuda\n' if gpu else '')+
             'swmminp network.inp\nswmmrpt network.rpt\nswmmout network.out\n'
             'floodconfig flood.config\nflooderror mapping.log\ncoupleint '+mode+'\n')
        if write: par+='hsfdt 5.25\nhsffloodtemp saved.hsf\n'
        if resume: par+='hsffloodread '+str(resume.resolve())+'\n'
        if args.rain:
            from netCDF4 import Dataset
            with Dataset(case/'rain.nc','w') as nc:
                nc.createDimension('time',4); nc.createDimension('x',10); nc.createDimension('y',10)
                t=nc.createVariable('time','f8',('time',));t.units='hours';t[:]=[0,4/3600,8/3600,12/3600]
                nc.createVariable('x','f8',('x',))[:]=list(range(5,100,10))
                nc.createVariable('y','f8',('y',))[:]=list(range(95,0,-10))
                rain=nc.createVariable('rainfall_depth','f4',('time','y','x'));rain.units='mm'
                for i,mm in enumerate((2,-1,0,0)):rain[i,:,:]=mm
            par+='dynamicrainfile rain.nc\n'
        (case/'case.par').write_text(par)

    def run(case,fail=False):
        started=time.perf_counter()
        with (case/'console.log').open('w') as log:
            result=subprocess.run([str(exe),'case.par'],cwd=case,env=env,
                stdout=log,stderr=subprocess.STDOUT,timeout=180)
        log=(case/'console.log').read_text(errors='replace')
        if fail:
            assert result.returncode!=0,(case,result.returncode)
            assert 'restart' in log.lower(),log[-2000:]
            report['rejections'].append({'case':case.name,'exit_code':result.returncode})
            return log
        assert result.returncode==0,(case,result.returncode,log[-5000:])
        audit_log=case/'results/check.audit.log'
        diagnostics=log+(audit_log.read_text(errors='replace') if audit_log.exists() else '')
        assert re.search(r'SWMM threads 4|swmm_threads=4',diagnostics),log[-3000:]
        files=sorted((case/'results').glob('*.wd'))
        assert files,case
        for p in files: grid(p)
        row={'case':case.name,'exit_code':result.returncode,'wall_seconds':time.perf_counter()-started}
        if (case/'results/check.water_balance.csv').exists():
            with (case/'results/check.water_balance.csv').open() as stream:
                rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(stream)]
            assert rows and all(math.isfinite(v) for r in rows for v in r.values())
            assert abs(rows[-1]['time_2d_s']-12)<1.e-5,rows[-1]
            row['max_balance_residual_m3']=max(abs(r['combined_balance_residual_m3']) for r in rows)
            assert row['max_balance_residual_m3']<.02,row
            interface=[float(x) for x in re.findall(r'INTERFACE [^\n]* residual_m3=([\deE+.\-]+)',diagnostics)]
            assert interface,case
            row['max_interface_balance_residual_m3']=max(map(abs,interface))
            assert row['max_interface_balance_residual_m3']<.002,row
            row['final_pending_to_2d_m3']=rows[-1]['pending_to_2d_m3']
            assert abs(row['final_pending_to_2d_m3'])<1.e-7,row
            if 'cuda' in case.name: assert 'Using CUDA device:' in log
            assert 'swmm_wall_s' in rows[-1] and rows[-1]['swmm_wall_s']>0,rows[-1]
            if args.rain:
                assert rows[-1]['netrain_removed_m3']>0,rows[-1]
        report['tests'].append(row)
        print(json.dumps(row),flush=True)
        return log

    pattern=r'RESTART (?:saved|restored) time_s=([\deE+.\-]+) pending_to_1d_m3=([\deE+.\-]+) pending_to_2d_m3=([\deE+.\-]+) exchange_elapsed_s=([\deE+.\-]+)'
    for gpu in ([False,True] if args.cuda else [False]):
        for solver in ('acceleration','fv1'):
            prefix=('cuda_' if gpu else 'cpu_')+solver
            full=out/(prefix+'_full'); setup(full,solver,gpu)
            log=run(full)
            state=tuple(map(float,re.search(pattern,log).groups()))
            assert state[2]>1.e-6,('Test did not create pending return water',state)
            assert state[1]>1.e-6 and state[3]>1.e-6,state
            report['tests'][-1]['saved_time_s']=state[0]
            report['tests'][-1]['saved_pending_to_1d_m3']=state[1]
            report['tests'][-1]['saved_pending_to_2d_m3']=state[2]
            report['tests'][-1]['saved_exchange_elapsed_s']=state[3]
            resumed=out/(prefix+'_resumed'); setup(resumed,solver,gpu,resume=full/'saved.hsf',write=False)
            rlog=run(resumed)
            restored=tuple(map(float,re.search(pattern,rlog).groups()))
            assert state==restored,(state,restored)
            for name in ('check.max',sorted(p.name for p in (full/'results').glob('*.wd'))[-1]):
                a,b=grid(full/'results'/name),grid(resumed/'results'/name)
                difference=max(abs(x-y) for x,y in zip(a,b))
                assert difference<2.e-6,(prefix,name,difference)
                report['tests'][-1][name+'_max_difference_m']=difference
            step=out/(prefix+'_step'); setup(step,solver,gpu,mode='step',write=False);run(step)
            wrong=out/(prefix+'_wrong_interval');setup(wrong,solver,gpu,mode='3',resume=full/'saved.hsf',write=False);run(wrong,True)
            other=out/(prefix+'_wrong_solver')
            setup(other,'fv1' if solver=='acceleration' else 'acceleration',gpu,resume=full/'saved.hsf',write=False)
            run(other,True)
            short=out/(prefix+'_truncated');setup(short,solver,gpu,resume=short/'bad.hsf',write=False)
            (short/'bad.hsf').write_bytes((full/'saved.hsf').read_bytes()[:-24]);run(short,True)
            unknown=out/(prefix+'_unsupported_format')
            setup(unknown,solver,gpu,resume=unknown/'bad.hsf',write=False)
            (unknown/'bad.hsf').write_bytes(b'LFSWM???'+(full/'saved.hsf').read_bytes()[8:])
            run(unknown,True)
    report['passed']=True
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print('All checks passed',flush=True)


if __name__=='__main__':main()
