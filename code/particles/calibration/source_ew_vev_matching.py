#!/usr/bin/env python3
"""Conditional source-to-Landau-VEV matching, with pinned optional SMDR replay.

Default --check validates metadata, source boundary and its independently checked
one-loop consequences without a compiler/network. --smdr-dir DIR --refresh
replays all loop orders. Add --build to build the pinned upstream source locally.
No target mass, measured GF, or default/reference SMDR model is read.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
SPEC=HERE/'source_ew_vev_matching_spec.json'
RECEIPT=ROOT/'code/particles/runs/calibration/source_ew_vev_matching.json'
C_SOURCE=HERE/'source_ew_vev_matching_smdr.c'
RGE=HERE/'sm_two_loop_rge_engine.py'


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_spec(): return json.loads(SPEC.read_text())

def build_inputs(source=None):
    spec=load_spec()
    if source is None: source=json.loads((ROOT/spec['source_receipt']).read_text())
    sc=source['source_scales_GeV']
    branch=source['fully_coupled_two_loop_same_low_gauge_anchor'][1]
    q0=float(sc['mz_run_gev']); v=float(sc['v_transmutation_gev']); uv=float(sc['log_midpoint_half_turn'])
    initial=[*map(float,branch['boundary_gauges']),float(branch['boundary_top_yukawa']),0.0]
    assert q0>0 and v>0 and uv>q0 and all(math.isfinite(x) and x>=0 for x in initial)
    module_spec=importlib.util.spec_from_file_location('_source_ew_rge',RGE)
    mod=importlib.util.module_from_spec(module_spec); module_spec.loader.exec_module(mod)
    sol=solve_ivp(lambda t,x:mod.beta_2loop(*x),(math.log(uv),math.log(q0)),initial,method='DOP853',rtol=2e-12,atol=2e-14)
    if not sol.success: raise ValueError('UV transport failed')
    g1,g2,g3,yt,lam=map(float,sol.y[:,-1])
    return {'schema':'oph.smdr_conditional_matching_input.v1','Q0_GeV':q0,'v_Q0_GeV':v,'gY_Q0':g1*math.sqrt(3/5),'g2_Q0':g2,'g3_Q0':g3,'yt_Q0':yt,'lambda_Q0':lam,'non_top_yukawas':0.0,'m2_minimum_order':2,'RG_loop_order':2,'scale_multipliers':[1,2,4],'SMDR_commit':spec['upstream']['commit'],'attachment':spec['hypotheses']['vev_attachment'],'criticality_input':{'scale_GeV':uv,'initial':initial},'source_file':spec['source_receipt'],'consumed_fields':spec['source_boundary_fields']}


def source_pins():
    paths=[Path(__file__),SPEC,C_SOURCE,RGE,ROOT/load_spec()['source_receipt']]
    return {p.relative_to(ROOT).as_posix():sha(p) for p in paths}


def check_upstream(path):
    path=Path(path).resolve(); spec=load_spec()['upstream']
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=path,text=True).strip()
    if commit!=spec['commit']: raise ValueError('SMDR commit mismatch')
    names=subprocess.check_output(['git','ls-files','-z'],cwd=path).split(b'\0')
    names=sorted(n for n in names if n)
    digest=hashlib.sha256()
    for name in names:
        digest.update(name+b'\0'+hashlib.sha256((path/os.fsdecode(name)).read_bytes()).digest())
    if digest.hexdigest()!=spec['tracked_tree_sha256']: raise ValueError('SMDR tracked source bytes changed')
    for archive,key in [('tsil-1.46.tar.gz','tsil_archive_sha256'),('3vil-v2.02.tar.gz','threevil_archive_sha256')]:
        if sha(path/archive)!=spec[key]: raise ValueError('SMDR dependency archive changed')
    return path


def build_upstream(path,cc='cc'):
    path=check_upstream(path)
    # These targets are sequential: SMDR consumes both completed dependency libs.
    for target in ['tsil','3vil','smdr']:
        subprocess.run(['make',target,f'CC={cc}','SMDR_OPT=-O2'],cwd=path,check=True)


def replay(smdr_dir,inputs,cc='cc'):
    smdr_dir=check_upstream(smdr_dir)
    with tempfile.TemporaryDirectory(prefix='oph-smdr-matching-') as name:
        tmp=Path(name)
        header=(smdr_dir/'src/smdr.h').read_text()
        names=re.findall(r'extern SMDR_REAL (SMDR_\w*EXPT\w*);',header)
        if len(names)!=41: raise ValueError('Experimental-global poison coverage changed')
        (tmp/'source_ew_poison.h').write_text('\n'.join(f'{n}=NAN;' for n in names)+'\n')
        binary=tmp/'source_ew_matching'
        subprocess.run([cc,'-O2','-fcommon','-I',str(smdr_dir),'-I',str(tmp),str(C_SOURCE),'-L',str(smdr_dir),'-lsmdr','-l3vil','-ltsil','-lm','-o',str(binary)],check=True)
        args=[format(inputs[k],'.18g') for k in ['Q0_GeV','v_Q0_GeV','gY_Q0','g2_Q0','g3_Q0','yt_Q0','lambda_Q0']]
        proc=subprocess.run([str(binary),*args],check=True,capture_output=True,text=True,timeout=900)
        result=json.loads(proc.stdout)
        if proc.stderr.strip(): result['upstream_stderr']=proc.stderr.strip()
    return result


def make_receipt(inputs,forward):
    spec=load_spec()
    return {'schema':'oph.source_ew_vev_matching.v1','inputs':inputs,'specification':spec,'scope':spec['scope'],'source_pins':source_pins(),'forward':forward}


def compare(expected,actual,path='root',*,matching_order=None):
    if isinstance(expected,dict):
        if not isinstance(actual,dict) or expected.keys()!=actual.keys(): raise ValueError(f'{path}: key mismatch')
        order=expected.get('matching_order',matching_order)
        for k in expected: compare(expected[k],actual[k],f'{path}/{k}',matching_order=order)
    elif isinstance(expected,list):
        if not isinstance(actual,list) or len(expected)!=len(actual):raise ValueError(f'{path}: list mismatch')
        for i,(a,b) in enumerate(zip(expected,actual)):compare(a,b,f'{path}/{i}',matching_order=matching_order)
    elif isinstance(expected,bool) or not isinstance(expected,(int,float)):
        if expected!=actual: raise ValueError(f'{path}: metadata mismatch')
    elif isinstance(actual,bool) or not isinstance(actual,(int,float)) or not math.isfinite(actual):
        raise ValueError(f'{path}: invalid numeric type or value')
    elif isinstance(expected,int):
        if type(actual) is not int or expected!=actual:raise ValueError(f'{path}: integer mismatch')
    else:
        # GeV, couplings and GeV^-2 cannot share an absolute tolerance.
        atol=2e-7 if ('_GeV' in path or '/m2' in path or '/mt_MSbar' in path) else 2e-10
        if 'GF_removed_' in path:atol=2e-20
        elif 'GF_' in path:atol=2e-14
        elif matching_order==2 and '/top_' in path and path.endswith('/width_GeV'):
            # Only the two-loop top absorptive parts show a sub-keV spread
            # between native 64- and 80-bit long-double implementations.
            # This is a replay budget, not a certified physical error bar.
            atol=load_spec()['calculation']['two_loop_top_width_replay_atol_GeV']
        if not math.isclose(expected,actual,rel_tol=2e-8,abs_tol=atol):
            raise ValueError(f'{path}: numerical mismatch {expected} versus {actual}')


def validate_cached(receipt):
    if receipt['schema']!='oph.source_ew_vev_matching.v1':raise ValueError('schema')
    if receipt['specification']!=load_spec() or receipt['scope']!=load_spec()['scope']:raise ValueError('scope/specification promotion or drift')
    if receipt['source_pins']!=source_pins():raise ValueError('source pin mismatch')
    payload=json.dumps(receipt['forward'],sort_keys=True,separators=(',',':'),allow_nan=False).encode()
    if hashlib.sha256(payload).hexdigest()!=load_spec()['retained_forward_payload']['sha256']:
        raise ValueError('retained forward payload changed; hash is integrity, not multiloop verification')
    compare(build_inputs(),receipt['inputs'],'inputs')
    from verify_source_ew_top_pole import verify as verify_top
    verify_top(receipt)
    # The independent GF/v validator owns the corresponding one-loop checks.
    import verify_source_ew_vev_matching as independent
    independent.verify(receipt)
    return True


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--smdr-dir',type=Path);p.add_argument('--build',action='store_true');p.add_argument('--refresh',action='store_true');p.add_argument('--check',action='store_true');p.add_argument('--cc',default=os.environ.get('CC','cc'));p.add_argument('--output',type=Path,default=RECEIPT)
    a=p.parse_args()
    if a.build:
        if a.smdr_dir is None:p.error('--build requires --smdr-dir')
        build_upstream(a.smdr_dir,a.cc)
    if a.refresh:
        if a.smdr_dir is None:p.error('--refresh requires --smdr-dir')
        generated=make_receipt(build_inputs(),replay(a.smdr_dir,build_inputs(),a.cc))
        if a.check:
            cached=json.loads(a.output.read_text());compare(cached,generated);validate_cached(cached)
        else:a.output.write_text(json.dumps(generated,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif a.check:validate_cached(json.loads(a.output.read_text()))
    else:p.error('choose --check or --refresh')
    print('Source EW VEV matching: PASS')

if __name__=='__main__':main()
