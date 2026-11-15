#!/usr/bin/env python3
"""Exact paired native current moments after current-only restoring censuses."""
from __future__ import annotations
from collections import Counter
from fractions import Fraction as F
import argparse
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INPUT = ROOT/'openwave/row_01_charge_quantization/work/2026-11-14_codex_relative_source'
SEED = 711701
ZERO = (F(0),F(0))
def q(a=0,b=0): return F(a),F(b)
def add(x,y): return x[0]+y[0],x[1]+y[1]
def neg(x): return -x[0],-x[1]
def sub(x,y): return add(x,neg(y))
def mul(x,y):
    a,b=x; c,d=y
    return a*c+b*d,a*d+b*c+b*d

def div(x,y):
    a,b=y; denominator=a*a+a*b-b*b
    if not denominator: raise ZeroDivisionError
    return mul(x,((a+b)/denominator,-b/denominator))

def total(xs):
    out=ZERO
    for x in xs: out=add(out,x)
    return out

def dot(x,y): return total(mul(a,b) for a,b in zip(x,y))
def outer(x): return [[mul(a,b) for b in x] for a in x]
def matvec(M,x): return [dot(row,x) for row in M]
def msub(A,B): return [[sub(x,y) for x,y in zip(ar,br)] for ar,br in zip(A,B)]
def iszero(obj):
    if isinstance(obj,tuple): return obj==ZERO
    return all(iszero(x) for x in obj)

def encode(obj):
    if isinstance(obj,tuple) and len(obj)==2 and all(isinstance(x,F) for x in obj):
        return [str(obj[0]),str(obj[1])]
    if isinstance(obj,(list,tuple)): return [encode(x) for x in obj]
    if isinstance(obj,dict): return {str(k):encode(v) for k,v in obj.items()}
    return obj

def canonical(x): return json.dumps(encode(x),sort_keys=True,separators=(',',':'))+'\n'
def sha(b): return hashlib.sha256(b).hexdigest()
def approx(x): return float(x[0])+float(x[1])*(1+5**.5)/2

def load_source():
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for name,digest in pins['files'].items():
        assert sha((ROOT/name).read_bytes())==digest, name
    spec=importlib.util.spec_from_file_location('paired_moment_charge',INPUT/'headline.py')
    mod=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=mod
    spec.loader.exec_module(mod)
    return mod


def geometry(seams):
    z=ZERO
    points=[(q(-1),q(0,1),z),(q(1),q(0,1),z),(q(-1),q(0,-1),z),(q(1),q(0,-1),z),
            (z,q(-1),q(0,1)),(z,q(1),q(0,1)),(z,q(-1),q(0,-1)),(z,q(1),q(0,-1)),
            (q(0,1),z,q(-1)),(q(0,1),z,q(1)),(q(0,-1),z,q(-1)),(q(0,-1),z,q(1))]
    actual=tuple((a,b) for a in range(12) for b in range(a+1,12) if dot(points[a],points[b])==q(0,1))
    assert actual==seams, 'geometry must use actual row1 seam labels'
    antipode=[points.index(tuple(neg(v) for v in p)) for p in points]
    rotated=[tuple(mul(q(s),v) for s,v in zip((-1,-1,1),p)) for p in points]
    permutation=[points.index(p) for p in rotated]
    assert permutation[0]==antipode[0]==3
    assert all(permutation[permutation[p]]==p for p in range(12))
    assert all(tuple(sorted((permutation[a],permutation[b]))) in seams for a,b in seams)
    # Same explicit port vectors as the Lean finite current quotient.
    lean=[(z,q(-1),q(0,-1)),(q(-1),q(0,-1),z),(q(0,-1),z,q(-1)),(q(1),q(0,-1),z),
          (z,q(1),q(0,-1)),(q(0,-1),z,q(1)),(q(0,1),z,q(-1)),(z,q(-1),q(0,1)),
          (q(-1),q(0,1),z),(q(0,1),z,q(1)),(q(1),q(0,1),z),(z,q(1),q(0,1))]
    return points,antipode,permutation,[lean.index(p) for p in points]


def moments(vectors):
    n=q(len(vectors))
    mean=[div(total(v[k] for v in vectors),n) for k in range(3)]
    second=[[div(total(mul(v[i],v[j]) for v in vectors),n) for j in range(3)] for i in range(3)]
    return {'mean':mean,'second_moment':second,'covariance':msub(second,outer(mean))}


def rotate_vector(v): return [mul(q(s),a) for s,a in zip((-1,-1,1),v)]
def rotate_tensor(C):
    return [[mul(q((-1,-1,1)[i]*(-1,-1,1)[j]),C[i][j]) for j in range(3)] for i in range(3)]


def analyze(stat,axis):
    m,C=stat['mean'],stat['covariance']; r2=dot(axis,axis)
    parallel=div(dot(axis,matvec(C,axis)),r2)
    transverse=div(sub(total(C[k][k] for k in range(3)),parallel),q(2))
    fit=[[add(mul(transverse,q(int(i==j))),mul(div(sub(parallel,transverse),r2),mul(axis[i],axis[j])))
          for j in range(3)] for i in range(3)]
    alignment=[sub(mul(r2,m[k]),mul(dot(m,axis),axis[k])) for k in range(3)]
    return {'mean_nonzero':not iszero(m),'covariance_nonzero':not iszero(C),
            'mean_axis_residual':alignment,'uniaxial_residual':msub(C,fit),
            'parallel_variance':parallel,'transverse_variance':transverse,
            'anisotropy':sub(parallel,transverse),
            'native_axis_predicate':not iszero(m) and sub(parallel,transverse)!=ZERO and iszero(alignment) and iszero(msub(C,fit))}


def rotation_detector(base,turned):
    """The identical comparison is used on native, controls and the plant."""
    covariance_equivariant=turned['covariance']==rotate_tensor(base['covariance'])
    tensor_returns=turned['covariance']==base['covariance']
    frame_rotates=turned['mean']==rotate_vector(base['mean'])
    frame_reverses=turned['mean']==[neg(v) for v in base['mean']]
    return {'covariance_equivariant':covariance_equivariant,'tensor_returns':tensor_returns,
            'frame_rotates':frame_rotates,'frame_reverses':frame_reverses,
            'frame_differs':turned['mean']!=base['mean'],
            'geometric_half_turn':covariance_equivariant and tensor_returns and frame_rotates and frame_reverses and turned['mean']!=base['mean']}


def branch_records(mod,x,y,points,law):
    records=[]; audit=[]; vectors=[]
    for e,(a,b) in enumerate(mod.SEAMS):
        for coin in (0,1):
            if law=='native':
                ox,(_,jx)=mod.native(x,(a,b),coin,'paired_moment_branches')
                oy,(_,jy)=mod.native(y,(a,b),coin,'paired_moment_branches')
            else:
                active=law=='fixed_seam' and (a,b)==(0,5)
                jx=(x[a]-x[b]) if active else 0
                jy=(y[a]-y[b]) if active else 0
                ox,oy=list(x),list(y)
                ox[a]-=jx; ox[b]+=jx; oy[a]-=jy; oy[b]+=jy
            delta=[mul(q(jx-jy),sub(points[b][k],points[a][k])) for k in range(3)]
            observed=list(mod.difference(x,y)); observed[a]-=jx-jy; observed[b]+=jx-jy
            truth=tuple(observed)==mod.difference(ox,oy)
            records.append([e,coin,jx,jy])
            vectors.append(delta)
            audit.append(truth)
    return records,vectors,audit


def run():
    mod=load_source(); mod.CALLS.clear()
    points,antipodes,rotation,lean_reindex=geometry(mod.SEAMS)
    public=[]; private=[]; by_orientation={law:{} for law in ('native','frozen','fixed_seam')}
    groups={law:[] for law in by_orientation}
    main_counts=Counter(); failures=Counter(); known_axes=set()
    native_vectors={}; audit_all=True; input_count=0
    for p in range(12):
        r=antipodes[p]
        for s in range(12):
            if s in (p,r): continue
            ident=input_count; input_count+=1
            x,y=mod.pair_state(p,s),mod.pair_state(r,s)
            rng=random.Random(SEED+ident)
            hx,cx=mod.census(x,rng); hy,cy=mod.census(y,rng)
            record={'trial':ident,'seed':SEED+ident,'census':[cx,cy]}
            truth=(cx['outcome']!='restored' or (tuple(cx['initial'])==x and hx==x)) and (cy['outcome']!='restored' or (tuple(cy['initial'])==y and hy==y))
            audit={'trial':ident,'input_x':x,'input_y':y,'census_truth':truth,'body_truth':{}}
            audit_all &= truth
            if cx['outcome']!='restored' or cy['outcome']!='restored':
                record['outcome']='census_failure'; main_counts['census_failure']+=1
                failures[cx['outcome']+'|'+cy['outcome']]+=1
            else:
                dx,dy=tuple(cx['initial']),tuple(cy['initial'])
                observed=mod.difference(dx,dy)
                assert sorted(observed)==[-1]+[0]*10+[1]
                a=observed.index(1); b=observed.index(-1)
                assert b==antipodes[a]
                record.update(outcome='restored',source=observed,source_axis=a,branches={})
                known_axes.add(a); main_counts['restored']+=1
                # Each branch starts from the same restored preparation. This is
                # a complete instrument calculation, not a sequential trajectory.
                for law in by_orientation:
                    rows,vec,truths=branch_records(mod,hx,hy,points,law)
                    record['branches'][law]=rows
                    audit['body_truth'][law]=truths
                    audit_all &= all(truths)
                    stat=moments(vec); analysis=analyze(stat,points[a])
                    groups[law].append({'trial':ident,'source_axis':a,**stat,**analysis})
                    by_orientation[law].setdefault(a,[]).append(stat)
                    if law=='native': native_vectors.setdefault(a,vec)
            public.append(record); private.append(audit)
    all_axes=len(known_axes)==12
    base={law:by_orientation[law][0][0] for law in by_orientation} if all_axes else {}
    turned={law:by_orientation[law][3][0] for law in by_orientation} if all_axes else {}
    detections={law:rotation_detector(base[law],turned[law]) for law in base}
    # Off-axis plant traverses the same raw-vector moment and rotation reader.
    plant_vectors=[[q(1),q(1),q(1)],[q(0),q(0),q(0)]]
    plant_base=moments(plant_vectors)
    plant_turned=moments([[neg(x) for x in v] for v in plant_vectors])
    plant_detector=rotation_detector(plant_base,plant_turned)
    pooled=moments([v for a in sorted(native_vectors) for v in native_vectors[a]])
    expected_mean=q(F(-3,15),F(1,15))
    expected_parallel=q(F(3,9),F(-1,9)); expected_transverse=q(F(2,15),F(1,15))
    checks={
        'all120_inputs_retained':input_count==len(public)==120,
        'at_least100_restored':main_counts['restored']>=100,
        'all12_source_orientations_read':all_axes,
        'all_census_and_continuity_audits':audit_all,
        'complete60branches_per_law_per_restored_input':all(len(t['branches'][law])==60 for t in public if t['outcome']=='restored' for law in by_orientation),
        'native_moments_independent_of_background':all(all(stat==stats[0] for stat in stats) for stats in by_orientation['native'].values()),
        'native_mean_exact_formula':all(d['mean']==[mul(expected_mean,v) for v in points[d['source_axis']]] for d in groups['native']),
        'native_covariance_exact_formula':all(d['parallel_variance']==expected_parallel and d['transverse_variance']==expected_transverse for d in groups['native']),
        'native_source_axis_predicate_all':all(d['native_axis_predicate'] for d in groups['native']),
        'native_geometric_half_turn':detections.get('native',{}).get('geometric_half_turn',False),
        'two_half_turns_return_full_frame':all(rotation[rotation[p]]==p for p in range(12)),
        'frozen_control_fails_nonzero_axis_predicate':all(not d['native_axis_predicate'] and iszero(d['mean']) and iszero(d['covariance']) for d in groups['frozen']),
        'fixed_seam_breaks_source_axis_alignment':any(not iszero(d['mean_axis_residual']) for d in groups['fixed_seam']),
        'fixed_seam_breaks_geometric_tensor_return':not detections.get('fixed_seam',{}).get('covariance_equivariant',True),
        'planted_off_axis_tensor_detected_same_reader':not plant_detector['covariance_equivariant'] and not plant_detector['geometric_half_turn'],
        'pooled_source_orientation_erases_axis':iszero(pooled['mean']) and pooled['covariance']==[[q(F(2,9)*int(i==j)) for j in range(3)] for i in range(3)]
    }
    records_text=''.join(canonical(t) for t in public)
    audit_text=''.join(canonical(t) for t in private)
    counts=dict(sorted(mod.CALLS.items()))
    result={'schema':'oph.paired_source_moments.v1','candidate':'A-paired-charge+moment','class':'N+',
            'checks':checks,'all_checks_pass':all(checks.values()),'counts':dict(main_counts),'census_failure_types':dict(failures),
            'native_calls':counts,'total_native_calls':sum(counts.values()),
            'geometry':{'points':points,'antipodes':antipodes,'half_turn':rotation,'lean_port_reindex':lean_reindex},
            'statistics':groups,'rotation_detector':detections,'pooled':pooled,
            'planted_control':{'base':plant_base,'turned':plant_turned,'detector':plant_detector},
            'exact_formula':{'mean_coefficient':expected_mean,'parallel_variance':expected_parallel,'transverse_variance':expected_transverse},
            'decimal_formula':{'mean_coefficient':approx(expected_mean),'parallel_variance':approx(expected_parallel),'transverse_variance':approx(expected_transverse)},
            'records_sha256':sha(records_text.encode()),'audit_sha256':sha(audit_text.encode()),
            'scope':'Static tensor and oriented-current records of the same measured primitive relative source as row1. Finite proper port rotation is executed; continuum angles describe its fitted representation. No physical angular momentum or particle dynamics claim.'}
    return result,records_text,audit_text


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true'); args=parser.parse_args()
    result,records,audit=run(); out=canonical(result)
    if args.write:
        for name,data in [('RECEIPT.json',out),('RECORDS.jsonl.gz',records),('PRIVATE_AUDIT.jsonl.gz',audit)]:
            p=HERE/name
            if name.endswith('.gz'):
                with p.open('xb') as f: f.write(gzip.compress(data.encode(),mtime=0))
            else:
                with p.open('x') as f: f.write(data)
    else:
        assert (HERE/'RECEIPT.json').read_text()==out, 'receipt mismatch'
        assert gzip.decompress((HERE/'RECORDS.jsonl.gz').read_bytes()).decode()==records
        assert gzip.decompress((HERE/'PRIVATE_AUDIT.jsonl.gz').read_bytes()).decode()==audit
    print(json.dumps({'all_checks_pass':result['all_checks_pass'],'checks':result['checks'],'counts':result['counts'],
                      'total_native_calls':result['total_native_calls'],'receipt_sha256':sha(out.encode())},indent=2))
    return 0 if result['all_checks_pass'] else 1

if __name__=='__main__': raise SystemExit(main())
