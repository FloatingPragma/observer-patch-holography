"""Exact integer-versus-real dual-face relaxation witness for RF1."""
from collections import deque
from functools import reduce
from math import lcm
from pathlib import Path
import hashlib
import importlib.util
import json

import sympy as s

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
GEOMETRY=ROOT/'openwave-tuning/next/ROW1_NATIVE/routes/auditor_flux_topology/audit_spatial_lift.py'
BRIEF=ROOT/'openwave-push/rows_native/RF1_BRIEF.md'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def run():
    pins={str(p.relative_to(ROOT)):sha(p) for p in (GEOMETRY,BRIEF)}
    spec=importlib.util.spec_from_file_location('frozen_integral_ico',GEOMETRY)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    edges,faces,A,B,_,certificate=module.geometry()
    assert certificate['cycle_quotient_Smith_factors']==[1]*19
    L=B*B.T;G=(L+s.ones(12)/12).inv()-s.ones(12)/12
    reduced=A[:,:19];face_inverse=(reduced.T*reduced).inv()*reduced.T
    adjacency=[[] for _ in range(12)]
    for e,(a,b) in enumerate(edges):adjacency[a].append((b,e,1));adjacency[b].append((a,e,-1))
    rows=[]
    for destination in range(1,12):
        todo=deque([0]);previous={0:None}
        while todo:
            x=todo.popleft()
            for y,e,sign in adjacency[x]:
                if y not in previous:previous[y]=(x,e,sign);todo.append(y)
        F=s.zeros(30,1);cursor=destination;distance=0
        while cursor:
            parent,e,sign=previous[cursor];F[e]=sign;cursor=parent;distance+=1
        q=B*F;optimal=face_inverse*F;real_field=F-reduced*optimal
        assert real_field==B.T*G*q
        real_energy=(real_field.T*real_field)[0]/2
        integer_energy=(F.T*F)[0]/2
        assert integer_energy==s.Rational(distance,2) and real_energy<integer_energy
        denominator=reduce(lcm,(int(v.q) for v in optimal),1)
        refinements=[]
        for k in sorted(set((1,2,4,8,16,denominator))):
            rounded=s.Matrix([s.Rational(s.floor(k*v+s.Rational(1,2)),k) for v in optimal])
            residual=F-reduced*rounded;energy=(residual.T*residual)[0]/2
            assert B*residual==q and energy>=real_energy
            assert energy-real_energy<=s.Rational(57,4*k*k)
            if k==denominator:assert energy==real_energy
            refinements.append({'face_step':'1/'+str(k),'energy':str(energy),'gap_from_real':str(energy-real_energy)})
        # Independent integer rounding with mean equal to optimal x exhibits
        # the distinction between energy of mean and mean energy exactly.
        fractions=[v-s.floor(v) for v in optimal]
        fluctuation=3*sum(v*(1-v) for v in fractions)/2
        mean_instant_energy=real_energy+fluctuation
        assert mean_instant_energy>=integer_energy
        rows.append({'source':0,'sink':destination,'distance':distance,
                     'integer_minimum_exact':str(integer_energy),'real_minimum_exact':str(real_energy),
                     'strict_gap':str(integer_energy-real_energy),
                     'integer_minimum_proof':'Any integer unit flow has squared norm at least its l1 norm, which is at least graph distance. The displayed shortest path attains distance; integral kerB=imageA reaches every other representative.',
                     'real_optimal_face_coordinates':[str(v) for v in optimal]+['0'],
                     'exact_fractional_face_denominator':denominator,'refinement_controls':refinements,
                     'independent_integer_mixture_mean_field_energy':str(real_energy),
                     'independent_integer_mixture_mean_instantaneous_energy':str(mean_instant_energy)})
    edge=next(r for r in rows if r['distance']==1)
    assert edge['integer_minimum_exact']=='1/2' and edge['real_minimum_exact']=='11/60'
    assert {r['distance']:(r['integer_minimum_exact'],r['real_minimum_exact']) for r in rows}=={
        1:('1/2','11/60'),2:('1','7/30'),3:('3/2','1/4')}
    assert pins=={str(p.relative_to(ROOT)):sha(p) for p in (GEOMETRY,BRIEF)}
    result={'schema':'oph.RF1.integer_relaxation.opposing.v1','checks':{
        'actual_oriented_icosahedral_chain_complex':True,'integral_cycle_lattice_certificate_reused':True,
        'all11_source_sink_real_minima_exact':True,'integer_minima_strictly_larger':True,
        'mean_energy_not_energy_of_mean':True,'fractional_refinement_explicit':True,'source_bytes_unchanged':True},
        'rows':rows,'scientific_consequence':'For integer face loads x, RF1 cannot identify min_x ||F-Ax||²/2 with q^T L+ q/2. That is the real-domain minimum. Time averaging can reproduce its mean field without reproducing its energy. A face-step refinement is an explicit different domain.',
        'scope':'Mathematical premise audit of the brief, not a run of an unavailable RF1 solver; no reserved acceptance cases.',
        'source_pins':pins,'auditor_sha256':sha(Path(__file__)),'official_rows_changed':[]}
    (HERE/'INTEGER_RELAXATION_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'by_distance':{r['distance']:[r['integer_minimum_exact'],r['real_minimum_exact']] for r in rows},'edge_refinement':edge['exact_fractional_face_denominator']}))


if __name__=='__main__':run()
