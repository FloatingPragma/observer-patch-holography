"""Executed source-law controls supplementing the frozen native row-1 package.

Development evidence only. Geometry and the actual endpoint implementation are
inputs; source histories are enumerated here, not read from archived outputs.
The contour samples are a diagnostic of a declared encoding of the current.
"""
from collections import deque
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import argparse
import cmath
import hashlib
import itertools
import json
import math
import random
import sys

HERE = Path(__file__).resolve().parent
PRIMARY = HERE.parent / 'openwave_row1_native'
sys.path.insert(0, str(PRIMARY))
import native_charge_winding as spatial
from native_kernel import native

# Chosen before execution, independent of the observed tape. Development
# schedules have <=3 steps, |j|<=1 in the native and random-integer arms.
SAMPLES_PER_EDGE = 32
DECLARED_TURN_BOUND = 3


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def schedules(edges):
    """One breadth-first shortest schedule per ordered pair, seam-index order."""
    neighbors = [[] for _ in range(12)]
    for e, (a, b) in enumerate(edges):
        neighbors[a].append((b, e)); neighbors[b].append((a, e))
    for origin in range(12):
        paths = {origin: ()}
        queue = deque([origin])
        while queue:
            a = queue.popleft()
            for b, e in neighbors[a]:
                if b not in paths:
                    paths[b] = paths[a] + (e,)
                    queue.append(b)
        for destination in range(12):
            if destination != origin:
                yield origin, destination, paths[destination]


@lru_cache(maxsize=None)
def curve_fixed(first, second, turns, deformed=False):
    ga = .31 * math.sin(first + 1) if deformed else 0.
    gb = .31 * math.sin(second + 1) if deformed else 0.
    amplitude = .23 * math.cos(first + second + 1) if deformed else 0.
    # This bound is the admitted domain, not an adaptive sampling density.
    require(abs(turns) <= DECLARED_TURN_BOUND, 'outside fixed sampler domain')
    derivative = abs(gb-ga) + 2*math.pi*(DECLARED_TURN_BOUND + abs(amplitude))
    require(derivative / SAMPLES_PER_EDGE < math.pi/2, 'insufficient fixed density')
    values = [cmath.exp(1j*(ga+(gb-ga)*k/SAMPLES_PER_EDGE
              +2*math.pi*float(turns)*k/SAMPLES_PER_EDGE
              +amplitude*math.sin(2*math.pi*k/SAMPLES_PER_EDGE)))
              for k in range(SAMPLES_PER_EDGE+1)]
    defect = abs(values[-1]-cmath.exp(1j*gb))
    return values, defect


def phase_readout(tape, geometry, deformed=False):
    dual_edges, contours, signs = geometry
    curves = [curve_fixed(a,b,signs[e]*tape[e],deformed)
              for e,(a,b) in enumerate(dual_edges)]
    defect = max(d for _,d in curves)
    if defect > 1e-10:
        return None, defect
    degrees = []
    for port in range(12):
        walk = contours[port]
        samples = []
        for a,b in zip(walk,walk[1:]+walk[:1]):
            pair = tuple(sorted((a,b)))
            values = curves[dual_edges.index(pair)][0]
            if a>b:
                values = list(reversed(values))
            require(not samples or abs(samples[-1]-values[0])<1e-10,
                    'unshared vertex phase')
            samples.extend(values if not samples else values[1:])
        degrees.append(spatial.winding(samples))
    return degrees, defect


def divergence(tape, edges):
    out = [F(0) for _ in range(12)]
    for j,(a,b) in zip(tape,edges):
        out[a] += j; out[b] -= j
    return out


def one_unit(loads):
    return (sum(loads)==1 and all(x.denominator==1 and x>=0 for x in loads)
            and sum(x==1 for x in loads)==1)


def regional_bound(q):
    # Singleton integrality and extremal subsets certify membership in
    # {-1,0,1} for all 4096 regions, not just a real interval bound.
    return (all(x.denominator==1 for x in q) and
            sum(x for x in q if x>0)<=1 and sum(x for x in q if x<0)>=-1)


def run_history(origin, path, coins, edges, law, rng):
    initial = [F(int(v==origin)) for v in range(12)]
    loads = initial.copy(); tape = [F(0) for _ in edges]
    held = tape.copy(); omitted = False
    prefixes = []; events = []
    for e,c in zip(path,coins):
        a,b = edges[e]; before = loads.copy()
        if law in ('native','omitted_record'):
            actual, (_,j) = native(tuple(int(x) for x in loads),(a,b),c)
            total = loads[a]+loads[b]
            left = total//2+c*(int(total)%2)
            require(actual[a]==left and actual[b]==total-left, 'endpoint law differs')
            j=F(j)
        elif law=='half_current':
            j=(loads[a]-loads[b])/2
        elif law=='random_integer':
            j=F(rng.choice((-1,0,1)))
        else:
            raise ValueError(law)
        loads[a]-=j; loads[b]+=j; tape[e]+=j
        if law=='omitted_record' and j and not omitted:
            omitted=True
        else:
            held[e]+=j
        q=[g-n for g,n in zip(initial,loads)]
        prefixes.append({'continuity':divergence(tape,edges)==q,
                         'one_unit':one_unit(loads), 'regional_bound':regional_bound(q)})
        events.append({'seam_index':e,'coin':c,'current':str(j),
                       'before':[str(x) for x in before], 'after':[str(x) for x in loads]})
    return initial,loads,tape,held,prefixes,events,omitted


def execute():
    packet=json.loads((PRIMARY/'source_packet.json').read_text())
    edges=[tuple(e) for e in packet['seams']]
    geometry=spatial.dual_boundaries(edges)
    chosen=list(schedules(edges))
    require(len(chosen)==132 and max(len(p) for _,_,p in chosen)==3,'schedule domain')
    rng=random.Random(104729)
    results={}; witnesses={}; native_events=native_holds=0
    for law,repeats in [('native',1),('half_current',1),('random_integer',3),('omitted_record',1)]:
        counts=dict(histories=0, continuity=0, closed_curves=0, degree_equals_current=0,
                    one_unit_every_step=0, regional_bound_every_step=0,
                    missing_record_against_loads=0, histories_with_omission=0)
        maximum_error=maximum_defect=0.
        for repeat in range(repeats):
            for origin,destination,path in chosen:
                for coins in itertools.product((0,1),repeat=len(path)):
                    initial,loads,tape,held,prefixes,events,omitted=run_history(origin,path,coins,edges,law,rng)
                    ordinary,defect=phase_readout(held,geometry)
                    maximum_defect=max(maximum_defect,defect)
                    period=divergence(held,edges)
                    q=[g-n for g,n in zip(initial,loads)]
                    counts['histories']+=1
                    counts['continuity']+=all(p['continuity'] for p in prefixes)
                    counts['one_unit_every_step']+=all(p['one_unit'] for p in prefixes)
                    counts['regional_bound_every_step']+=all(p['regional_bound'] for p in prefixes)
                    counts['histories_with_omission']+=omitted
                    counts['missing_record_against_loads']+=period!=q
                    if ordinary is not None:
                        deformed,_=phase_readout(held,geometry,True)
                        error=max(abs(float(w)-d) for w,d in zip(period,ordinary))
                        error=max(error,max(abs(float(w)-d) for w,d in zip(period,deformed)))
                        maximum_error=max(maximum_error,error)
                        counts['closed_curves']+=1
                        counts['degree_equals_current']+=error<1e-10
                    if law=='native':
                        native_events+=len(events)
                        native_holds+=sum(F(e['current'])==0 for e in events)
                    if law not in witnesses and (law=='native' or defect>1e-10 or period!=q
                            or not all(p['one_unit'] for p in prefixes)):
                        witnesses[law]={'origin':origin,'destination':destination,'schedule':list(path),
                            'events':events,'prefix_gates':prefixes,'terminal_gluing_defect':defect,
                            'held_period':[str(x) for x in period],'load_difference':[str(x) for x in q]}
        results[law]=dict(counts,maximum_phase_error=maximum_error,maximum_gluing_defect=maximum_defect)
    n,h,r,o=(results[k] for k in ('native','half_current','random_integer','omitted_record'))
    checks={
        'native_histories_enumerated_from_law':n['histories']==456 and native_events==888 and native_holds==564,
        'native_whole_unit_at_every_step':n['one_unit_every_step']==n['histories'],
        'native_all_regional_bounds_at_every_step':n['regional_bound_every_step']==n['histories'],
        'native_fixed_sample_contours':n['degree_equals_current']==n['histories'],
        'half_law_preserves_continuity_but_fails_closure':h['continuity']==h['histories'] and h['closed_curves']==0,
        'integer_fragmenting_law_passes_topology_but_fails_unit_gate':r['degree_equals_current']==r['histories'] and r['one_unit_every_step']<r['histories'],
        'integer_fragmenting_law_can_exceed_regional_domain':r['regional_bound_every_step']<r['histories'],
        'omission_detected_against_loads_for_every_affected_history':o['histories_with_omission']>0 and o['missing_record_against_loads']==o['histories_with_omission'],
        'damaged_integer_tape_still_has_valid_spatial_degree':o['degree_equals_current']==o['histories'],
    }
    require(all(checks.values()),'law mutation or positive control failed')
    pins={str(p.relative_to(HERE.parent)):sha(p) for p in (
        Path(__file__),PRIMARY/'native_kernel.py',PRIMARY/'native_charge_winding.py',
        PRIMARY/'source_packet.json',PRIMARY/'manifest.json')}
    return {'schema':'oph.row1.executed_law_controls.v1','development_only':True,
            'official_rows_changed':[],'samples_per_edge':SAMPLES_PER_EDGE,
            'declared_absolute_turn_bound':DECLARED_TURN_BOUND,'checks':checks,
            'arms':results,'witnesses':witnesses,'files_sha256':pins,
            'interpretation':'Topology diagnoses integer record classes; the separate native integer-unit law supplies nonfragmentation. Fixed phase samples test an observer encoding. The pointer reads counters, not these samples. Missing records are detected against source loads, not internally from the tape.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=HERE/'LAW_CONTROLS.json')
    args=parser.parse_args()
    result=execute()
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'arms':result['arms']},indent=2))
