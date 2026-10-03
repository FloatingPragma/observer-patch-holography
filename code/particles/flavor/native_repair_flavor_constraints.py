"""Exact native repair constraints on a finite one-record-deviation sector.

The integer source relation supplies the dynamics. Its composition with the
declared uniform scheduler is an explicitly conditional count-level operator,
not a selected physical clock, vacuum, Yukawa map, or quark-mass prediction.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import combinations
import json
from pathlib import Path
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SPEC = HERE / 'native_repair_flavor_constraints_spec.json'
OUTPUT = HERE.parent / 'runs/flavor/native_repair_flavor_constraints.json'
sys.path.insert(0, str(ROOT / 'code/a5_closure'))
import record_counting_mechanism_certificate as native
import source_repair_generator_certificate as geometry

PINS = (
    'code/a5_closure/manifests/echosahedral_federation_reference.json',
    'code/a5_closure/record_counting_mechanism_certificate.py',
    'code/a5_closure/a3_scheduler_kernel_certificate.py',
    'code/a5_closure/source_repair_generator_certificate.py',
    'code/particles/flavor/native_repair_flavor_constraints_spec.json',
    'code/particles/flavor/native_repair_flavor_constraints.py',
)


def canonical_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode()


def _states():
    states = []
    for m in range(7):
        for high in combinations(range(12), m):
            remaining = sorted(set(range(12)) - set(high))
            for low in combinations(remaining, m):
                z = [1]*12
                for i in high:
                    z[i] = 2
                for i in low:
                    z[i] = 0
                states.append(tuple(z))
    return states


def _orbit(z, group):
    result = set()
    for g in group:
        y = [0]*12
        for i, j in enumerate(g):
            y[j] = z[i]
        result.add(tuple(y))
    return result


@lru_cache(None)
def _maximal_matchings(edges):
    if not edges:
        return frozenset({()})
    result = set()
    for e in edges:
        residual = tuple(f for f in edges if not set(e) & set(f))
        for matching in _maximal_matchings(residual):
            result.add(tuple(sorted((e,)+matching)))
    return frozenset(result)


def _quadrupole_rows(orbits, group, edges):
    phi = (1+sp.sqrt(5))/2
    vertices = sp.Matrix([(0,1,phi),(0,1,-phi),(0,-1,phi),(0,-1,-phi),
                          (1,phi,0),(1,-phi,0),(-1,phi,0),(-1,-phi,0),
                          (phi,0,1),(phi,0,-1),(-phi,0,1),(-phi,0,-1)])
    assert (vertices.T*vertices-4*(phi+2)*sp.eye(3)).applyfunc(sp.simplify) == sp.zeros(3)
    derived = [(i,j) for i in range(12) for j in range(i+1,12)
               if sp.simplify((vertices[i,:]-vertices[j,:]).dot(vertices[i,:]-vertices[j,:])) == 4]
    assert derived == edges
    antipodes = [next(j for j in range(12) if vertices[j,:] == -vertices[i,:]) for i in range(12)]
    pairs = sorted({tuple(sorted((i,j))) for i,j in enumerate(antipodes)})
    pair_index = {i:a for a,pair in enumerate(pairs) for i in pair}
    proposal_gram = sp.zeros(6)
    for i,j in edges:
        step = sp.eye(6)[:,pair_index[i]]-sp.eye(6)[:,pair_index[j]]
        proposal_gram += step*step.T/30
    assert proposal_gram == sp.Rational(2,5)*(sp.eye(6)-sp.ones(6)/6)
    projectors = [vertices[i,:].T*vertices[i,:]/(phi+2)-sp.eye(3)/3 for i,j in pairs]
    gram = sp.Matrix(6,6,lambda i,j:sp.simplify(sp.trace(projectors[i]*projectors[j])))
    assert gram == sp.Rational(4,5)*(sp.eye(6)-sp.ones(6)/6)
    rows = []
    for orbit in orbits:
        if orbit['active']:
            continue
        delta = sp.Matrix([x-1 for x in orbit['state']])
        q = (vertices.T*sp.diag(*delta)*vertices/(phi+2)).applyfunc(sp.simplify)
        t2, t3 = (sp.simplify(sp.trace(q**k)) for k in (2,3))
        covariance = sp.zeros(12)
        for y in _orbit(orbit['state'], group):
            d = sp.Matrix([x-1 for x in y])
            covariance += d*d.T/orbit['size']
        c0 = covariance[0,0]
        c1 = covariance[0,next(j for i,j in edges if i == 0)]
        c3 = covariance[0,antipodes[0]]
        neighbors = {j for i,j in edges if i == 0} | {i for i,j in edges if j == 0}
        c2 = covariance[0,next(j for j in range(12) if j != 0 and j not in neighbors and j != antipodes[0])]
        b3 = sp.simplify(c0+sp.sqrt(5)*(c1-c2)-c3)
        b3p = sp.simplify(c0-sp.sqrt(5)*(c1-c2)-c3)
        b5 = c0-c1-c2+c3
        assert sp.simplify(3*b3+3*b3p+5*b5-delta.dot(delta)) == 0
        assert sp.simplify(t2-8*b5) == 0
        rows.append({'state':list(orbit['state']),'orbit_size':orbit['size'],
                     'excess_missing_pairs':orbit['state'].count(2),
                     'trace_Q2':str(t2),'trace_Q3':str(t3),
                     'trace_Q4':str(sp.simplify(t2*t2/2)),
                     'chi':str(sp.simplify(t3*t3/t2**3)) if t2 else None,
                     'port_covariance_P3':str(b3),'port_covariance_P3prime':str(b3p),
                     'port_covariance_W5':str(b5),
                     'quadrupole_covariance_scalar':str(t2/5)})
    return rows, pairs, vertices


def _global_stationary_menu(edges, group, vertices):
    """All integer one-Lipschitz heights, modulo their common shift.

    Fixing the height of port zero to zero removes the translation freedom.
    Constrained graph traversal is independent of a ternary-box enumeration.
    """
    adjacency = [set() for _ in range(12)]
    for i,j in edges:
        adjacency[i].add(j)
        adjacency[j].add(i)
    distances = []
    for start in range(12):
        d = {start:0}
        queue = [start]
        for i in queue:
            for j in adjacency[i]:
                if j not in d:
                    d[j] = d[i]+1
                    queue.append(j)
        distances.append([d[i] for i in range(12)])
    assert all(Counter(d) == {0:1,1:5,2:5,3:1} and sum(d) == 18 for d in distances)
    for i in range(12):
        antipode = distances[i].index(3)
        assert all(distances[i][j]+distances[antipode][j] == 3 for j in range(12))
    order = sorted(range(12),key=lambda i:(distances[0][i],i))
    height = {0:0}
    states = set()
    def visit(index):
        if index == 12:
            low = min(height.values())
            states.add(tuple(height[i]-low for i in range(12)))
            return
        i = order[index]
        known = [height[j] for j in adjacency[i] if j in height]
        lo = max([-distances[0][i]]+[v-1 for v in known])
        hi = min([distances[0][i]]+[v+1 for v in known])
        for value in range(lo,hi+1):
            height[i] = value
            visit(index+1)
        height.pop(i,None)
    visit(1)
    range_three = {tuple(3-d for d in row) for row in distances}
    assert {z for z in states if max(z)==3} == range_three
    digest = hashlib.sha256(b''.join(bytes(z) for z in sorted(states))).hexdigest()
    by_residue = {str(i):Counter() for i in range(12)}
    types = Counter()
    totals = Counter(str(sum(z)) for z in states)
    seen = set()
    orbit_count = 0
    phi = (1+sp.sqrt(5))/2
    for z in sorted(states):
        if z in seen:
            continue
        members = _orbit(z,group)
        assert members <= states
        seen.update(members)
        orbit_count += 1
        q = (vertices.T*sp.diag(*z)*vertices/(phi+2)-sp.Rational(sum(z),3)*sp.eye(3)).applyfunc(sp.simplify)
        t2 = sp.simplify(sp.trace(q*q))
        t3 = sp.expand(sp.simplify(sp.trace(q**3)))
        b = t3.coeff(sp.sqrt(5))
        a = sp.simplify(t3-b*sp.sqrt(5))
        assert t2.is_Rational and a.is_Rational and b.is_Rational
        key = ';'.join(map(str,(t2,a,b)))
        types[key] += len(members)
        by_residue[str(sum(z)%12)][key] += len(members)
    assert len(states) == 6077 and len(types) == 31
    chi_types = Counter()
    zero_types = repeated_types = 0
    for key in types:
        t2,a,b = map(sp.sympify,key.split(';'))
        if t2 == 0:
            zero_types += 1
            continue
        chi = sp.simplify((a+b*sp.sqrt(5))**2/t2**3)
        chi_types[str(chi)] += 1
        repeated_types += int(chi == sp.Rational(1,6))
    gap = sp.Symbol('g',positive=True)
    q = [-(2+gap)/3,(1-gap)/3,(1+2*gap)/3]
    chi_gap = (gap-1)**2*(gap+2)**2*(2*gap+1)**2/(24*(gap*gap+gap+1)**3)
    assert sp.simplify(sum(x**3 for x in q)**2/sum(x*x for x in q)**3-chi_gap) == 0
    assert sp.simplify(chi_gap.subs(gap,1/gap)-chi_gap) == 0
    derivative = sp.factor(sp.diff(chi_gap,gap))
    assert sp.simplify(derivative/(gap-1)).is_positive
    assert sp.limit(chi_gap,gap,sp.oo) == sp.Rational(1,6)
    assert len(chi_types) == 11 and zero_types == 1 and repeated_types == 8
    pairs = sorted({tuple(sorted((i,distances[i].index(3)))) for i in range(12)})
    for ai, left in enumerate(pairs):
        for right in pairs[ai+1:]:
            between = [(i,j) for i,j in edges if (i in left and j in right) or (j in left and i in right)]
            assert len(between) == 2
            assert Counter(v for edge in between for v in edge) == Counter(left+right)
    # Popoviciu's bound follows by averaging (b-m)(M-b)>=0:
    # variance <= (mean-m)(M-mean) <= (M-m)^2/4 <= 1.
    x,m,M,mean = sp.symbols('x m M mean',real=True)
    assert sp.expand((x-m)*(M-x)+x*x-(m+M)*x+m*M) == 0
    assert sp.expand((M-m)**2/4-(mean-m)*(M-mean)-(mean-(m+M)/2)**2) == 0
    sharp = tuple(int(i in set(sum((list(p) for p in pairs[:3]),[]))) for i in range(12))
    assert all(abs(sharp[i]-sharp[j]) <= 1 for i,j in edges)
    sharp_q = (vertices.T*sp.diag(*sharp)*vertices/(phi+2)-sp.Rational(sum(sharp),3)*sp.eye(3)).applyfunc(sp.simplify)
    assert sp.simplify(sp.trace(sharp_q**2)) == sp.Rational(24,5)
    assert max(sp.sympify(key.split(';')[0]) for key in types) == sp.Rational(24,5)
    # A source-admitted signed preparation outside the transient ternary box.
    initial = (-1,13)+(0,)*10
    endpoint = initial
    schedule = []
    while legal := native.admissible_moves(endpoint,edges):
        move = legal[0]
        endpoint = native.apply_move(endpoint,move)
        schedule.append(list(move))
    assert sum(endpoint) == 12 and set(endpoint) <= {0,1,2}
    return {
        'distance_profile':[1,5,5,1],'sum_distances':18,
        'all_uniform_shift_classes':len(states),'A5_orbits_modulo_uniform_shift':orbit_count,
        'ternary_minimum_zero_classes':len(states)-len(range_three),'range_three_classes':len(range_three),
        'individual_quadrupole_spectral_types':len(types),
        'complete_height_sha256':digest,
        'digest_format':'Lexicographic normalized settled heights, minimum0, as12 bytes each.',
        'spectral_key':'tr(Q^2); rational_part(tr(Q^3)); sqrt(5)_coefficient(tr(Q^3))',
        'spectral_type_multiplicities':dict(types),
        'chi_spectral_type_counts':dict(chi_types),'zero_Q_spectral_types':zero_types,
        'repeated_eigenvalue_spectral_types':repeated_types,
        'nondegenerate_spectral_types':len(types)-zero_types-repeated_types,
        'nondegenerate_chi_values':len(chi_types)-1,
        'conditional_positive_log_gap_ratios':2*(len(chi_types)-2)+1,
        'chi_gap_derivative':str(derivative),
        'sharp_quadrupole_bound':{'pair_matchings_checked':15,'pair_total_range_bound':2,
            'centered_pair_total_square_bound':6,'trace_Q2_upper_bound':'24/5',
            'attaining_state':list(sharp),'attaining_trace_Q2':'24/5'},
        'spectral_types_by_total_residue':{k:dict(v) for k,v in by_residue.items()},
        'normalized_total_histogram':dict(totals),
        'signed_outside_sector_witness':{'initial':list(initial),'accepted_schedule':schedule,'endpoint':list(endpoint)},
    }


def build_payload():
    spec = json.loads(SPEC.read_text())
    g = geometry.classify()
    edges, group = g['edges'], g['group']
    states = _states()
    state_orbit, orbits = {}, []
    for z in sorted(states, key=lambda z:(sum(x*x for x in z),z)):
        if z in state_orbit:
            continue
        members = _orbit(z, group)
        index = len(orbits)
        for y in members:
            state_orbit[y] = index
        orbits.append({'state':z,'size':len(members)})
    digest = hashlib.sha256()
    spectrum = Counter()
    for z in sorted(states):
        legal = native.admissible_moves(z, edges)
        children = sorted(native.apply_move(z,e) for e in legal)
        for y in children:
            assert sum(y) == 12 and all(x in (0,1,2) for x in y)
            assert sum(x*x for x in y) == sum(x*x for x in z)-2
        digest.update(bytes(z)+bytes([len(legal)]))
        for y in children:
            digest.update(bytes(y))
        spectrum[str(Fraction(30-len(legal),30))] += 1
    absorption, minimum, maximum = [], [], []
    matching_digest = hashlib.sha256()
    matching_count = 0
    cost_histogram = Counter()
    branch_witness = None
    for index, orbit in enumerate(orbits):
        z = orbit['state']
        legal = native.admissible_moves(z, edges)
        orbit['active'] = len(legal)
        successors = Counter(state_orbit[native.apply_move(z,e)] for e in legal)
        if not legal:
            absorption.append({index:Fraction(1)})
            minimum.append(0)
            maximum.append(0)
        else:
            assert all(j < index for j in successors)
            weights = Counter()
            for j, count in successors.items():
                for sink, probability in absorption[j].items():
                    weights[sink] += Fraction(count,len(legal))*probability
            assert sum(weights.values()) == 1
            absorption.append(dict(weights))
            minimum.append(1+min(minimum[j] for j in successors))
            maximum.append(1+max(maximum[j] for j in successors))
            if branch_witness is None and len(weights) > 1:
                branch_witness = {'state':list(z),'terminal_orbits':[
                    {'representative':list(orbits[j]['state']),'probability':str(probability)}
                    for j,probability in sorted(weights.items())]}
        matched = _maximal_matchings(tuple(sorted(legal)))
        matching_count += len(matched)
        assert min(map(len,matched)) == minimum[index]
        assert max(map(len,matched)) == maximum[index]
        terminal_orbits = set()
        for matching in matched:
            y = list(z)
            for i,j in matching:
                y[i] = y[j] = 1
            assert not native.admissible_moves(y,edges)
            terminal_orbits.add(state_orbit[tuple(y)])
        assert terminal_orbits == set(absorption[index])
        matching_digest.update(canonical_bytes([z,sorted(matched)]))
        cost_histogram[f'{minimum[index]},{maximum[index]}'] += orbit['size']
    sinks = [o for o in orbits if not o['active']]
    minimal = [o for o in orbits if o['state'].count(2) <= 1]
    minimal_sinks = [o for o in minimal if not o['active']]
    minimal_spectrum = Counter()
    for orbit in minimal:
        minimal_spectrum[str(Fraction(30-orbit['active'],30))] += orbit['size']
    rows, pairs, vertices = _quadrupole_rows(orbits,group,edges)
    global_menu = _global_stationary_menu(edges,group,vertices)
    types = Counter()
    for row in rows:
        types[row['trace_Q2']+';'+row['trace_Q3']] += row['orbit_size']
    return {
        'schema':'oph.native_repair_flavor_constraints.v1',
        'specification':spec,'scope':spec['scope'],
        'source_pins':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in PINS},
        'enumeration':{'states':len(states),'A5_orbits':len(orbits),
                       'settled_states':sum(o['size'] for o in sinks),'settled_A5_orbits':len(sinks),
                       'A5_invariant_stationary_simplex_dimension':len(sinks)-1,
                       'complete_transition_sha256':digest.hexdigest(),
                       'digest_format':'Lexicographic states as12 bytes, active-seam count as1 byte, then lexicographic accepted successors as12 bytes each.'},
        'spectrum_K_multiplicities':dict(spectrum),
        'matching':{'maximal_matchings_checked':matching_count,
                    'complete_matching_sha256':matching_digest.hexdigest(),
                    'digest_format':'Potential-then-lexicographic orbit representatives; concatenated compact JSON [state,sorted maximal directed-edge matchings].',
                    'minimum_maximum_move_histogram':dict(cost_histogram),
                    'largest_minimum_cost':max(minimum),'largest_maximum_steps':max(maximum)},
        'minimal133':{'states':sum(o['size'] for o in minimal),
                      'settled_states':sum(o['size'] for o in minimal_sinks),
                      'settled_A5_orbits':len(minimal_sinks),'stationary_simplex_dimension':len(minimal_sinks)-1,
                      'spectrum_K_multiplicities':dict(minimal_spectrum),
                      'stationary_orbit_sizes':sorted(o['size'] for o in minimal_sinks)},
        'nonunique_absorption_witness':branch_witness,
        'settled_orbits':rows,'quadrupole_type_multiplicities':dict(types),
        'global_stationary_menu':global_menu,
        'accepted_response':{'antipodal_pairs':[list(x) for x in pairs],
                             'full_pair_increment_Gram_coefficient':'2/5',
                             'quadrupole_map_squared_norm':'4/5',
                             'quadrupole_covariance_upper_coefficient':'8/25',
                             'stationary_accepted_moments':'zero'},
        'conditional_log_gap_ratios':spec['conditional_affine_log_readout']['allowed_set'],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    rendered = json.dumps(build_payload(),indent=2,sort_keys=True)+'\n'
    if args.check:
        if OUTPUT.read_text() != rendered:
            raise SystemExit('Native repair receipt differs; regenerate')
        print('Native repair flavor constraints verified')
    else:
        OUTPUT.write_text(rendered)
        print(OUTPUT)


if __name__ == '__main__':
    main()
