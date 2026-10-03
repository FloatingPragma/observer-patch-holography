"""Independent exact verifier: ternary enumeration, Burnside and matchings.

No producer or native repair/geometry Python module is imported. The source
manifest is read directly; NetworkX supplies a separate graph-isomorphism
enumeration. Every settled quadrupole is checked in the exact field QQ(sqrt5).
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import product
import json
from pathlib import Path

import networkx as nx
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SPEC = HERE / 'native_repair_flavor_constraints_spec.json'
OUTPUT = HERE.parent / 'runs/flavor/native_repair_flavor_constraints.json'
REVIEWED_SPEC_SHA256 = 'a429cf26529ffc8750e8a12acbadc4d324c7fe3da7e44c7a977bf8f34e5b60ee'
PINS = (
    'code/a5_closure/manifests/echosahedral_federation_reference.json',
    'code/a5_closure/record_counting_mechanism_certificate.py',
    'code/a5_closure/a3_scheduler_kernel_certificate.py',
    'code/a5_closure/source_repair_generator_certificate.py',
    'code/particles/flavor/native_repair_flavor_constraints_spec.json',
    'code/particles/flavor/native_repair_flavor_constraints.py',
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def strict_load(path):
    def pairs(items):
        result = {}
        for key,value in items:
            require(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    def nonfinite(value):
        raise ValueError('nonfinite JSON value')
    return json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_constant=nonfinite)


def serialized(value):
    return json.dumps(value,sort_keys=True,separators=(',', ':'),allow_nan=False).encode()


def reviewed_spec():
    require(hashlib.sha256(SPEC.read_bytes()).hexdigest() == REVIEWED_SPEC_SHA256,
            'unreviewed specification')
    return strict_load(SPEC)


def maximal_matchings(edges):
    """Enumerate subsets by include/exclude; check maximality at leaves."""
    result = set()
    def visit(index, occupied, chosen):
        if index == len(edges):
            if all(i in occupied or j in occupied for i,j in edges):
                result.add(tuple(chosen))
            return
        i,j = edges[index]
        visit(index+1,occupied,chosen)
        if i not in occupied and j not in occupied:
            visit(index+1,occupied | {i,j},chosen+[(i,j)])
    visit(0,set(),[])
    return result


def global_menu(edges,group,vertices,distances):
    distance = np.array([[distances[i][j] for j in range(12)] for i in range(12)],dtype=np.int64)
    require(all(Counter(row)=={0:1,1:5,2:5,3:1} for row in distance),'global distance profile')
    require(np.all(distance.sum(axis=1)==18),'distance sum')
    for i in range(12):
        antipode = int(np.argmax(distance[i]))
        require(np.all(distance[i]+distance[antipode]==3),'antipodal distances')
    # Independent box enumeration, including the separately proved range-three
    # cones. The producer instead traverses signed graph-constrained heights.
    ternary = np.array(list(product(range(3),repeat=12)),dtype=np.int8)
    keep = np.min(ternary,axis=1)==0
    for i,j in edges:
        keep &= np.abs(ternary[:,i]-ternary[:,j])<=1
    cones = 3-distance
    heights = np.vstack([ternary[keep],cones]).astype(np.int64)
    require(len(heights)==6077 and np.all(cones.sum(axis=1)==18),'global class count')
    for row in cones:
        require(min(row)==0 and max(row)==3 and all(abs(row[i]-row[j])<=1 for i,j in edges),'range-three cone')
    states = sorted(tuple(map(int,row)) for row in heights)
    digest = hashlib.sha256(b''.join(bytes(z) for z in states)).hexdigest()
    fixed = sum(int(np.count_nonzero(np.all(heights==heights[:,p],axis=1))) for p in group)
    require(fixed%60 == 0,'global Burnside count')
    # All6077 spectral invariants are evaluated with integer matrix arithmetic:
    # 30Q=A+sqrt(5)B. Bounds on the small entries exclude int64 overflow.
    sqrt5 = sp.sqrt(5)
    phi = (1+sqrt5)/2
    aa,bb = [],[]
    for v in vertices:
        q = (30*(v*v.T/(phi+2)-sp.eye(3)/3)).applyfunc(sp.simplify)
        ar,br = [],[]
        for x in q:
            p = sp.Poly(x,sqrt5)
            a,b = p.coeff_monomial(1),p.coeff_monomial(sqrt5)
            require(a.is_Integer and b.is_Integer,'integer quadrupole coefficients')
            ar.append(int(a));br.append(int(b))
        aa.append(ar);bb.append(br)
    aa,bb = np.array(aa,dtype=np.int64),np.array(bb,dtype=np.int64)
    require(np.all(aa.sum(axis=0)==0) and np.all(bb.sum(axis=0)==0),'uniform-shift invariance')
    a = (heights@aa).reshape((-1,3,3))
    b = (heights@bb).reshape((-1,3,3))
    require(max(int(np.max(abs(a))),int(np.max(abs(b))))<1000,'integer arithmetic bound')
    t2 = np.einsum('nij,nji->n',a,a)+5*np.einsum('nij,nji->n',b,b)
    require(np.all(np.einsum('nij,nji->n',a,b)==0),'rational quadratic trace')
    t3a = np.einsum('nij,njk,nki->n',a,a,a)+15*np.einsum('nij,njk,nki->n',a,b,b)
    t3b = 3*np.einsum('nij,njk,nki->n',a,a,b)+5*np.einsum('nij,njk,nki->n',b,b,b)
    require(np.all(a[-12:]==0) and np.all(b[-12:]==0),'range-three zero quadrupole')
    types,totals = Counter(),Counter()
    residues = {str(i):Counter() for i in range(12)}
    for i,z in enumerate(heights):
        key = ';'.join(map(str,(Fraction(int(t2[i]),900),Fraction(int(t3a[i]),27000),Fraction(int(t3b[i]),27000))))
        types[key] += 1
        residues[str(int(sum(z))%12)][key] += 1
        totals[str(int(sum(z)))] += 1
    require(len(types)==31,'global spectral types')
    chis = Counter()
    zero_types = repeated_types = 0
    for key in types:
        t,a,b = map(sp.sympify,key.split(';'))
        if t == 0:
            zero_types += 1
            continue
        chi = sp.simplify((a+b*sqrt5)**2/t**3)
        chis[str(chi)] += 1
        repeated_types += int(chi==sp.Rational(1,6))
    gap = sp.Symbol('g',positive=True)
    chi_gap = (gap-1)**2*(gap+2)**2*(2*gap+1)**2/(24*(gap*gap+gap+1)**3)
    derivative = 9*gap*(gap-1)*(gap+1)*(gap+2)*(2*gap+1)/(8*(gap*gap+gap+1)**4)
    require(sp.simplify(sp.diff(chi_gap,gap)-derivative)==0,'log-gap monotonicity')
    require(sp.simplify(chi_gap.subs(gap,1/gap)-chi_gap)==0,'reciprocal gap symmetry')
    require(len(chis)==11 and zero_types==1 and repeated_types==8,'shape type counts')
    pairs = [(i,j) for i in range(12) for j in range(i+1,12) if distance[i,j]==3]
    for k,left in enumerate(pairs):
        for right in pairs[k+1:]:
            adjacency = [[int(tuple(sorted((i,j))) in edges) for j in right] for i in left]
            require(all(sum(row)==1 for row in adjacency) and
                    all(sum(row[j] for row in adjacency)==1 for j in range(2)),
                    'antipodal-pair seams are perfect matchings')
    sharp = [int(any(i in pair for pair in pairs[:3])) for i in range(12)]
    sharp_index = states.index(tuple(sharp))
    require(all(abs(sharp[i]-sharp[j])<=1 for i,j in edges),'sharp-bound settled witness')
    # Separate exact centered-pair calculation checks both the finite maximum
    # and the witness, complementing the source-free graph/variance proof.
    pair_heights = np.array([heights[:,i]+heights[:,j] for i,j in pairs]).T
    centered_six = 6*pair_heights-pair_heights.sum(axis=1)[:,None]
    require(np.all(np.ptp(pair_heights,axis=1)<=2),'pair-total range')
    require(np.all((centered_six**2).sum(axis=1)<=216),'centered variance bound')
    require(np.all(t2==20*(centered_six**2).sum(axis=1)), 'quadrupole pair normalization')
    require(max(t2)==4320,'sharp quadratic maximum')
    require(Fraction(sum((6*sum(sharp[i] for i in pair)-sum(sharp))**2 for pair in pairs),36)==6,
            'sharp variance attainment')
    require(tuple(sharp) in states and sharp_index>=0,'sharp witness in full census')
    x,m,M,mean = sp.symbols('x m M mean',real=True)
    require(sp.expand((x-m)*(M-x)+x*x-(m+M)*x+m*M)==0,'range variance identity')
    require(sp.expand((M-m)**2/4-(mean-m)*(M-mean)-(mean-(m+M)/2)**2)==0,
            'Popoviciu square identity')
    initial = [-1,13]+[0]*10
    endpoint = list(initial)
    schedule = []
    while True:
        chosen = next(((i,j) if endpoint[i]>endpoint[j] else (j,i)
                       for i,j in edges if abs(endpoint[i]-endpoint[j])>=2),None)
        if chosen is None:
            break
        i,j = chosen
        endpoint[i] -= 1
        endpoint[j] += 1
        schedule.append([i,j])
    require(sum(endpoint)==12 and set(endpoint)<={0,1,2},'signed endpoint bound')
    return {
        'distance_profile':[1,5,5,1],'sum_distances':18,
        'all_uniform_shift_classes':len(heights),'A5_orbits_modulo_uniform_shift':fixed//60,
        'ternary_minimum_zero_classes':len(heights)-12,'range_three_classes':12,
        'individual_quadrupole_spectral_types':len(types),'complete_height_sha256':digest,
        'digest_format':'Lexicographic normalized settled heights, minimum0, as12 bytes each.',
        'spectral_key':'tr(Q^2); rational_part(tr(Q^3)); sqrt(5)_coefficient(tr(Q^3))',
        'spectral_type_multiplicities':dict(types),'chi_spectral_type_counts':dict(chis),
        'zero_Q_spectral_types':zero_types,'repeated_eigenvalue_spectral_types':repeated_types,
        'nondegenerate_spectral_types':len(types)-zero_types-repeated_types,
        'nondegenerate_chi_values':len(chis)-1,'conditional_positive_log_gap_ratios':2*(len(chis)-2)+1,
        'chi_gap_derivative':str(sp.factor(derivative)),
        'sharp_quadrupole_bound':{'pair_matchings_checked':15,'pair_total_range_bound':2,
            'centered_pair_total_square_bound':6,'trace_Q2_upper_bound':'24/5',
            'attaining_state':sharp,'attaining_trace_Q2':'24/5'},
        'spectral_types_by_total_residue':{k:dict(v) for k,v in residues.items()},
        'normalized_total_histogram':dict(totals),
        'signed_outside_sector_witness':{'initial':initial,'accepted_schedule':schedule,'endpoint':endpoint},
    }


@lru_cache(None)
def independent_expectation():
    spec = reviewed_spec()
    carrier = strict_load(ROOT/spec['source']['carrier'])['carrier']
    labels = carrier['ports']
    edges = sorted(tuple(sorted((labels.index(i),labels.index(j)))) for i,j in carrier['edges'])
    graph = nx.Graph()
    graph.add_nodes_from(range(12))
    graph.add_edges_from(edges)
    oriented = [tuple(labels.index(p) for p in f) for f in carrier['oriented_faces']]
    positive = {f[i:]+f[:i] for f in oriented for i in range(3)}
    group = []
    automorphisms = 0
    for p in nx.algorithms.isomorphism.GraphMatcher(graph,graph).isomorphisms_iter():
        automorphisms += 1
        if all(tuple(p[i] for i in f) in positive for f in oriented):
            group.append(tuple(p[i] for i in range(12)))
    require(automorphisms == 120 and len(group) == 60, 'proper rotation group')

    def legal(z):
        return sorted((i,j) if z[i]>z[j] else (j,i)
                      for i,j in edges if abs(z[i]-z[j]) == 2)

    def child(z, move):
        i,j = move
        y = list(z)
        y[i] = y[j] = 1
        return tuple(y)

    states, settled = [], []
    digest = hashlib.sha256()
    spectrum = Counter()
    # A different enumeration from the producer's disjoint high/low sets.
    for z in product(range(3),repeat=12):
        if sum(z) != 12:
            continue
        states.append(z)
        moves = legal(z)
        children = sorted(child(z,e) for e in moves)
        digest.update(bytes(z)+bytes([len(moves)]))
        for y in children:
            require(sum(y) == 12 and sum(v*v for v in y) == sum(v*v for v in z)-2,
                    'conservation/strict descent')
            digest.update(bytes(y))
        spectrum[str(Fraction(30-len(moves),30))] += 1
        if not moves:
            settled.append(z)

    def image(z, p):
        return tuple(z[p.index(i)] for i in range(12))

    orbit_index, orbits = {}, []
    for z in sorted(states,key=lambda z:(z.count(2),z)):
        if z in orbit_index:
            continue
        members = {image(z,p) for p in group}
        index = len(orbits)
        for y in members:
            orbit_index[y] = index
        orbits.append((z,len(members)))
    # Independent Burnside generating functions check the orbit traversal.
    fixed = 0
    sink_fixed = 0
    for p in group:
        unused, lengths = set(range(12)), []
        while unused:
            i = min(unused)
            cycle = []
            while i not in cycle:
                cycle.append(i)
                unused.remove(i)
                i = p[i]
            lengths.append(len(cycle))
        polynomial = [1]+[0]*12
        for length in lengths:
            polynomial = [sum(polynomial[j-m*length] for m in range(3) if j>=m*length)
                          for j in range(13)]
        fixed += polynomial[12]
        sink_fixed += sum(all(z[i] == z[p[i]] for i in range(12)) for z in settled)
    require(fixed == 60*len(orbits) and sink_fixed == 60*9,'Burnside counts')

    cost_histogram = Counter()
    matching_digest = hashlib.sha256()
    matching_count = 0
    max_minimum = max_maximum = 0
    witness = None
    for z,size in orbits:
        moves = legal(z)
        matched = maximal_matchings(moves)
        matching_count += len(matched)
        matching_digest.update(serialized([z,sorted(matched)]))
        minimum, maximum = min(map(len,matched)), max(map(len,matched))
        max_minimum, max_maximum = max(max_minimum,minimum), max(max_maximum,maximum)
        cost_histogram[f'{minimum},{maximum}'] += size
        sinks = set()
        for matching in matched:
            endpoint = list(z)
            for i,j in matching:
                endpoint[i] = endpoint[j] = 1
            require(not legal(endpoint),'matching endpoint not settled')
            sinks.add(orbit_index[tuple(endpoint)])
        if witness is None and len(sinks)>1:
            # Direct ordered-word tree rather than the producer's orbit DP.
            weights = Counter()
            def visit(state,probability):
                next_moves = legal(state)
                if not next_moves:
                    weights[orbit_index[state]] += probability
                else:
                    for e in next_moves:
                        visit(child(state,e),probability/len(next_moves))
            visit(z,Fraction(1))
            witness = {'state':list(z),'terminal_orbits':[
                {'representative':list(orbits[j][0]),'probability':str(p)}
                for j,p in sorted(weights.items())]}

    sqrt5 = sp.sqrt(5)
    phi = (1+sqrt5)/2
    vertices = [sp.Matrix(v) for v in [(0,1,phi),(0,1,-phi),(0,-1,phi),(0,-1,-phi),
                (1,phi,0),(1,-phi,0),(-1,phi,0),(-1,-phi,0),
                (phi,0,1),(phi,0,-1),(-phi,0,1),(-phi,0,-1)]]
    require([(i,j) for i in range(12) for j in range(i+1,12)
             if sp.simplify((vertices[i]-vertices[j]).dot(vertices[i]-vertices[j])) == 4] == edges,
            'independent geometric frame')
    field = sp.QQ.algebraic_field(sqrt5)
    zero = field.zero
    basis = [[field.from_sympy(sp.simplify(x)) for x in v*v.T/(phi+2)] for v in vertices]
    invariants = {}
    types = Counter()
    for z in settled:
        q = [sum(((z[i]-1)*basis[i][j] for i in range(12)),zero) for j in range(9)]
        squared = [sum((q[3*i+k]*q[3*k+j] for k in range(3)),zero) for i in range(3) for j in range(3)]
        t2 = field.to_sympy(squared[0]+squared[4]+squared[8])
        t3 = field.to_sympy(sum((squared[3*i+k]*q[3*k+i] for i in range(3) for k in range(3)),zero))
        invariants[z] = t2,t3
        types[str(t2)+';'+str(t3)] += 1
    require(dict(types) == {'0;0':63,'8/5;0':180,'16/5;-24*sqrt(5)/25':30,
                            '16/5;24*sqrt(5)/25':30},'four exact quadrupole types')
    distances = dict(nx.all_pairs_shortest_path_length(graph))
    pairs = sorted((i,j) for i in range(12) for j in range(i+1,12) if distances[i][j] == 3)
    pair_index = {i:a for a,p in enumerate(pairs) for i in p}
    gram = [[Fraction(0) for _ in range(6)] for _ in range(6)]
    for i,j in edges:
        a,b = pair_index[i],pair_index[j]
        for u,v,sign in ((a,a,1),(b,b,1),(a,b,-1),(b,a,-1)):
            gram[u][v] += Fraction(sign,30)
    require(all(gram[i][j] == Fraction(2,5)*(int(i==j)-Fraction(1,6))
                for i in range(6) for j in range(6)), 'proposal Gram coefficient')
    for a,(i,_) in enumerate(pairs):
        for b,(j,_) in enumerate(pairs):
            inner = sum((basis[i][3*u+v]*basis[j][3*v+u] for u in range(3) for v in range(3)),zero)
            require(sp.simplify(field.to_sympy(inner)-sp.Rational(1,3)-sp.Rational(4,5)*(int(a==b)-sp.Rational(1,6))) == 0,
                    'quadrupole normalization')
    rows = []
    for z,size in orbits:
        if z not in invariants:
            continue
        t2,t3 = invariants[z]
        members = {image(z,p) for p in group}
        def cov(i,j):
            return sp.Rational(sum((y[i]-1)*(y[j]-1) for y in members),size)
        c = [cov(0,next(j for j in range(12) if distances[0][j] == d)) for d in range(4)]
        b3 = sp.simplify(c[0]+sqrt5*(c[1]-c[2])-c[3])
        b3p = sp.simplify(c[0]-sqrt5*(c[1]-c[2])-c[3])
        b5 = c[0]-c[1]-c[2]+c[3]
        require(all(cov(i,j) == c[distances[i][j]] for i in range(12) for j in range(12)), 'invariant covariance')
        require(sp.simplify(t2-8*b5) == 0, 'W5 covariance normalization')
        rows.append({'state':list(z),'orbit_size':size,'excess_missing_pairs':z.count(2),
                     'trace_Q2':str(t2),'trace_Q3':str(t3),'trace_Q4':str(sp.simplify(t2*t2/2)),
                     'chi':str(sp.simplify(t3*t3/t2**3)) if t2 else None,
                     'port_covariance_P3':str(b3),'port_covariance_P3prime':str(b3p),
                     'port_covariance_W5':str(b5),'quadrupole_covariance_scalar':str(t2/5)})
    # Derive the conditional gap-ratio set from the exact four spectra.
    allowed = [sp.sympify(x) for x in spec['conditional_affine_log_readout']['allowed_set']]
    ratios = []
    for expressions in spec['quadrupole']['stationary_spectra']:
        values = [sp.sympify(x) for x in expressions]
        require(sp.simplify(sum(values)) == 0, 'quadrupole trace')
        if all(x == 0 for x in values):
            continue
        require(all(sp.simplify(values[i+1]-values[i]).is_positive for i in (0,1)), 'ordered spectrum')
        ratio = sp.simplify((values[2]-values[1])/(values[1]-values[0]))
        ratios.append(ratio)
    require(len(ratios) == 3 and all(any(sp.simplify(a-r)==0 for r in ratios) for a in allowed),
            'conditional log-gap ratios')
    minimal = [z for z in states if z.count(2)<=1]
    minimal_sinks = [row for row in rows if row['excess_missing_pairs']<=1]
    global_result = global_menu(edges,group,vertices,distances)
    return {
        'schema':'oph.native_repair_flavor_constraints.v1','specification':spec,'scope':spec['scope'],
        'source_pins':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in PINS},
        'enumeration':{'states':len(states),'A5_orbits':len(orbits),'settled_states':len(settled),
                       'settled_A5_orbits':len(rows),'A5_invariant_stationary_simplex_dimension':len(rows)-1,
                       'complete_transition_sha256':digest.hexdigest(),
                       'digest_format':'Lexicographic states as12 bytes, active-seam count as1 byte, then lexicographic accepted successors as12 bytes each.'},
        'spectrum_K_multiplicities':dict(spectrum),
        'matching':{'maximal_matchings_checked':matching_count,'complete_matching_sha256':matching_digest.hexdigest(),
                    'digest_format':'Potential-then-lexicographic orbit representatives; concatenated compact JSON [state,sorted maximal directed-edge matchings].',
                    'minimum_maximum_move_histogram':dict(cost_histogram),
                    'largest_minimum_cost':max_minimum,'largest_maximum_steps':max_maximum},
        'minimal133':{'states':len(minimal),'settled_states':sum(x['orbit_size'] for x in minimal_sinks),
                      'settled_A5_orbits':len(minimal_sinks),'stationary_simplex_dimension':len(minimal_sinks)-1,
                      'spectrum_K_multiplicities':dict(Counter(str(Fraction(30-len(legal(z)),30)) for z in minimal)),
                      'stationary_orbit_sizes':sorted(x['orbit_size'] for x in minimal_sinks)},
        'nonunique_absorption_witness':witness,'settled_orbits':rows,
        'global_stationary_menu':global_result,
        'quadrupole_type_multiplicities':dict(types),
        'accepted_response':{'antipodal_pairs':[list(x) for x in pairs],
                             'full_pair_increment_Gram_coefficient':'2/5','quadrupole_map_squared_norm':'4/5',
                             'quadrupole_covariance_upper_coefficient':'8/25','stationary_accepted_moments':'zero'},
        'conditional_log_gap_ratios':spec['conditional_affine_log_readout']['allowed_set'],
    }


def verify_payload(payload):
    reviewed_spec()
    require(payload.get('source_pins') == {
        p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in PINS},
        'source pin changed since independent replay')
    expected = independent_expectation()
    require(serialized(payload) == serialized(expected), 'native repair receipt mismatch')
    return {'verified':True,'independent_states':73789,'independent_settled_states':303,
            'independent_A5_orbits':1275,'global_stationary_classes':6077,
            'individual_quadrupole_spectral_types':31,'conditional_positive_log_gap_ratios':19,
            'physical_mass_identification':False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt',nargs='?',type=Path,default=OUTPUT)
    args = parser.parse_args()
    print(json.dumps(verify_payload(strict_load(args.receipt)),indent=2))


if __name__ == '__main__':
    main()
