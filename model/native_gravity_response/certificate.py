"""Exact regional response of the native closed-history quantum state.

The finite source law and its count-derived density are inputs from the
native reconstruction. Regional restriction, exterior multiplicity and
midpoint-filter state variations are calculated on that same carrier.
No geometric area, physical stress or Einstein balance is defined here.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from math import isqrt, lcm
from pathlib import Path

import numpy as np
import sympy as sy

from tools import strict_json
from ..native_quantum_reconstruction import certificate as source
from ..native_thermodynamic_bridge.certificate import combine, log_rational

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
STATES = source.STATES
DIM = len(STATES)
SEAM = source.quantum_mixing.seams()[0]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def trace(matrix):
    return sum((matrix[i, i] for i in range(len(matrix))), Q())


def hs_squared(matrix):
    return sum((x*x for x in matrix.flat), Q())


@lru_cache(maxsize=1)
def source_data():
    counts, denominator = source.source_counts()
    C = np.array(counts, dtype=object)
    square = C@C
    Z = int(trace(square))
    sigma = np.array([[Q(x, Z) for x in row] for row in square], dtype=object)
    require(tuple(map(tuple, sigma)) == source.reconstructed_density(2),
            'regional state differs from the native count reconstruction')
    # The source Gershgorin bound C >= 40 I is independent of diagonalization.
    lower = min(2*int(C[i, i])-denominator for i in range(DIM))
    require(lower == 40 and Z == 171720, 'source spectral data changed')
    return C, Z, sigma, Q(lower**2, Z)


@lru_cache(maxsize=64)
def layout(region):
    require(type(region) is tuple and 0 < len(region) < 12
            and tuple(sorted(set(region))) == region
            and all(type(i) is int and 0 <= i < 12 for i in region),
            'a sorted proper nonempty native port region is required')
    groups = {}
    for i, state in enumerate(STATES):
        inside = tuple(p for p in region if state[p])
        outside = tuple(p for p in range(12) if p not in region and state[p])
        groups.setdefault(len(inside), {})[inside, outside] = i
    result = {}
    for q, lookup in sorted(groups.items()):
        local = tuple(sorted({a for a, _ in lookup}))
        exterior = tuple(sorted({b for _, b in lookup}))
        require(len(lookup) == len(local)*len(exterior), 'sector does not factor')
        result[q] = local, exterior, lookup
    return result


def restrict(matrix, region):
    """Intrinsic regional density: direct sum of unnormalized partial traces."""
    result = {}
    for q, (local, exterior, lookup) in layout(region).items():
        result[q] = np.array([[sum((matrix[lookup[a,b], lookup[c,b]]
                                  for b in exterior), Q()) for c in local]
                              for a in local], dtype=object)
    return result


def complete(matrix, region):
    """Tracial conditional expectation, with counted maximally mixed exterior."""
    result = np.full((DIM, DIM), Q(), dtype=object)
    blocks = restrict(matrix, region)
    for q, (local, exterior, lookup) in layout(region).items():
        for i, a in enumerate(local):
            for j, c in enumerate(local):
                for b in exterior:
                    result[lookup[a,b], lookup[c,b]] = blocks[q][i,j]/len(exterior)
    return result


def filter_tangent(values):
    """C diag(f) C / Z from weighted two-attempt native paths."""
    require(type(values) is tuple and len(values) == DIM
            and all(type(x) in (int, Q) for x in values),
            'one exact rational filter value per native configuration is required')
    values = tuple(map(Q, values))
    require(max(map(abs, values)) <= 1, 'filter exceeds the unit interval bound')
    C, Z, _, _ = source_data()
    denominator = lcm(*(x.denominator for x in values))
    numerator = np.array([int(x*denominator) for x in values], dtype=object)
    product = (C*numerator)@C
    tangent = np.array([[Q(x, denominator*Z) for x in row] for row in product], dtype=object)
    require(trace(tangent) == 0, 'filter changes normalization of the return event')
    return tangent


def filtered_state(values, parameter):
    require(type(parameter) in (int, Q) and abs(parameter) < 1,
            'an exact filter parameter strictly between minus one and one is required')
    return source_data()[2] + Q(parameter)*filter_tangent(values)


def complete_return_weights(values, parameter):
    """Retain accepted returns, rejected returns and every nonreturn event."""
    filtered_state(values,parameter)
    C,_,_,_ = source_data()
    accepted,rejected = Q(),Q()
    for initial in range(DIM):
        for middle in range(DIM):
            mass = Q(int(C[initial,middle]*C[middle,initial]),DIM*60**2)
            acceptance = (1+parameter*values[middle])/2
            accepted += mass*acceptance
            rejected += mass*(1-acceptance)
    nonreturn = 1-source.return_probability(2)
    require(accepted+rejected+nonreturn == 1,'return filter drops complete outcomes')
    require(accepted == rejected == source.return_probability(2)/2,
            'centered filter changed the return normalization')
    return {'accepted_return':str(accepted),'rejected_return':str(rejected),
            'nonreturn':str(nonreturn)}


def rational_sqrt(value):
    value = Q(value)
    require(value >= 0, 'negative spectral discriminant')
    a, b = isqrt(value.numerator), isqrt(value.denominator)
    require(a*a == value.numerator and b*b == value.denominator,
            'this exact spectral receipt requires rational eigenvalues')
    return Q(a,b)


def spectrum(block):
    """Exact eigenvalues and projectors of the one/two-dimensional seam blocks."""
    if len(block) == 1:
        require(block[0,0] > 0, 'nonfaithful regional state')
        return ((block[0,0], np.array([[Q(1)]], dtype=object)),)
    require(block.shape == (2,2) and block[0,1] == block[1,0],
            'exact spectrum requires a real symmetric seam block')
    a, b, d = block[0,0], block[0,1], block[1,1]
    gap = rational_sqrt((a-d)**2+4*b*b)
    lo, hi = (a+d-gap)/2, (a+d+gap)/2
    require(lo > 0, 'nonfaithful regional state')
    eye = np.array([[Q(1),Q()],[Q(),Q(1)]],dtype=object)
    if gap == 0:
        return ((lo,eye),)
    P = (block-lo*eye)/gap
    require(np.array_equal(P@P,P), 'spectral projector is not idempotent')
    return ((lo,eye-P),(hi,P))


def bulk_entropy(blocks):
    return combine(*[(-value*trace(projector),log_rational(value))
                     for block in blocks.values() for value,projector in spectrum(block)])


def edge_entropy(blocks, region):
    return combine(*[(trace(block),log_rational(Q(len(layout(region)[q][1]))))
                     for q,block in blocks.items()])


def completion_entropy(blocks, region):
    # Independently use each full completed eigenvalue and its multiplicity.
    return combine(*[(-value*trace(projector),log_rational(value/len(layout(region)[q][1])))
                     for q,block in blocks.items() for value,projector in spectrum(block)])


def bulk_sector_split(blocks):
    shannon = combine(*[(-trace(B),log_rational(trace(B))) for B in blocks.values()])
    conditional = combine(*[(-value*trace(P),log_rational(value/trace(B)))
                            for B in blocks.values() for value,P in spectrum(B)])
    return combine((1,shannon),(1,conditional))


def modular_pairing(variation, reference):
    return combine(*[(-trace(variation[q]@projector),log_rational(value))
                     for q,block in reference.items() for value,projector in spectrum(block)])


def relative_entropy(state, reference):
    return combine((1,modular_pairing(state,reference)),(-1,bulk_entropy(state)))


@lru_cache(maxsize=2048)
def log_interval(value, terms=40):
    """Exact atanh series after binary range reduction; complete tail included."""
    value = Q(value)
    require(value > 0, 'positive logarithm argument required')
    k = value.numerator.bit_length()-value.denominator.bit_length()
    r = value/Q(2)**k
    if r < 1:
        r *= 2
        k -= 1
    def local(q):
        x = (q-1)/(q+1)
        lo = 2*sum((x**(2*j+1)/Q(2*j+1) for j in range(terms)),Q())
        return lo,lo+2*x**(2*terms+1)/((2*terms+1)*(1-x*x))
    a,b = local(r)
    c,d = sorted(k*t for t in local(Q(2)))
    return a+c,b+d


def expression_interval(expression):
    lo,hi = Q(),Q()
    for base,coefficient in expression:
        a,b = sorted(coefficient*x for x in log_interval(base))
        lo += a
        hi += b
    return lo,hi


def expression_receipt(expression):
    lo,hi = expression_interval(expression)
    scale = 10**12
    lower = Q((lo*scale).__floor__(),scale)
    upper = Q((hi*scale).__ceil__(),scale)
    return {'prime_log_coefficients':{str(p):str(c) for p,c in expression},
            'outward_interval':[str(lower),str(upper)]}


@lru_cache(maxsize=1)
def occupation_fixed_filter():
    """Project the native pair occupation onto the exact moment-null space."""
    C,_,_,_ = source_data()
    squared = C*C  # Entrywise square, for diagonal response moments.
    rows = [[int(sum(squared[:,x])) for x in range(DIM)]]
    rows += [[int(sum(STATES[y][p]*squared[y,x] for y in range(DIM)))
              for x in range(DIM)] for p in range(12)]
    matrix = sy.Matrix(rows)
    independent_rows = matrix.T.rref()[1]
    require(len(independent_rows) == 12, 'occupation response rank changed')
    basis = matrix[list(independent_rows),:]
    pair = sy.Matrix([int(all(state[p] for p in SEAM)) for state in STATES])
    vector = pair-basis.T*(basis*basis.T).inv()*basis*pair
    require(matrix*vector == sy.zeros(13,1), 'pair projection violates a fixed moment')
    largest = max(abs(v) for v in vector)
    require(largest > 0, 'pair occupation lies in the constrained response space')
    values = tuple(Q(v/largest) for v in vector)
    return values,DIM-len(independent_rows)


def response_receipt(values, parameter):
    _,_,sigma,floor = source_data()
    tangent = filter_tangent(values)
    reference = restrict(sigma,SEAM)
    delta = restrict(tangent,SEAM)
    state = restrict(filtered_state(values,parameter),SEAM)
    slope = modular_pairing(delta,reference)
    entropy_change = combine((1,bulk_entropy(state)),(-1,bulk_entropy(reference)))
    relative = relative_entropy(state,reference)
    require(combine((1,entropy_change),(-parameter,slope),(1,relative)) == (),
            'finite modular first-law identity failed')
    lower,upper = expression_interval(relative)
    bound = parameter**2*sum(hs_squared(B) for B in delta.values())/(2*(1-abs(parameter))*floor)
    require(lower > 0 and upper <= bound, 'certified entropy remainder violated')
    edge_change = combine((1,edge_entropy(state,SEAM)),(-1,edge_entropy(reference,SEAM)))
    full_change = combine((1,completion_entropy(state,SEAM)),(-1,completion_entropy(reference,SEAM)))
    require(combine((1,full_change),(-1,entropy_change),(-1,edge_change)) == (),
            'same-state bulk/edge response split failed')
    return {'parameter':str(parameter), 'filter':[str(x) for x in values],
            'complete_terminal_weights':complete_return_weights(values,parameter),
            'bulk_entropy_derivative':expression_receipt(slope),
            'bulk_entropy_change':expression_receipt(entropy_change),
            'edge_entropy_change':expression_receipt(edge_change),
            'relative_entropy':expression_receipt(relative),
            'quadratic_remainder_upper':str(bound),
            'normalized_trace_tangent':str(trace(tangent)),
            'regional_commutator_hs_squared':str(sum(hs_squared(reference[q]@B-B@reference[q])
                                                     for q,B in delta.items())),
            'complete_spectrum':[str(value) for B in state.values() for value,P in spectrum(B)]}


def report():
    C,Z,sigma,floor = source_data()
    reference = restrict(sigma,SEAM)
    regions = []
    for region in source.quantum_mixing.seams():
        blocks = restrict(sigma,region)
        bulk,edge,total = bulk_entropy(blocks),edge_entropy(blocks,region),completion_entropy(blocks,region)
        require(bulk == bulk_sector_split(blocks), 'central Shannon term is misplaced')
        require(combine((1,total),(-1,bulk),(-1,edge)) == (), 'entropy split failed')
        completion = complete(sigma,region)
        require(trace(completion) == 1 and np.array_equal(complete(completion,region),completion),
                'regional map is not a trace-preserving projection')
        charge = np.array([sum(state[p] for p in region) for state in STATES],dtype=object)
        obstruction = hs_squared(sigma*charge-charge[:,None]*sigma)
        require(obstruction > 0, 'source center-mixing obstruction disappeared')
        regions.append({'region':list(region),'sector_dimensions':
                        {str(q):[len(a),len(b)] for q,(a,b,_) in layout(region).items()},
                        'sector_weights':{str(q):str(trace(B)) for q,B in blocks.items()},
                        'intrinsic_spectrum':[str(v) for B in blocks.values() for v,P in spectrum(B)],
                        'completion_defect_hs_squared':str(hs_squared(sigma-completion)),
                        'center_commutator_hs_squared':str(obstruction),
                        'bulk_entropy':expression_receipt(bulk),'edge_entropy':expression_receipt(edge)})
    small,large = (SEAM[0],),SEAM
    require(np.array_equal(complete(complete(sigma,large),small),complete(sigma,small)),
            'nested regional expectation tower failed')
    values,nullity = occupation_fixed_filter()
    tangent = filter_tangent(values)
    moments = [sum(STATES[i][p]*tangent[i,i] for i in range(DIM)) for p in range(12)]
    require(all(x == 0 for x in moments), 'occupation-fixed response changes a port mean')
    fixed = response_receipt(values,Q(1,2))
    require(expression_interval(modular_pairing(restrict(tangent,SEAM),reference))[0] > 0,
            'static-cut nonstationarity witness has no positive bulk slope')
    fixed.update({'all_twelve_occupation_derivatives':[str(x) for x in moments],
                  'filter_nullity':nullity,'fixed_bare_cut_edges':sum((a in SEAM)!=(b in SEAM)
                                       for a,b in source.quantum_mixing.seams()),
                  'bare_cut_area_derivative':'0','fixed_proper_volume_claimed':False})
    contrast = tuple(state[SEAM[0]]-state[SEAM[1]] for state in STATES)
    Y = restrict(filter_tangent(contrast),SEAM)[1]
    B = reference[1]
    parameter = 3*B[0,1]/(4*abs(Y[0,0]))
    noncommuting = response_receipt(contrast,parameter)
    require(Q(noncommuting['regional_commutator_hs_squared']) > 0,
            'quantum response accidentally commutes with the reference')
    names = ['model/native_gravity_response/certificate.py',
             'model/native_gravity_response/PROOF.md','model/native_gravity_response/test_certificate.py',
             'model/native_quantum_reconstruction/certificate.py',
             'model/native_quantum_reconstruction/PROOF.md','model/observer/native.py',
             'model/observer/quantum_mixing.py','model/native_thermodynamic_bridge/certificate.py',
             'model/common_source/slow_geometry.py',
             'evidence/native_geometric_source_20260925/spec.json',
             'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
             'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
             'Lean/ObserverPatchHolography/EinsteinBranch/Entropy.lean',
             'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean']
    return {'schema':'oph.native_gravity_response.v1','all_pass':True,
            'source':{'states':DIM,'seams':30,'branches':DIM*60,'return_square_count':Z,
                      'reference':'rho_2 = C^2 / Tr(C^2)','reference_spectral_floor':str(floor),
                      'return_probability':str(Q(Z,DIM*60**2)),
                      'midpoint_filter_acceptance_on_return':'1/2',
                      'unconditional_acceptance_probability':str(Q(Z,2*DIM*60**2)),
                      'source_law_changed_by_response':False},
            'regional_seams':regions,'nested_completion_tower_checked':True,
            'occupation_fixed_response':fixed,'noncommuting_response':noncommuting,
            'scope':{'count_derived_common_reference':True,'regional_entropy_split':True,
                     'quantitative_modular_first_law':True,'native_quantum_state_prepared':False,
                     'response_is_postselected_history_family':True,
                     'causal_quantum_state_intervention_executed':False,
                     'filter_reads_simulator_configuration':True,'filter_failed_events_retained':True,
                     'exterior_multiplicity_identified_with_area':False,
                     'geometric_modular_normalization':False,'physical_stress_or_Ward_identity':False,
                     'fixed_proper_volume_stationarity':False,'positive_Newton_normalization':False,
                     'controlled_small_ball_refinement':False,'Einstein_premises_G1_to_G6':False},
            'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    value=report()
    path=HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    if args.check:
        require(strict_json.load(path)==value,'native gravity response receipt is stale')
    print(json.dumps({'all_pass':True,'regional_seams':len(value['regional_seams']),
                      'Einstein_premises_G1_to_G6':False}))


if __name__ == '__main__':
    main()
