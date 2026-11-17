"""Exact anomaly multiplicities and residual family symmetry on the clock source.

These are conditional representation obstructions, not a matter/clock splice.
The fixed endpoint law and its stationary reference have a measured-clock
protocol. The separate native current construction chooses a spin/exterior
matter lift. Neither choice alone makes a propagating physical particle.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import combinations_with_replacement, product
import hashlib
import json
from math import gcd, lcm
from pathlib import Path

import sympy as s

from model.common_source import native_gauge_object as geometry
from model.common_source import native_matter_attachment as matter
from model.local_repair_family import local_bias, modular_clock
from model.observer import record_recovery
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / 'results.json'


def _clean(m):
    return s.ImmutableMatrix(m.applyfunc(lambda x: s.expand(s.radsimp(x))))


def _matrix(m):
    return [[str(x) for x in row] for row in m.tolist()]


def _primitive(row):
    """Normalize a rational linear equation without changing its nullspace."""
    values = tuple(s.Rational(x) for x in row)
    denominator = lcm(*(int(x.q) for x in values))
    integers = [int(x * denominator) for x in values]
    divisor = gcd(*integers)
    if divisor == 0:
        return tuple(integers)
    sign = next(1 if x > 0 else -1 for x in integers if x)
    return tuple(sign * x // divisor for x in integers)


@lru_cache(maxsize=1)
def current_data():
    """Use actual admitted native currents, then the existing declared lift."""
    native = matter.canonical_currents()
    # The native center is i diag(1,1,1,0,0,0). Its chosen lift is
    # i diag(1,1,1,-3/2,-3/2), so -2 gives primitive q=6Y.
    selected = tuple(mask for mask in range(32) if mask.bit_count() in (2, 4))
    currents = []
    for j, generator in enumerate(native):
        lifted = matter.spin_lift((-2 if j == 11 else 1) * generator)
        hermitian = -s.I * matter.exterior_derivative(lifted)
        currents.append(_clean(hermitian.extract(selected, selected)))
    groups = {}
    for k, mask in enumerate(selected):
        key = ((mask & 7).bit_count(), (mask >> 3).bit_count())
        groups.setdefault(key, []).append(k)
    ordered = tuple(sorted(groups))
    blocks = tuple(tuple(groups[key]) for key in ordered)
    for h in currents:
        if h != h.H:
            raise ArithmeticError('native current is not Hermitian in the chosen convention')
        for block in blocks:
            if any(h[i, j] != 0 for i in block for j in range(15) if j not in block):
                raise ArithmeticError('occupation block is not invariant under every source current')
    return tuple(currents), selected, ordered, blocks


@lru_cache(maxsize=1)
def anomaly_certificate():
    """All 364 symmetric cubic and twelve linear traces, before labels."""
    currents, selected, counts, blocks = current_data()
    equations = []
    for a in range(12):
        row = tuple(s.trace(currents[a].extract(block, block)) for block in blocks)
        equations.append({'kind': 'mixed_gravitational', 'indices': [a],
                          'coefficients': list(map(str, row))})
    for a, b, c in combinations_with_replacement(range(12), 3):
        tensor = currents[a] * (currents[b] * currents[c] + currents[c] * currents[b])
        row = tuple(s.expand(sum(tensor[k, k] for k in block)) for block in blocks)
        equations.append({'kind': 'symmetric_gauge', 'indices': [a, b, c],
                          'coefficients': list(map(str, row))})
    matrix = s.Matrix([row['coefficients'] for row in equations])
    independent = sorted({_primitive(row) for row in matrix.tolist() if any(row)})
    constraints = s.Matrix(independent)
    reduced, pivots = constraints.rref()
    kernel = constraints.nullspace()
    if len(kernel) != 1 or kernel[0] != s.ones(5, 1) or len(pivots) != 4:
        raise ArithmeticError('the actual current multiplets lost the complete-generation ray')
    if matrix * s.ones(5, 1) != s.zeros(376, 1):
        raise ArithmeticError('a complete generation has a nonzero full anomaly tensor')
    # Enumerate finite deletions as controls. The nullspace proof above,
    # rather than this finite enumeration, proves the unbounded statement.
    surviving = []
    excluded = []
    for keep in product((0, 1), repeat=5):
        residual = matrix * s.Matrix(keep)
        witness = next((j for j, value in enumerate(residual) if value != 0), None)
        if witness is None:
            surviving.append(list(keep))
        else:
            excluded.append({'retained_blocks': list(keep), 'equation_index': witness,
                             'nonzero_residual': str(residual[witness])})
    rows = []
    for key, block in zip(counts, blocks, strict=True):
        q = {currents[11][k, k] for k in block}
        if len(q) != 1:
            raise ArithmeticError('charge is not constant on an invariant block')
        weak_doublets = len(block)//2 if key[1] == 1 else 0
        rows.append({'color_occupations': key[0], 'weak_occupations': key[1],
                     'dimension': len(block), 'primitive_q6Y': str(next(iter(q))),
                     'source_masks': [selected[k] for k in block],
                     'weak_doublet_count': weak_doublets})
    doublets = sum(row['weak_doublet_count'] for row in rows)
    if doublets != 4:
        raise ArithmeticError('source weak-doublet count changed')
    return {'multiplet_order': rows, 'all_trace_equations': equations,
            'primitive_nonzero_equations': [list(row) for row in independent],
            'constraint_rank': 4, 'reduced_constraints': _matrix(reduced),
            'integer_kernel_generator': [1]*5,
            'general_solution': 'n_L - n_R = k*(1,1,1,1,1), k integer',
            'all_left_nonnegative_solution': 'n = k*(1,1,1,1,1), k nonnegative integer',
            'binary_submodule_controls': {'admitted': surviving, 'excluded_with_witness': excluded},
            'weak_doublets_per_complete_generation': doublets,
            'original_SU2_doublet_parity_per_generation': doublets % 2,
            'number_of_generations_selected': False,
            'quantum_gauge_measure_or_physical_chirality_constructed': False}


def rotation_representation():
    data = geometry.geometry()
    frame = data['frame']
    gram = _clean(frame.T * frame)
    left = _clean(gram.inv() * frame.T)
    result = []
    for g in data['rotations']:
        inverse = tuple(g.index(i) for i in range(12))
        image = frame[list(inverse), :]
        r = _clean(left * image)
        if _clean(frame*r-image) != s.zeros(12, 3) or _clean(r.T*r) != s.eye(3):
            raise ArithmeticError('source port permutation has no stated orthogonal band action')
        result.append((g, r))
    return tuple(result)


def _parameters(parameters):
    if type(parameters) is not local_bias.Parameters:
        raise ValueError('original bounded endpoint parameters required')
    return local_bias.Parameters(parameters.port_weights, parameters.seam_weights,
                                 parameters.c_weight, parameters.repair_probability)


def fixed_reference_stabilizer(parameters):
    """State symmetry must preserve the actual fixed weights and scheduler."""
    parameters = _parameters(parameters)
    seams = local_bias.SEAMS
    weights = dict(zip(seams, parameters.seam_weights, strict=True))
    result = []
    for g, r in rotation_representation():
        if (all(parameters.port_weights[g[i]] == parameters.port_weights[i] for i in range(12))
            and all(weights[tuple(sorted((g[a], g[b])))] == weights[(a, b)] for a, b in seams)):
            result.append((g, r))
    return tuple(result)


def branch_covariance(parameters, g):
    """Original complete Kraus columns, with coin and signed-record recharting."""
    parameters = _parameters(parameters)
    if g not in tuple(p for p, _ in fixed_reference_stabilizer(parameters)):
        raise ValueError('the port permutation must preserve this fixed source and scheduler')
    def permute(state):
        return tuple(state[g.index(j)] for j in range(12))
    count = 0
    for law, edge, coin in product(('C', 'E'), local_bias.SEAMS, (0, 1)):
        target = tuple(sorted((g[edge[0]], g[edge[1]])))
        orientation = 1 if g[edge[0]] < g[edge[1]] else -1
        new_coin = coin if orientation == 1 else 1-coin
        if local_bias.coin_probability(parameters, edge, coin) != local_bias.coin_probability(parameters, target, new_coin):
            raise ArithmeticError('residual action changes the original Kraus amplitude')
        for state in local_bias.basis():
            output, record = record_recovery.recorded_step(law, state, edge, coin)
            image, new_record = record_recovery.recorded_step(law, permute(state), target, new_coin)
            if image != permute(output) or new_record != (*target, new_coin, orientation*record[-1]):
                raise ArithmeticError('complete native branch lost its stated record-covariance map')
            count += 1
    return count


@lru_cache(maxsize=1)
def family_certificate():
    parameters = modular_clock.parameters()
    source = rotation_representation()
    stabilizer = fixed_reference_stabilizer(parameters)
    identity = tuple(range(12))
    if len(stabilizer) != 2:
        raise ArithmeticError('canonical clock source residual group changed')
    g, involution = next((g, r) for g, r in stabilizer if g != identity)
    if involution != s.diag(-1, 1, -1):
        raise ArithmeticError('the actual source frame has a different residual involution')
    plus = (s.eye(3) + involution)/2
    minus = (s.eye(3) - involution)/2
    units = tuple(s.eye(3)[:, i] * s.eye(3)[j, :] for i, j in product(range(3), repeat=2))
    equations = s.Matrix.hstack(*(s.Matrix(9, 1, list(involution*e-e*involution)) for e in units))
    kernel = equations.nullspace()
    full_dimension = s.simplify(sum(s.trace(r)**2 for _, r in source)/60)
    if len(kernel) != 5 or full_dimension != 1:
        raise ArithmeticError('full or residual family commutant dimension changed')
    # Check fixed-reference and complete population-law covariance, rather
    # than inferring invariance from a geometric label alone.
    pi = local_bias.reference(parameters)
    states = local_bias.basis()
    index = {state: j for j, state in enumerate(states)}
    transition = local_bias.kernel(parameters)
    pair_permutation = []
    for state in states:
        image = tuple(state[g.index(j)] for j in range(12))
        pair_permutation.append(index[image])
    if any(pi[pair_permutation[i]] != pi[i] for i in range(66)):
        raise ArithmeticError('residual action fails actual stationary reference')
    if any(transition[pair_permutation[i]][pair_permutation[j]] != transition[i][j]
           for i, j in product(range(66), repeat=2)):
        raise ArithmeticError('residual action fails actual complete population law')
    branch_columns = branch_covariance(parameters, g)
    # A source-computed family response anisotropy, before any Yukawa map.
    # The conditional identification of this response triplet with families
    # remains a premise. The actual stationary occupations are not fitted.
    occupation = tuple(sum(p*state[j] for p, state in zip(pi, states, strict=True))
                       for j in range(12))
    frame = geometry.geometry()['frame']
    gram = _clean(frame.T*frame)
    anisotropy = _clean(gram.inv()*frame.T*s.diag(*occupation)*frame)
    expected = s.diag((27765-1649*s.sqrt(5))/150100,
                      (27765+1649*s.sqrt(5))/150100, s.Rational(976, 7505))
    if anisotropy != expected or len(anisotropy.eigenvals()) != 3:
        raise ArithmeticError('stationary source response anisotropy changed')
    if not all(value.is_positive is True for value in anisotropy.diagonal()):
        raise ArithmeticError('stationary source response is not strictly positive')
    if anisotropy*involution != involution*anisotropy:
        raise ArithmeticError('derived source moment violates its residual symmetry')
    # Retain a positive, nondegenerate example with genuine two-family
    # mixing, and a symmetry-breaking example with a nonzero CP-odd trace.
    yu = s.diag(1, 2, 3)
    yd = s.Matrix([[2, 0, 1], [0, 5, 0], [1, 0, 2]])
    hu, hd = yu*yu.H, yd*yd.H
    commutator = hu*hd-hd*hu
    if yu*involution != involution*yu or yd*involution != involution*yd:
        raise ArithmeticError('positive split control violates its fixed source symmetry')
    if commutator == s.zeros(3) or commutator.det() != 0 or s.trace(commutator**3) != 0:
        raise ArithmeticError('two-family mixing control has the wrong CP-odd invariant')
    broken = s.Matrix([[2, 1, 1], [1, 3, s.I], [1, -s.I, 4]])
    hb = broken*broken.H
    cb = hu*hb-hb*hu
    if broken*involution == involution*broken or cb.det() == 0:
        raise ArithmeticError('symmetry-breaking control failed to escape the obstruction')
    # A symbolic block proof covers every complex invariant Yukawa,
    # beyond the exhibited nondegenerate examples.
    a, b, c, d, e = s.symbols('a b c d e', complex=True)
    u, v, w, x, y = s.symbols('u v w x y', complex=True)
    first = s.Matrix([[a, 0, b], [0, c, 0], [d, 0, e]])
    second = s.Matrix([[u, 0, v], [0, w, 0], [x, 0, y]])
    h1, h2 = first*first.H, second*second.H
    bracket = h1*h2-h2*h1
    if bracket[1, :] != s.zeros(1, 3) or bracket[:, 1] != s.zeros(3, 1):
        raise ArithmeticError('general residual intertwiner did not preserve the common one-dimensional block')
    return {'clock_parameters': parameters.description(),
            'fixed_source_proper_stabilizer': [list(g) for g, _ in stabilizer],
            'residual_family_involution': _matrix(involution),
            'positive_projector': _matrix(plus), 'negative_projector': _matrix(minus),
            'full_A5_commutant_complex_dimension': int(full_dimension),
            'fixed_clock_C2_commutant_complex_dimension': len(kernel),
            'invariant_matrix_pattern': [['a', '0', 'b'], ['0', 'c', '0'], ['d', '0', 'e']],
            'derived_source_anisotropy': {
                'stationary_one_port_occupations': list(map(str, occupation)),
                'formula': '(F^T F)^(-1) F^T diag(E_pi[x]) F',
                'matrix': _matrix(anisotropy),
                'distinct_positive_eigenvalues_in_source_axis_order': list(map(str, anisotropy.diagonal())),
                'trace': str(s.trace(anisotropy)),
                'physical_family_or_Yukawa_map_selected': False,
                'all_Yukawas_functions_of_this_one_matrix_implies_no_mixing': True},
            'exact_reference_entries_checked': 66, 'exact_population_kernel_entries_checked': 66**2,
            'complete_C_E_branch_columns_covariant_with_record_recharting': branch_columns,
            'distinct_mass_and_two_family_mixing_control': {
                'Y_u': _matrix(yu), 'Y_d': _matrix(yd),
                'up_squared_singular_values': ['1', '4', '9'],
                'down_squared_singular_values': ['1', '9', '25'],
                'commutator': _matrix(commutator), 'commutator_determinant': '0'},
            'symmetry_breaking_control': {'Y_d': _matrix(broken),
                'commutator_determinant': str(s.expand(cb.det())),
                'trace_commutator_cubed': str(s.expand(s.trace(cb**3))),
                'selected_by_record_clock': False},
            'universal_same_action_conclusion': 'common 1+2 family decomposition; det[Hu,Hd]=tr([Hu,Hd]^3)=0',
            'nondegenerate_three_family_unitary_mixing': 'one isolated family and a 2x2 unitary block up to permutations; Jarlskog invariant zero',
            'degenerate_mass_case': 'basis-independent commutator obstruction remains; mixing angles within degenerate spaces are not uniquely identified',
            'fixed_source_symmetry_supplies_full_A5_degree_obstruction_hypotheses': False,
            'inferred_finite_record_estimates_alone_prove_exact_group': False,
            'family_identification_or_Yukawa_selection_derived': False}


def validate_multiplicities(left, right=None):
    """Exact unbounded family test on this five-block representation only."""
    if right is None:
        right = (0,)*5
    for value in (left, right):
        if (type(value) is not tuple or len(value) != 5
                or any(type(n) is not int or n < 0 for n in value)):
            raise ValueError('five original nonnegative integer block multiplicities required')
    signed = tuple(a-b for a, b in zip(left, right, strict=True))
    constraints = s.Matrix(anomaly_certificate()['primitive_nonzero_equations'])
    residual = constraints*s.Matrix(signed)
    doublets = sum(row['weak_doublet_count']*(a+b)
                   for row, a, b in zip(anomaly_certificate()['multiplet_order'], left, right, strict=True))
    return {'signed_multiplicities': list(signed), 'all_local_traces_zero': not any(residual),
            'original_SU2_doublet_parity': doublets % 2,
            'complete_signed_generation_count': signed[0] if len(set(signed)) == 1 else None}


def source_bindings():
    paths = ['model/particle_obstruction/classification.py',
        'model/common_source/native_gauge_object.py', 'model/common_source/slow_geometry.py',
        'model/common_source/native_sm_current.py', 'model/common_source/native_matter_attachment.py',
        'model/local_repair_family/local_bias.py', 'model/local_repair_family/modular_clock.py',
        'model/observer/native.py', 'model/observer/record_recovery.py',
        'model/local_repair_family/inference.py', 'model/local_repair_family/record_clock_campaign.json',
        'model/local_repair_family/record_clock_results.json',
        'model/observer_particles/family_naturality.py',
        'Lean/Screen/A5FamilyBand.lean', 'Lean/Screen/WeylYukawaConventions.lean',
        'evidence/native_geometric_source_20260925/spec.json',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py']
    return bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})


def report():
    return {'schema': 'oph.particle-obstruction.v1',
            'campaign': 'OPH-CURRENT-AND-CLOCK-FAMILY-2026-10-29-v1',
            'issue': 'https://github.com/FloatingPragma/observer-patch-holography/issues/1141',
            'pr1151_baseline': '67aa7b6f4a904a1a90ea6d77bf2f66c741f2ba0c',
            'source_sha256': bound_pins(source_bindings()), 'anomaly_classification': anomaly_certificate(),
            'clock_family_symmetry': family_certificate(),
            'scope': {'integrated_gate_accepted': False, 'physical_comparison_count_added': 0,
                      'anomaly_quantifier': 'vanishing fermionic trace contributions for arbitrary nonnegative integer copies of the five actual source-current blocks, with explicit left/right Weyl signs; no extra anomaly-canceling fields, Green-Schwarz terms or inflow',
                      'family_quantifier': 'all complex invariant Yukawa intertwiners for the same residual source triplet on left and right; family-singlet Higgs and invariant positive kinetic metrics',
                      'remaining': ['same moving defect carrying this matter lift and current',
                                    'physical Weyl/chirality and gauge measure',
                                    'selected common kinetic action and calibrated physical clock',
                                    'family assignment and symmetry-breaking dynamics',
                                    'refinement and complete detector correspondence']}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = report()
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2)+'\n')
    if args.check and json.loads(OUTPUT.read_text()) != result:
        raise SystemExit('retained particle-obstruction receipt differs from fresh exact replay')
    print(json.dumps({'anomaly_constraint_rank': result['anomaly_classification']['constraint_rank'],
                      'clock_source_stabilizer_order': len(result['clock_family_symmetry']['fixed_source_proper_stabilizer']),
                      'clock_family_commutant_dimension': result['clock_family_symmetry']['fixed_clock_C2_commutant_complex_dimension'],
                      'integrated_gate_accepted': False}))


if __name__ == '__main__':
    main()
