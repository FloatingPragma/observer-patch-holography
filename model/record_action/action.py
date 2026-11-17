"""CP-04: a finite stochastic-matrix parent, not a physical Hamiltonian claim.

The exact population law gives a positive local kinetic operator. Its energy
ledger, modular-reference mismatch and record-only approximation are checked
separately. Coherent execution of this new operator is not the repair channel.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

from ..local_repair_family import inference, local_bias as local, mixing, modular_clock
from ..observer import record_recovery as records
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BASELINE = '67aa7b6f4a904a1a90ea6d77bf2f66c741f2ba0c'


def _parameters(parameters):
    # Validate original fields, including mutations of frozen dataclasses.
    mixing._parameters(parameters)
    return parameters


def _weights(weights):
    if (type(weights) is not tuple or len(weights) != 12
            or any(type(w) is not F or w <= 0 for w in weights)):
        raise ValueError('twelve positive exact Fraction weights required')
    return tuple(F(w.numerator, w.denominator) for w in weights)


def local_projectors(weights):
    """Q_e on H66: ten independent rank-one blocks per native seam.

Each block has diagonal (1-p,p), off diagonal -sqrt(p(1-p)),
and null vector (sqrt(theta_a),sqrt(theta_b)). Numerical matrices
are diagnostics; projector, positivity and ground-state claims have an
algebraic proof in README.md, with exact native-kernel binding below.
"""
    weights = _weights(weights)
    states = local.basis()
    index = {x: i for i, x in enumerate(states)}
    projectors = []
    for a, b in local.SEAMS:
        p = weights[a]/(weights[a]+weights[b])
        matrix = np.zeros((66, 66))
        for k in range(12):
            if k in (a, b):
                continue
            ia = index[tuple(int(v in (a, k)) for v in range(12))]
            ib = index[tuple(int(v in (b, k)) for v in range(12))]
            matrix[ia, ia], matrix[ib, ib] = float(1-p), float(p)
            matrix[ia, ib] = matrix[ib, ia] = -np.sqrt(float(p*(1-p)))
        projectors.append(matrix)
    return tuple(projectors)


def kinetic_parent(weights, *, seam_weights=(1,)*30, repair_probability=F(1)):
    _weights(weights)
    if (type(seam_weights) is not tuple or len(seam_weights) != 30
            or any(type(x) is not int or x <= 0 for x in seam_weights)
            or type(repair_probability) is not F or not 0 < repair_probability <= 1):
        raise ValueError('positive fixed seam weights and exact repair probability required')
    q = tuple(repair_probability*F(w, sum(seam_weights)) for w in seam_weights)
    return sum((float(w)*p for w, p in zip(q, local_projectors(weights), strict=True)),
               np.zeros((66, 66)))


def native_binding(parameters):
    """Exact all-row detailed balance, branch energy ratios and flux identities."""
    p = _parameters(parameters)
    pi = local.reference(p)
    kernel = mixing.probability_kernel(p)  # Full C/E native branches, with holds.
    if any(sum(row) != 1 for row in kernel):
        raise ArithmeticError('nonstochastic native kernel')
    if any(pi[i]*kernel[i][j] != pi[j]*kernel[j][i]
           for i in range(66) for j in range(66)):
        raise ArithmeticError('native detailed balance failed')
    if kernel != local.kernel(p):
        raise ArithmeticError('complete C/E source differs from endpoint kernel')
    visited, todo = {0}, [0]
    while todo:
        i = todo.pop()
        for j, value in enumerate(kernel[i]):
            if value and j not in visited:
                visited.add(j); todo.append(j)
    if len(visited) != 66:
        raise ArithmeticError('native support is disconnected')
    states = local.basis()
    index = {x: i for i, x in enumerate(states)}
    flux_rows, raw_rows = [], {}
    for name in ('C', 'E'):
        columns = 0
        for edge_index, (a, b) in enumerate(local.SEAMS):
            forward = backward = F()
            for i, state in enumerate(states):
                for coin in (0, 1):
                    output, rec = records.recorded_step(name, state, (a, b), coin)
                    j = index[output]
                    current = records.record_current(name, rec)
                    if tuple(y-x for x, y in zip(state, output, strict=True)) != tuple(
                            current*(int(v == b)-int(v == a)) for v in range(12)):
                        raise ArithmeticError('native current continuity changed')
                    exp_delta_k = F(p.port_weights[a], p.port_weights[b])**current
                    if pi[i]/pi[j] != exp_delta_k:
                        raise ArithmeticError('record does not determine modular energy increment')
                    mass = pi[i]*local.coin_probability(p, (a, b), coin)
                    if current == 1:
                        forward += mass
                    elif current == -1:
                        backward += mass
                    columns += 1
            if forward != backward or forward <= 0:
                raise ArithmeticError('stationary forward/reverse flux mismatch')
            if name == 'C':
                schedule = p.repair_probability*F(p.seam_weights[edge_index], sum(p.seam_weights))
                flux_rows.append({'seam': (a, b), 'forward_per_seam_attempt': forward,
                                  'backward_per_seam_attempt': backward,
                                  'per_attempt_current_second_moment': schedule*(forward+backward),
                                  'exp_modular_increment_for_positive_current': F(p.port_weights[a], p.port_weights[b])})
        raw_rows[name] = columns
    # Every off-diagonal transition is one moving endpoint on exactly one seam.
    moving = sum(kernel[i][j] != 0 for i in range(66) for j in range(i))
    if moving != 300:
        raise ArithmeticError('kinetic support no longer has ten pairs per seam')
    return {'parameters': p.description(), 'states': 66,
            'connected_support': True, 'all_rows_detailed_balance': True,
            'complete_native_columns': raw_rows, 'nonzero_undirected_hops': moving,
            'local_rank_one_blocks': 300, 'flux': flux_rows,
            'mean_modular_increment_per_attempt': '0',
            'stationary_population_entropy_production': '0',
            'positive_movement_is_not_positive_entropy_production': True,
            'archive_erasure_cost_or_physical_heat_inferred': False}


def parent_diagnostics(parameters):
    """Numerical diagnostics alongside exact support/proof certificates."""
    p = _parameters(parameters)
    weights = tuple(map(F, p.port_weights))
    pi = np.array(list(map(float, local.reference(p))))
    exact_kernel = mixing.probability_kernel(p)
    kernel = np.array([[float(x) for x in row] for row in exact_kernel])
    parent = kinetic_parent(weights, seam_weights=p.seam_weights,
                            repair_probability=p.repair_probability)
    similarity = np.eye(66)-np.sqrt(pi[:, None]/pi[None, :])*kernel
    if not np.allclose(parent, similarity, atol=2e-15, rtol=0):
        raise ArithmeticError('local kinetic parent fails native similarity map')
    values = np.linalg.eigvalsh(parent)
    ground = np.sqrt(pi)
    k = -np.log(pi)
    commutator = parent*(pi[None, :]-pi[:, None])
    witness = next(((i, j) for i in range(66) for j in range(i)
                    if parent[i, j] != 0 and local.reference(p)[i] != local.reference(p)[j]), None)
    # Pair-density term absent from the one-particle hopping lift. This is a
    # fixed hard-core interaction, not an adjustable physical coupling.
    onsite = [sum(F(p.port_weights[b if v == a else a], p.port_weights[a]+p.port_weights[b])
                  * p.repair_probability*F(w, sum(p.seam_weights))
                  for (a, b), w in zip(local.SEAMS, p.seam_weights, strict=True)
                  if v in (a, b)) for v in range(12)]
    diagonal_residuals = []
    for i, state in enumerate(local.basis()):
        occupied = tuple(v for v, n in enumerate(state) if n)
        exact_parent_diagonal = 1-exact_kernel[i][i]
        remainder = exact_parent_diagonal-sum(onsite[v] for v in occupied)
        expected = -sum(p.repair_probability*F(w, sum(p.seam_weights))
                        for edge, w in zip(local.SEAMS, p.seam_weights, strict=True)
                        if tuple(edge) == occupied)
        if remainder != expected:
            raise ArithmeticError('hard-core interaction residual changed')
        if remainder:
            diagonal_residuals.append({'pair': occupied, 'density_interaction': remainder})
    return {'local_similarity_max_error': float(np.max(abs(parent-similarity))),
            'ground_residual': float(np.linalg.norm(parent@ground)),
            'first_six_eigenvalues_diagnostic': values[:6].tolist(),
            'modular_reference_commutator_norm_diagnostic': float(np.linalg.norm(commutator, 2)),
            'modular_generator_commutator_norm_diagnostic': float(np.linalg.norm(parent*(k[None, :]-k[:, None]), 2)),
            'unequal_reference_hop_witness': None if witness is None else {
                'state_indices': witness,
                'reference_weights': [local.reference(p)[i] for i in witness]},
            'affine_modular_identification_excluded': True,
            'affine_exclusion_reason': 'K is diagonal; the parent has 300 strictly nonzero off-diagonal hops for every admitted positive law.',
            'same_reference_gibbs_identification_excluded': witness is not None,
            'exact_density_interactions': diagonal_residuals,
            'kinetic_operator_source_executed': False,
            'spectral_parameter': 'normalized native attempted events for relaxation; independently supplied parent-action parameter for coherent evolution'}


def infer_parent(counts):
    """Use only the PR1151 observer's counts, never private theta or loads.

Conditional on the already-declared simultaneous confidence event, each
true p lies in its returned interval. On [1/8,7/8], ||Q'(p)|| <= 2.
This gives a rational operator bound for the uniform-seam, no-hold parent.
No new experiment, acceptance conditioning or physical time is introduced.
"""
    learned = inference.infer_counts(counts)
    if learned['status'] != 'accept':
        return {'status': 'stop', 'reason': learned['stop_reason']}
    weights = _weights(learned['relative_weights'])
    rows = []
    for edge in local.SEAMS:
        a, b = edge
        center = weights[a]/(weights[a]+weights[b])
        lo, hi = learned['probability_intervals'][edge]
        if not F(1, 8) <= min(lo, center) <= max(hi, center) <= F(7, 8):
            return {'status': 'stop', 'reason': 'projector_lipschitz_domain'}
        bound = 2*max(abs(center-lo), abs(center-hi))
        rows.append({'seam': edge, 'inferred_probability': center,
                     'true_probability_interval': (lo, hi), 'projector_operator_error_upper': bound})
    error = sum(row['projector_operator_error_upper'] for row in rows)/30
    return {'status': 'accept', 'relative_weights': weights, 'edges': rows,
            'parent_operator_error_upper': error,
            'each_ordered_eigenvalue_absolute_error_upper': error,
            'unitary_halfdiamond_at_parent_parameter_one_upper': min(F(1), error),
            'simultaneous_statistical_failure_upper': inference.fixed_budget()['statistical_failure_upper'],
            'bound_is_on_confidence_event_without_conditioning_outputs': True,
            'ideal_comparison_retains_same_record_determined_stops': True,
            'all_stops_retained': True, 'physical_or_modular_time_identification': False,
            'private_weights_or_source_loads_consumed': False}


def _json(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(v) for v in value]
    return value


def source_bindings():
    paths = ['model/record_action/action.py', 'model/record_action/test_action.py',
             'model/record_action/README.md', 'model/local_repair_family/local_bias.py',
             'model/local_repair_family/mixing.py', 'model/local_repair_family/inference.py',
             'model/local_repair_family/modular_clock.py',
             'model/local_repair_family/record_clock_results.json',
             'model/local_repair_family/record_clock_campaign.json',
             'model/observer/native.py', 'model/observer/record_recovery.py']
    return bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})


def report():
    inherited = json.loads((ROOT/'model/local_repair_family/record_clock_results.json').read_text())
    counts = {tuple(b['seam']): (b['trials'], b['coin1_count']) for b in inherited['acquisition']['blocks']}
    learned = infer_parent(counts)
    if learned['status'] != 'accept':
        raise ArithmeticError('retained count input cannot certify kinetic parent')
    p = modular_clock.parameters()
    actual = kinetic_parent(tuple(map(F, p.port_weights)))
    estimated = kinetic_parent(learned['relative_weights'])
    error = float(np.linalg.norm(actual-estimated, 2))
    if error > float(learned['parent_operator_error_upper']):
        raise ArithmeticError('retained parent estimate exceeds its rational bound')
    return _json({'schema': 'oph.record_action.v1', 'issue': 1140, 'gate': 'CP-04',
                  'working_baseline': BASELINE, 'source_sha256': bound_pins(source_bindings()),
                  'native_binding': native_binding(p), 'parent': parent_diagnostics(p),
                  'inferred_parent': learned, 'retained_sample_operator_error_diagnostic': error,
                  'input_custody': {'uses_pinned_existing_acquisition': True,
                                   'fresh_acquisition_claimed': False,
                                   'inherited_stream_sha256': inherited['acquisition']['stream_sha256']},
                  'scope': {'integrated_gate_accepted': False, 'physical_mass_or_stress_derived': False,
                            'new_physical_comparisons': 0, 'coherent_parent_hardware_supplied': False},
                  'remaining': ['coherent execution and its complete detector/resource contract',
                                'compatible common observer clock/action and physical normalization',
                                'same localized matter state, independent dispersion and stress/Ward map',
                                'refinement/locality domain beyond the finite source']})


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = report()
    path = HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    if args.check:
        # Diagnostic floats are platform dependent; exact data/bindings are not.
        retained = json.loads(path.read_text())
        for candidate in (result, retained):
            candidate.pop('retained_sample_operator_error_diagnostic')
            for key in tuple(candidate['parent']):
                if key.endswith('_diagnostic') or key in ('local_similarity_max_error', 'ground_residual'):
                    candidate['parent'].pop(key)
        if result != retained:
            raise ArithmeticError('record-action receipt changed')
    print(json.dumps({'verified': True, 'gate_accepted': False,
                      'parent_error_upper_diagnostic': float(F(result['inferred_parent']['parent_operator_error_upper']))}))


if __name__ == '__main__':
    main()
