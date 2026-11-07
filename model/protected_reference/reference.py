"""Exact protected-charge state, preparation and current identities.

The domain composes the declared native endpoint weights with the existing
ordered-cover memory law. It is distinct from the unextended native member.
No interacting field Hamiltonian is an input or an output of this module.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from ..defect_protection import protection
from ..exchange_sector import statistics as exchange
from ..gauge_defect_transport import transport as gauge_transport
from ..local_repair_family import local_bias as local, mixing
from ..native_program import interpreter

ROOT = Path(__file__).resolve().parents[2]
SOURCE_FILES = (
    'model/protected_reference/__init__.py',
    'model/protected_reference/reference.py',
    'model/protected_reference/test_reference.py',
    'model/protected_reference/PROOF.md',
    'model/protected_reference/gauge_join.py',
    'model/protected_reference/test_gauge_join.py',
    'model/protected_reference/GAUGE_JOIN.md',
    'model/protected_reference/reversal.py',
    'model/protected_reference/test_reversal.py',
    'model/protected_reference/REVERSAL.md',
    'Lean/Thermodynamics/ProtectedGaugeCorner.lean',
    'Lean/Thermodynamics/CofinalSpectralTailFamily.lean',
    'Lean/Thermodynamics/FourLawAdequacySurface.lean',
    'model/native_program/interpreter.py',
    'model/local_repair_family/local_bias.py',
    'model/local_repair_family/mixing.py',
    'model/defect_protection/protection.py',
    'model/defect_protection/service.py',
    'model/exchange_sector/statistics.py',
    *gauge_transport.SOURCE_FILES,
)


def parameters():
    """Read the named member's endpoint, seam and law weights."""
    source = interpreter.program()
    raw = source['parameters']
    if source['carrier'] != {'kind': 'local', 'shell': 'H66'}:
        raise ArithmeticError('the named native carrier changed')
    result = local.Parameters(tuple(raw['port_weights']), tuple(raw['seam_weights']),
                              F(raw['c_weight']), F(raw['repair_probability']))
    if result.port_weights != (25, 25)+(8,)*10 or result.repair_probability != 1:
        raise ArithmeticError('the declared biased native member changed')
    if not 0 < result.c_weight < 1:
        raise ArithmeticError('the declared member requires both tagged laws')
    return result


@dataclass(frozen=True)
class Branch:
    law: str
    seam: tuple[int, int]
    coin: int
    raw: int
    abort: bool
    probability: F
    inputs: tuple[int, ...]
    outputs: tuple[int, ...]
    signs: tuple[int, ...]

    def cover_pairs(self):
        return tuple((r*66+i, (r ^ int(sign == -1))*66+j)
                     for i, j, sign in zip(self.inputs, self.outputs, self.signs, strict=True)
                     for r in (0, 1))


def branches():
    """All 720 raw slots, with their biased squared Kraus amplitudes."""
    p = parameters()
    if local.basis() != exchange.STATES:
        raise ArithmeticError('the native and ordered-cover source bases differ')
    seams = dict(zip(local.SEAMS, p.seam_weights, strict=True))
    total = sum(p.seam_weights)
    result = []
    for law in ('C', 'E'):
        for row in protection.branches(law):
            weight = (F(seams[row.seam], total)*(p.c_weight if law == 'C' else 1-p.c_weight)
                      * local.coin_probability(p, row.seam, row.coin))
            result.append(Branch(law, row.seam, row.coin, row.raw, row.abort, weight,
                                 row.inputs, row.outputs, row.signs))
    return tuple(result)


def _probability(value):
    if type(value) is not F or not 0 <= value <= 1:
        raise ValueError('an exact Fraction in [0,1] is required')
    return value


def stationary_entries(plus_weight):
    """Sparse exact rho_* tensor (a P_plus + (1-a) P_minus).

    Ordering is (orientation, source), matching the native cover. Entries
    include zero coherence at a=1/2 so no support is inferred from rounding.
    """
    a = _probability(plus_weight)
    result = {}
    for i, weight in enumerate(local.reference(parameters())):
        for r in (0, 1):
            result[r*66+i, r*66+i] = weight/2
            result[r*66+i, (1-r)*66+i] = weight*(2*a-1)/2
    return result


def apply_exact(entries):
    """Complete unread channel on exact sparse operators, without conditioning."""
    if type(entries) is not dict or any(
            type(key) is not tuple or len(key) != 2
            or any(type(i) is not int or not 0 <= i < 132 for i in key)
            or type(value) is not F for key, value in entries.items()):
        raise ValueError('an exact sparse 132 by 132 operator is required')
    output = {}
    for row in branches():
        mapping = dict(row.cover_pairs())
        for (i, j), value in entries.items():
            if i in mapping and j in mapping:
                target = (mapping[i], mapping[j])
                output[target] = output.get(target, F())+row.probability*value
    return {key: value for key, value in output.items() if value}


def modular_gap_ratio(initial_index, final_index, initial_charge, final_charge, plus_weight):
    """exp(K_final-K_initial), retaining exact rational source data.

    Logarithms exist only on the support. At a=0 or a=1 a cross-charge
    transition has no finite full-state modular energy assignment.
    """
    a = _probability(plus_weight)
    for index in (initial_index, final_index):
        if type(index) is not int or not 0 <= index < 66:
            raise ValueError('a native basis index is required')
    for charge in (initial_charge, final_charge):
        if type(charge) is not int or charge not in (-1, 1):
            raise ValueError('an exact charge +1 or -1 is required')
    charge_weight = {1: a, -1: 1-a}
    if not charge_weight[initial_charge] or not charge_weight[final_charge]:
        raise ValueError('the requested transition leaves the modular support')
    pi = local.reference(parameters())
    return pi[initial_index]*charge_weight[initial_charge]/(
        pi[final_index]*charge_weight[final_charge])


def current_certificate():
    """Join every native moving column to its CAR current and modular change."""
    p = parameters()
    pi = local.reference(p)
    effects, stationary = [F()]*66, [F()]*66
    movements = Counter()
    direction_maps = {edge: exchange.car_bilinear(edge[1], edge[0]) for edge in local.SEAMS}
    reverse_maps = {edge: exchange.car_bilinear(edge[0], edge[1]) for edge in local.SEAMS}
    columns = 0
    for row in branches():
        if len(set(row.inputs)) != len(row.inputs) or len(set(row.outputs)) != len(row.outputs):
            raise ArithmeticError('a raw source branch is not a partial injection')
        a, b = row.seam
        for i, j, sign in zip(row.inputs, row.outputs, row.signs, strict=True):
            x, y = local.basis()[i], local.basis()[j]
            endpoint, record = local.records.recorded_step(row.law, x, row.seam, row.coin)
            current = local.records.record_current(row.law, record)
            if endpoint != y:
                raise ArithmeticError('the ordered lift changed the native endpoint law')
            raw = {0: 0, 1: 1, -1: 2}[current]
            if row.law == 'E':
                raw += 4*((x[a]+x[b]) % 2)
            if raw != row.raw:
                raise ArithmeticError('complete native raw-record meaning changed')
            if any(y[k]-x[k] != current*(int(k == b)-int(k == a)) for k in range(12)):
                raise ArithmeticError('the exact native current continuity identity failed')
            if current:
                expected = (direction_maps if current == 1 else reverse_maps)[row.seam][i]
                if (j, sign) != expected:
                    raise ArithmeticError('moving negative-character column is not its CAR bilinear')
            elif sign != 1:
                raise ArithmeticError('a holding column rotates the ordering')
            ratio = F(p.port_weights[a], p.port_weights[b])**current
            if pi[i]/pi[j] != ratio:
                raise ArithmeticError('record current lost the same-state modular energy ratio')
            effects[i] += row.probability
            stationary[j] += row.probability*pi[i]
            movements[current] += 1
            columns += 1
    if effects != [F(1)]*66 or tuple(stationary) != pi:
        raise ArithmeticError('complete biased effects or stationary populations failed')
    return {'raw_slots': len(branches()), 'zero_slots': sum(not row.inputs for row in branches()),
            'abort_slots': sum(row.abort for row in branches()), 'source_columns': columns,
            'current_columns': {str(k): v for k, v in sorted(movements.items())},
            'continuity': 'n_k(output)-n_k(input) = j*(1[k=b]-1[k=a])',
            'negative_character_moving_operator': 'c_b^dagger c_a for j=+1; c_a^dagger c_b for j=-1',
            'modular_energy_change': 'j*log(theta_a/theta_b)',
            'branch_modular_eigenoperator': '[K_a,K_h] = j*log(theta_a/theta_b)*K_h',
            'complete_record_instrument_covariant_under_same_modular_clock': True,
            'one_complete_attempt_counted_including_holds': True,
            'native_motion_identified_with_modular_unitary_evolution': False}


def preparation_certificate():
    """Existing replacement bound survives unit-modulus exchange signs.

    Input is any state in one sharp charge sector, with arbitrary untouched
    reference. The bound compares the source/memory marginal after tracing
    the preparation tape. It does not erase or independently reset that tape.
    """
    p = parameters()
    minimum = mixing.population_minorization(p)
    rate = mixing.coherence_rate(p)
    if minimum['beta'] < F(33, 200) or rate > F(893, 990):
        raise ArithmeticError('the inherited native mixing constants changed')
    bound = F(65, 2)*rate**1600+F(167, 200)**80
    upper = F(-(-bound.numerator*10**12//bound.denominator), 10**12)
    if upper >= F(1, 10**6):
        raise ArithmeticError('the sharp protected-sector preparation bound failed')
    return {'attempts': 1600, 'population_block': 20,
            'population_minorization_lower': '33/200',
            'common_record_survival_upper': str(rate),
            'source_and_untouched_reference_half_trace_distance_upper': str(upper),
            'domain': 'arbitrary input in one fixed Q=+1 or Q=-1 sector',
            'target': 'rho_* tensor P_charge tensor input_reference_marginal',
            'preparation_archive_traced_only_for_this_comparison': True,
            'fresh_full_biased_source_pulse_execution_claimed': False}


def finite_thermodynamic_certificate(join):
    """Bind the proved Fin 2 heat bath to this exact declared singlet corner.

    The energy is the supported modular energy. This is not a physical
    calibration or a capacity-refinement construction.
    """
    from . import gauge_join
    edge = gauge_join.native_edge()
    mass = edge['target_reference_weight']+edge['source_reference_weight']
    weights = (edge['target_reference_weight']/mass, edge['source_reference_weight']/mass)
    if (weights != (F(25, 33), F(8, 33))
            or weights != (edge['forward_probability'], edge['reverse_probability'])):
        raise ArithmeticError('the native corner does not supply the Lean heat-bath weights')
    rows = join['Gauss_singlet_stationary_corners']
    if (len(rows) != len(gauge_transport.NAMES)
            or {row['primitive_charge'] for row in rows} != set(gauge_transport.NAMES)):
        raise ArithmeticError('the same-corner thermodynamic binding omitted a matter block')
    for row in rows:
        if ((F(row['target_probability']), F(row['source_probability'])) != weights
                or not row['two_orthonormal_singlets_form_reducing_corner']
                or not row['one_unread_attempt_exact_preparation_in_each_sharp_charge_corner']):
            raise ArithmeticError('the declared singlet channel differs from the Lean heat bath')
    # The four subnormalized outcomes have two supports for each input.
    outcomes = ((1, 0, weights[0]), (0, 0, weights[0]),
                (0, 1, weights[1]), (1, 1, weights[1]))
    kernel = tuple(tuple(sum((w for i, j, w in outcomes if (i, j) == (x, y)), F())
                         for y in range(2)) for x in range(2))
    if kernel != (weights, weights):
        raise ArithmeticError('retained corner outcomes do not give complete replacement')
    for x in range(2):
        if sum(kernel[x]) != 1:
            raise ArithmeticError('corner outcome completeness failed')
        for y in range(2):
            if weights[x]*kernel[x][y] != weights[y]*kernel[y][x]:
                raise ArithmeticError('corner detailed balance failed')
    return {
        'Lean_module': 'Lean/Thermodynamics/ProtectedGaugeCorner.lean',
        'Lean_theorem': 'OPH.Thermodynamics.ProtectedGaugeCorner.corner_four_laws',
        'same_kernel_Lean_theorem': 'OPH.Thermodynamics.ProtectedGaugeCorner.corner_same_kernel_four_laws',
        'source_coin_counts': [25, 8], 'source_coin_denominator': 33,
        'retained_outcomes': [{'input': i, 'output': j, 'probability': str(w)}
                              for i, j, w in outcomes],
        'kernel': [[str(w) for w in row] for row in kernel],
        'energy': ['log(33)-log(25)', 'log(33)-log(8)'],
        'positive_modular_gap': 'log(25/8)',
        'inverse_temperature_in_modular_units': '1', 'log_partition': '0',
        'finite_four_law_core_inputs_constructed': True,
        'nonconstant_modular_energy_is_conserved_by_replacement': False,
        'zeroth_law_domain': 'multiplier identifiability for this fixed nondegenerate spectrum',
        'fixed_two_state_Gibbs_entropy_limit': '0',
        'native_inverse_temperature_variation_implemented': False,
        'domain': 'diagonal marginal of each declared two-singlet corner at one fixed Q character',
        'record_archive_traced_only_for_this_marginal': True,
        'physical_energy_calibration_proved': False,
        'refinement_uniform_low_temperature_family_produced': False,
        'canonical_physical_premises_discharged': [],
    }


def certificate():
    from . import gauge_join, reversal
    current = current_certificate()
    algebras = {law: protection.algebra_certificate(law) for law in ('C', 'E')}
    for value in algebras.values():
        if value['commutant_complex_dimension'] != 2 or not value['complete_branch_X_commutation']:
            raise ArithmeticError('the same-carrier protected branch algebra changed')
    for a in (F(), F(1, 3), F(1, 2), F(1)):
        state = stationary_entries(a)
        if apply_exact(state) != {key: value for key, value in state.items() if value}:
            raise ArithmeticError('an exact protected-cover stationary state failed')
    ratios = {str(a): str(modular_gap_ratio(0, 0, 1, -1, a))
              for a in (F(1, 4), F(1, 2), F(3, 4))}
    joined = gauge_join.certificate()
    return {
        'schema': 'oph.protected_reference.v1',
        'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_FILES},
        'domain': 'native biased H66 endpoint law with the explicit persistent ordered-cover bit',
        'dimension': 132, 'parameters': parameters().description(),
        'current_and_energy': current,
        'stationary_family': 'rho_* tensor (a P_plus + (1-a) P_minus), 0<=a<=1',
        'complete_stationary_space_complex_dimension': 2,
        'full_stationary_classification_proof': 'faithful invariant state, Schwarz equality, exact existing branch commutant, equal adjoint nullities',
        'preserved_plus_weight': 'a=Tr(P_plus rho_input)',
        'faithful_reference_for_fixed_algebra_proof': 'rho_* tensor I2/2',
        'sharp_sector_preparation': preparation_certificate(),
        'modular_generator': '-sum_i log(theta_i)n_i + log(Z)I - log(a)P_plus - log(1-a)P_minus',
        'charge_flip': 'I66 tensor Z_orientation; QF=-FQ; F rho_plus F=rho_minus',
        'plus_to_minus_modular_energy': 'log(a/(1-a))',
        'exp_plus_to_minus_energy_controls': ratios,
        'endpoint_reference_domain': 'a=0 or 1 admits only its supported charge corner',
        'source_and_complete_classical_tape_distance_between_sharp_references': '0',
        'joint_source_memory_trace_distance_between_sharp_references': '1',
        'closed_native_exchange_loop': {
            law: {'history': [[name, list(seam), coin, raw]
                              for name, seam, coin, raw in protection.exchange_loop(law)['history']],
                  'source_return': True, 'orientation_action': 'X',
                  'modular_energy_change': '0'} for law in ('C', 'E')},
        'gauge_transport_join': joined,
        'same_corner_finite_thermodynamics': finite_thermodynamic_certificate(joined),
        'weighted_quantum_record_reversal': reversal.certificate(),
        'scope': {'full_native_record_slots_retained': True,
                  'original_record_clock_member_contains_this_memory_extension': False,
                  'sharp_charge_sector_is_a_spatially_localized_source_defect': False,
                  'charge_protection_selects_nonzero_mass': False,
                  'modular_clock_equals_repair_attempt_time_proved': False,
                  'interacting_Higgs_action_or_state_selected': False,
                  'physical_particle_identification_proved': False},
    }


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = certificate()
    path = Path(__file__).with_name('results.json')
    if args.check:
        if json.loads(path.read_text()) != result:
            raise SystemExit('protected reference receipt drift')
        print('protected reference receipt OK')
    elif args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
        print(f'wrote {path}')
    else:
        print(json.dumps(result, indent=2))
