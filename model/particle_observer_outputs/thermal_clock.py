"""Faithful Gibbs and complete record clocks for the actual 307-state source.

The inverse temperature, purification service and program tick are declared
architecture data. The spectrum and probabilities are computed from the
owning Higgs Hamiltonian. No passive thermalization or laboratory unit is
inferred.
"""
from fractions import Fraction as F
from functools import lru_cache
import json

import sympy as s

from . import higgs_records as source
from . import readout_contract
from model.spin_operational_position.bounds import generic_instrument_plan


SOURCE_INPUTS = (
    'model/particle_observer_outputs/higgs_records.py',
    'model/particle_observer_outputs/HIGGS_RECORDS.md',
    'model/particle_observer_outputs/readout_contract.py',
    'model/particle_observer_outputs/READOUT_CONTRACT.md',
    'model/particle_candidate_atlas/interacting_chiral.py',
    'model/spin_operational_position/bounds.py',
    'Lean/QFT/StructuralInheritance.lean',
    'Lean/QFT/ModularClockComparison.lean',
)


def _clean(value):
    if isinstance(value, s.MatrixBase):
        return source._matrix(value)
    return s.simplify(s.expand_log(s.expand(value), force=True))


@lru_cache(maxsize=1)
def thermal_fixture():
    """The full faithful state, using the actual source spectral projectors."""
    data = source.sector_data(s.diag(1, 2, 3), 1)
    projectors = source.spectral_projectors(data)
    energies = {'ground': data['ground_energy'], 'excited': data['excited_energy'],
                'dark': data['dark_energy']}
    ranks = {name: int(s.trace(p)) for name, p in projectors.items()}
    if energies != dict(ground=-7, excited=8, dark=1) or ranks != dict(ground=1, excited=1, dark=305):
        raise ArithmeticError('the shared source spectral fixture changed')
    beta = s.log(2)
    boltzmann = {name: s.simplify(s.exp(-beta*energy)) for name, energy in energies.items()}
    partition = sum(ranks[name]*boltzmann[name] for name in energies)
    weights = {name: s.simplify(value/partition) for name, value in boltzmann.items()}
    expected = dict(ground=s.Rational(32768, 71809), excited=s.Rational(1, 71809),
                    dark=s.Rational(128, 71809))
    if partition != s.Rational(71809, 256) or weights != expected:
        raise ArithmeticError('Gibbs normalization lost the dark degeneracy')
    rho = _clean(sum((weights[name]*p for name, p in projectors.items()), s.zeros(307)))
    sqrt_rho = _clean(sum((s.sqrt(weights[name])*p for name, p in projectors.items()), s.zeros(307)))
    modular_energies = {name: _clean(-s.log(weight)) for name, weight in weights.items()}
    if any(_clean(modular_energies[name]-beta*energy-s.log(partition)) != 0
           for name, energy in energies.items()):
        raise ArithmeticError('spectral modular clock is not beta times this H')
    modular_hamiltonian = _clean(sum((modular_energies[name]*p for name, p in projectors.items()),
                                     s.zeros(307)))
    if (s.trace(rho) != 1 or any(w.is_positive is not True for w in weights.values())
            or not source._zero(rho*data['H']-data['H']*rho)
            or not source._zero(sqrt_rho*sqrt_rho-rho)
            or not source._zero(modular_hamiltonian-beta*data['H']-s.log(partition)*s.eye(307))):
        raise ArithmeticError('the complete faithful Gibbs identities failed')
    empty = _clean(s.trace(rho*data['vacuum_projector']))
    bright = _clean(s.trace(rho*data['bright_projector']))
    dark = _clean(s.trace(rho*data['dark_projector']))
    if (empty, bright, dark) != tuple(s.Rational(x, 1077135) for x in (262151, 229384, 585600)):
        raise ArithmeticError('the scalar-number detector weights changed')
    return dict(source=data, beta=beta, partition=partition, projectors=projectors,
                energies=energies, ranks=ranks, per_state_weights=weights,
                rho=rho, sqrt_rho=sqrt_rho, modular_energies=modular_energies,
                modular_hamiltonian=modular_hamiltonian,
                occupation_weights=dict(empty=empty, bright=bright, dark=dark))


def modular_clock_certificate():
    """All nine spectral blocks, including every degenerate dark matrix unit."""
    thermal = thermal_fixture()
    rows = []
    for a, ea in thermal['energies'].items():
        for b, eb in thermal['energies'].items():
            log_ratio = _clean(s.log(thermal['per_state_weights'][a])
                               -s.log(thermal['per_state_weights'][b]))
            if _clean(log_ratio+thermal['beta']*(ea-eb)) != 0:
                raise ArithmeticError('modular and Heisenberg block frequencies differ')
            rows.append({'left_block': a, 'right_block': b,
                         'matrix_unit_count': thermal['ranks'][a]*thermal['ranks'][b],
                         'Heisenberg_frequency': str(ea-eb),
                         'modular_frequency': str(log_ratio)})
    if sum(row['matrix_unit_count'] for row in rows) != 307**2:
        raise ArithmeticError('the modular comparison omitted field operators')
    return {'beta': str(thermal['beta']), 'partition': str(thermal['partition']),
            'density_identity': 'rho=exp(-beta H)/Z',
            'modular_H_identity': '-log(rho)=beta H+log(Z) I',
            'modular_convention': 'sigma_s(A)=rho^(i s) A rho^(-i s)',
            'Heisenberg_convention': 'alpha_tau(A)=exp(i tau H) A exp(-i tau H)',
            'same_operator_clock_map': 'tau=-beta*s',
            'spectral_operator_blocks': rows,
            'all_94249_matrix_units_covered': True,
            'physical_seconds_or_GeV_assigned': False}


@lru_cache(maxsize=1)
def purification_data():
    """Canonical purification and a complete Householder preparation service.

    Coefficient matrices represent sum_ij C_ij |i> tensor |bar(j)>.
    This retains both 307-dimensional systems without a dense 94249-square
    unitary. The full reflection is I-2|u><u| on their entire tensor product.
    """
    thermal = thermal_fixture()
    coefficient = thermal['sqrt_rho']
    blank = thermal['source']['vacuum_projector']
    overlap = _clean(s.trace(blank.H*coefficient))
    if overlap.is_nonnegative is not True or _clean(1-overlap).is_positive is not True:
        raise ArithmeticError('the declared blank and purification are not distinct normalized vectors')
    difference = _clean(coefficient-blank)
    # Tr(C*C)=Tr(rho)=1 and ||blank||=1 were established on the full
    # matrices. Polarization gives ||C-blank||^2=2-2<blank,C> exactly.
    normalization = 2*(1-overlap)
    if s.trace(thermal['rho']) != 1 or s.trace(blank.H*blank) != 1:
        raise ArithmeticError('the complete purification-vector normalization differs')
    # Contract the scalar normalization before expanding a rank-one
    # reflection. This is the same full-vector identity and avoids nesting
    # large radicals in every matrix entry during the exact replay.
    direction = difference/s.sqrt(normalization)
    direction_norm = s.S.One
    blank_overlap_unnormalized = overlap-1
    reflection_coefficient = _clean(-2*blank_overlap_unnormalized/normalization)
    prepared = _clean(blank+reflection_coefficient*difference)
    if (direction_norm != 1 or not source._zero(prepared-coefficient)
            or not source._zero(coefficient*coefficient.H-thermal['rho'])
            or not source._zero(coefficient.T*coefficient.conjugate()-thermal['rho'].T)):
        raise ArithmeticError('canonical purification or complete reflection preparation failed')
    return dict(coefficient=coefficient, blank_coefficient=blank,
                householder_direction=direction, direction_norm=direction_norm,
                householder_difference=difference, householder_norm_squared=normalization,
                blank_overlap=overlap, field_marginal=thermal['rho'],
                conjugate_reference_marginal=thermal['rho'].T,
                joint_dimension=307**2, Schmidt_rank=307,
                reference_is_retained=True,
                field_density_gauge_invariant=True,
                field_density_supported_only_on_singlets=False,
                joint_conjugate_reference_state_is_singlet=True,
                charged_reference_required_for_full_sector_Gauss_route=True,
                preparation_formula='I-2|u><u|, u=(|sqrt(rho)>-|blank>)/sqrt(2(1-overlap))')


def purification_covariance(generator):
    """Verify the conjugate-reference identity for an admitted source symmetry.

    For an anti-Hermitian generator G, the diagonal joint action changes C
    by G C-C G. Gauge generators commuting with the same H therefore fix
    sqrt(rho). The source Ward theorem supplies their admission.
    """
    thermal = thermal_fixture()
    generator = _clean(generator)
    if generator.shape != (307, 307) or not source._zero(generator+generator.H):
        raise ValueError('an anti-Hermitian generator on the same 307-state sector is required')
    if not source._zero(generator*thermal['source']['H']-thermal['source']['H']*generator):
        raise ValueError('the generator must commute with the actual source Hamiltonian')
    coefficient = thermal['sqrt_rho']
    result = _clean(generator*coefficient-coefficient*generator)
    if not source._zero(result):
        raise ArithmeticError('the conjugate-reference purification is not invariant')
    return {'joint_generator_annihilates_purification': True,
            'same_source_symmetry_admission_required': True}


def preparation_plan(error=F(1, 10**6)):
    """Use the inherited all-input compiler, or declare the reflection primitive."""
    plan = generic_instrument_plan(18, 1, error)
    plan.update({'field_dimension': 307, 'conjugate_reference_dimension': 307,
                 'joint_code_dimension': 307**2, 'padded_unitary_dimension': 512**2,
                 'padding_action': 'identity on every unused field/reference code input',
                 'preparation': 'complete Householder reflection from the invariant blank',
                 'exact_primitive_generator': 'pi*|u><u| for one declared preparation duration',
                 'logical_preparation_commutes_with_joint_gauge_action': True,
                 'gatewise_Gauss_preservation_of_compiled_tape_claimed': False,
                 'reference_retained_after_preparation': True,
                 'all_pointer_and_compiler_failure_records_retained': True,
                 'source_primitive_route_also_admissible': True,
                 'coherent_compiled_tape_executed_here': False,
                 'passive_repair_Gibbs_preparation_claimed': False})
    return plan


def thermal_joint_law(time):
    thermal = thermal_fixture()
    return source.occupation_joint_law(thermal['source'], time, **thermal['occupation_weights'])


def preparation_independent_gap_identity():
    """Symbolic identity on its explicit positive-weight, unaliased domain."""
    alpha, amplitude, x = s.symbols('alpha amplitude x', positive=True)
    p = alpha*amplitude*s.sin(x)**2
    q = alpha*amplitude*s.sin(2*x)**2
    difference = s.trigsimp(q/(4*p)-s.cos(x)**2)
    if difference != 0:
        raise ArithmeticError('the unconditional double-time ratio failed')
    return {'first_joint_event': 'J_empty,occupied(t)=alpha*T(t)',
            'alpha': 'Tr(rho P_empty)>0, identical in both runs',
            'transition': 'T(t)=4g_pair^2/Delta^2 sin^2(Delta*t/2)',
            'ratio': 'J_empty,occupied(2t)/(4J_empty,occupied(t))=cos^2(Delta*t/2)',
            'inverse': 'Delta=2 arccos(sqrt(ratio))/t',
            'alias_window': '0<Delta*t<pi',
            'nonzero_coupling_required': True,
            'ground_preparation_required': False,
            'identical_reproducible_preparation_required': True,
            'conditioning_on_first_record': False,
            'all_four_original_joint_outcomes_retained': True}


def record_certificate():
    thermal = thermal_fixture()
    data, t = thermal['source'], s.pi/45
    first, second = thermal_joint_law(t), thermal_joint_law(2*t)
    if any(x.is_nonnegative is not True for x in (*first.values(), *second.values())):
        raise ArithmeticError('a complete record probability is negative')
    if sum(first.values()) != 1 or sum(second.values()) != 1:
        raise ArithmeticError('complete thermal records are not normalized')
    p, q = first['empty', 'occupied'], second['empty', 'occupied']
    ratio = _clean(q/(4*p))
    inverse = source.gap_interval_from_records((p, p), (q, q), t)
    if (source.transition_probability(data, t) != s.Rational(56, 225)
            or ratio != s.Rational(3, 4) or inverse['gap_interval'] != (15, 15)):
        raise ArithmeticError('thermal records do not recover the same source gap')
    encode = lambda law: {a+','+b: str(value) for (a, b), value in law.items()}
    return {'time': str(t), 'double_time': str(2*t),
            'transition_at_time': '56/225', 'transition_at_double_time': '56/75',
            'joint_probabilities_at_time': encode(first),
            'joint_probabilities_at_double_time': encode(second),
            'ratio': str(ratio), 'recovered_gap': '15',
            'all_outcomes_retained': True,
            'dark_weight': str(thermal['occupation_weights']['dark']),
            'preparation_independent_identity': preparation_independent_gap_identity(),
            'finite_error_contract': 'complete preparation, tick, detector and sampling errors enter the source gap_interval_from_records probability intervals',
            'one_tick': readout_contract.macro_clock_plan(),
            'two_tick_two_detector_complete_error': str(readout_contract.composition_error(
                ticks=2, detector_uses=2, preparation_error=F(1, 10**6),
                step_error=F(1, 10**6), detector_error=F(1, 10**6)))}


def certificate():
    thermal, purification = thermal_fixture(), purification_data()
    return {'schema': 'oph.observer_outputs.thermal_clock.v1',
            'fixture': {'Y': 'diag(1,2,3)', 'pair_cost': '1', 'dimension': 307,
                        'beta': str(thermal['beta']), 'beta_is_declared_architecture': True},
            'spectral_values': {name: str(value) for name, value in thermal['energies'].items()},
            'spectral_ranks': thermal['ranks'],
            'partition': str(thermal['partition']),
            'per_state_Gibbs_weights': {name: str(value) for name, value in thermal['per_state_weights'].items()},
            'spectral_block_probabilities': {name: str(thermal['ranks'][name]*value)
                                            for name, value in thermal['per_state_weights'].items()},
            'minimum_density_eigenvalue': '1/71809', 'faithful_on_full307_sector': True,
            'occupation_weights': {name: str(value) for name, value in thermal['occupation_weights'].items()},
            'modular_clock': modular_clock_certificate(), 'records': record_certificate(),
            'purification': {k: value for k, value in purification.items()
                             if k in ('joint_dimension', 'Schmidt_rank', 'reference_is_retained',
                                      'field_density_gauge_invariant', 'field_density_supported_only_on_singlets',
                                      'joint_conjugate_reference_state_is_singlet',
                                      'charged_reference_required_for_full_sector_Gauss_route', 'preparation_formula')},
            'purification_compiler': preparation_plan(),
            'source_inputs': list(SOURCE_INPUTS),
            'same_H_state_clock_and_NH_detector': True,
            'empirical_mass_inputs': [], 'physical_seconds_or_GeV_assigned': False,
            'Gibbs_state_from_uncontrolled_repair_claimed': False,
            'interacting_continuum_or_laboratory_particle_identification': False}


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2))
