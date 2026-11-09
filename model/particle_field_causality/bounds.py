"""Exact local influence bounds for the retained gauge-coordinate program.

All bounds concern its full Fock/coordinate parent at fixed lattice spacing
and polynomial order. A total representation-degree compression inherits
record estimates through its separate cutoff error, not exact locality.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.particle_observer_outputs import field_route


SOURCE_INPUTS = (
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/interacting_extension/chiral_lattice/source_data.py',
    'model/common_source/overlap_regulator.py',
    'model/full_fock_field/README.md',
    'model/source_family_compatibility/LOCAL_MATTER_READOUT.md',
    'Lean/QFT/FiniteCausalObserverNet.lean',
    'Lean/Geometry/RefiningLatticeCausalCone.lean',
)


def integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(name+' must be an integer at least '+str(minimum))
    return value


def rational(value, name, positive=False):
    if type(value) not in (int, F, str):
        raise ValueError(name+' requires exact rational data')
    result = F(value)
    if result < 0 or positive and result == 0:
        raise ValueError(name+' has the wrong sign')
    return result


def ball_size(radius):
    """Exact cardinality of a three-dimensional integer Manhattan ball."""
    r = integer(radius, 'radius')
    return (4*r**3+6*r**2+8*r+3)//3


def interaction_constants(degree, spacing='1', plaquette='0', sharpen_zero=True):
    """Weighted local interaction sum sup_x sum_Z |Z| ||Phi_Z|| 2^diam(Z).

    Outgoing link Hilbert spaces are attached to their tail cells. The
    general polynomial uses one rooted Hermitian star per site. At m=0
    the actual rank-two spin projectors give a sharper nearest-bond norm.
    Coordinate-independent onsite terms and electric Casimirs are removed
    in a support-preserving interaction picture and do not enter this sum.
    """
    m = integer(degree, 'Legendre degree')
    a = rational(spacing, 'spacing', True)
    kappa = rational(plaquette, 'absolute plaquette coefficient')
    inherited = field_route.polynomial_bounds(m, a)
    d = 2*m+1
    if m == 0 and sharpen_zero:
        bond = F(120, 17)/a
        rate = 24*bond+192*kappa
        rooted = 3*bond
        kind = 'actual nearest-neighbor rank30 Dirac/M15 hopping'
    else:
        bond = None
        rooted = 60*inherited['global_one_body_H_norm_upper']
        rate = rooted*ball_size(d)**2*4**d+192*kappa
        kind = 'Hermitian polynomial stars with local trace-norm bound'
    return {
        'Legendre_degree': m, 'spacing': a, 'plaquette_norm': kappa,
        'coordinate_word_degree': d, 'root_support_radius': d,
        'maximum_term_diameter': max(2*d, 2),
        'CAR_modes_per_cell': 60, 'rooted_kinetic_norm_upper': rooted,
        'nearest_bond_norm_upper': bond,
        'weighted_local_interaction_upper': rate,
        'decomposition': kind, 'outer_volume_in_bound': False,
        'electric_Casimir_norm_in_bound': False,
        'unrestricted_onsite_even_terms_in_interaction_picture': True,
        'uniform_in_polynomial_degree_claimed': False,
    }


def commutator_bound(constants, duration, distance, sites=1):
    """Norm bound divided by ||A||||B||, for disjoint homogeneous fields.

    Use the ordinary commutator if at least one field is even; for two
    odd fields use the anticommutator. Extra untouched references are kept.
    """
    t = rational(duration, 'duration')
    r = integer(distance, 'positive support distance', 1)
    count = integer(sites, 'observed cell count', 1)
    rate = rational(constants['weighted_local_interaction_upper'], 'local rate', True)
    exponential = field_route.exponential_tail_upper(2*rate*t, 1)
    raw = 2*count*exponential/F(2**r)
    return {'distance': r, 'duration': t, 'observed_cells': count,
            'exponential_argument': 2*rate*t,
            'uncapped_normalized_bound': raw,
            'normalized_graded_commutator_upper': min(F(2), raw),
            'complete_local_record_TV_upper': min(F(1), raw)}


def volume_bound(constants, duration, source_radius, observed_radius=0,
                 sites=1, interventions=0):
    """Complete local history comparison with an open finite Wilson box.

    Compare the actual polynomial of each box's Wilson operator. Their
    boundary-star differences are included. Arbitrary references and all
    outcomes are retained. The initial joint marginal on the entire
    radius-R simulation box and untouched reference is the same, and
    local apparatus acts under the same fixed policy. Matching only the
    observed cells is insufficient.
    """
    t = rational(duration, 'duration')
    r = integer(observed_radius, 'observed radius')
    R = integer(source_radius, 'source radius', r+1)
    count = integer(sites, 'observed cell count', 1)
    k = integer(interventions, 'local intervention count')
    if count > (2*r+1)**3:
        raise ValueError('observed cell count exceeds its cube')
    d = constants['root_support_radius']
    support = max(d, 2)
    distance = R-r-2*support
    if distance < 1:
        raise ValueError('source collar must exceed twice the term support radius')
    # Each differing kinetic star is the difference of two norm-bounded
    # stars. All potentially differing roots are counted, including paths
    # that leave and return through a smaller Wilson box's boundary.
    crossing_norm = (2*constants['rooted_kinetic_norm_upper']*(2*(R+d)+1)**3
                     +3*constants['plaquette_norm']*(2*(R+2)+1)**3)
    factor = commutator_bound(constants, t, distance, count)
    raw = crossing_norm*t*factor['uncapped_normalized_bound']
    local = min(F(2), raw)
    return {'source_radius': R, 'observed_radius': r, 'duration': t,
            'source_cells': (2*R+1)**3, 'boundary_separation_lower': distance,
            'boundary_perturbation_norm_upper': crossing_norm,
            'local_operator_difference_upper': local,
            'local_interventions': k,
            'complete_local_history_half_trace_upper': min(F(1), (2*k+1)*local),
            'unbounded_exterior_volume_in_bound': False,
            'initial_state_requirement': 'same entire radius-R box and untouched-reference joint marginal',
            'initial_state_error_included': False,
            'representation_cutoff_error_included': False}


def choose_box(constants, duration, accuracy, observed_radius=0, sites=1,
               interventions=0):
    accuracy = rational(accuracy, 'accuracy', True)
    if accuracy >= 1:
        raise ValueError('accuracy must be below one')
    radius = observed_radius+2*max(constants['root_support_radius'], 2)+1
    steps = 0
    while True:
        result = volume_bound(constants, duration, radius, observed_radius,
                              sites, interventions)
        steps += 1
        if result['complete_local_history_half_trace_upper'] <= accuracy:
            return {**result, 'requested_accuracy': accuracy,
                    'doubling_steps': steps, 'large_Fock_allocation_used': False}
        radius *= 2


def cutoff_record_bound(parent_bound, first_cutoff_error, second_cutoff_error):
    """Triangle inequality for two compressed complete experiments."""
    values = [rational(x, 'complete record error')
              for x in (parent_bound, first_cutoff_error, second_cutoff_error)]
    if any(x > 1 for x in values):
        raise ValueError('TV/half-trace errors must be at most one')
    return min(F(1), sum(values))


@lru_cache(maxsize=1)
def native_bond_check():
    """Exact native spin projectors and all-occupation bond norm."""
    data = field_route.source_data.load()
    gamma, _ = data['gamma']
    beta = gamma[3]
    beta = s.Matrix(beta).applyfunc(lambda x: s.nsimplify(x, rational=True))
    result = []
    for matrix in gamma[:3]:
        g = s.Matrix(matrix).applyfunc(lambda x: s.nsimplify(x, rational=True))
        p = (s.eye(4)-g)/2
        b = -s.Rational(4, 17)*beta*p
        assert p*p == p and p.H == p and p.rank() == 2
        assert b.H*b == s.Rational(16, 289)*p
        one_body = s.zeros(8)
        one_body[:4, 4:] = b
        one_body[4:, :4] = b.H
        spectrum = one_body.eigenvals()
        positive_sum = sum(e*n for e, n in spectrum.items() if e > 0)
        assert positive_sum == s.Rational(8, 17)
        result.append({'spin_projector_rank': 2,
                       'native_internal_dimension': len(data['charges']),
                       'single_internal_Fock_norm': str(positive_sum),
                       'full_Dirac_M15_Fock_bond_norm': str(15*positive_sum)})
    return result


def global_cutoff_counterexample():
    """Two distant commuting coordinate shifts fail after total-degree cut."""
    data = field_route.source_data.load()
    assert data['charges'].count(6) == 1
    basis = ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1))
    def shift(dx, dy):
        mat = s.zeros(len(basis))
        for j, (x, y) in enumerate(basis):
            target = (x+dx, y+dy)
            if target in basis:
                mat[basis.index(target), j] = 1
        return mat
    a, b = shift(1, 0), shift(0, -1)
    commutator = a*b-b*a
    eigenvalues = (commutator.H*commutator).eigenvals()
    assert max(eigenvalues) == 1
    assert commutator[basis.index((1, 0)), basis.index((0, 1))] == 1
    return {'rotor_basis': [list(x) for x in basis],
            'total_word_degree_cutoff': 1,
            'parent_disjoint_shift_commutator': '0',
            'compressed_commutator_norm': '1',
            'actual_M15_total_degree_commutator_norm_lower': '1',
            'actual_M15_witness': 'nontrivial native q6 singlet character on each of two links',
            'witness_input': [0, 1], 'witness_output': [1, 0],
            'exact_locality_of_total_degree_compression_claimed': False}


def _json(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k: _json(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(v) for v in value]
    return value


def certificate():
    source = interaction_constants(0, '1', '1/100')
    remote = commutator_bound(source, '1/64', 32)
    local = choose_box(source, '1/64', '1/1000000', interventions=2)
    assert remote['complete_local_record_TV_upper'] < F(1, 10000000)
    return _json({
        'schema': 'oph.same_field_causality.v1',
        'source_program': 'particle_observer_outputs.field_route full gauge-coordinate parent',
        'fixed_polynomial_control': source,
        'general_degree_one_control': interaction_constants(1, '1', '1/100'),
        'native_bond_norm_replay': native_bond_check(),
        'remote_complete_detector_comparison': {
            **remote,
            'prepared_field': 'all60 modes filled at origin; other cells empty; Haarconstant links',
            'detector': 'complete origin-occupation PVM, or any admitted local complete instrument',
            'remote_intervention': 'gauge-invariant unitary Wilson-plaquette phase at distance32',
            'clock': 'same Hamiltonian duration1/64 in both experiments',
            'retained_failures_and_zero_outcomes': True,
            'large_field_experiment_executed': False},
        'complete_two_intervention_volume_plan': local,
        'global_cutoff_locality_counterexample': global_cutoff_counterexample(),
        'cutoff_composition_example': str(cutoff_record_bound(
            remote['complete_local_record_TV_upper'], F(1, 10**8), F(1, 10**8))),
        'cutoff_composition_example_errors_are_assumed': True,
        'fixed_spacing_local_dynamics_constructed': True,
        'norm_continuous_dynamics_on_all_link_BH_claimed': False,
        'strict_Hamiltonian_light_cone_claimed': False,
        'uniform_refinement_or_Lorentz_boost_claimed': False,
        'physical_particle_assignments': [], 'source_predicted_GeV_masses': []})
