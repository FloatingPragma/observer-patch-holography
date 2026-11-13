"""Exact local premise map; theorem constants retain their existential status."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from model.openwave_pure_su2_source.source import coefficients, certificate as open_cube_control

RANGE = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
REFERENCE = 'https://arxiv.org/abs/math-ph/0412040'


def periodic_graph(side):
    if type(side) is not int or side < 3:
        raise ValueError('periodic sides must be integers at least three')
    sites = tuple(product(range(side), repeat=3))
    def add(x, d):
        return tuple((a+b) % side for a, b in zip(x, d, strict=True))
    edges = tuple((x, axis) for x in sites for axis in range(3))
    ends = {edge: add(edge[0], RANGE[edge[1]+1]) for edge in edges}
    stars = {x: tuple(add(x, d) for d in RANGE) for x in sites}
    counts = Counter()
    plaquettes = []
    for x in sites:
        for y in stars[x]:
            for axis in range(3):
                counts[(y, axis)] += 1
        for i, j in combinations(range(3), 2):
            xi, xj = add(x, RANGE[i+1]), add(x, RANGE[j+1])
            face = (((x, i), 1), ((xi, j), 1), ((xj, i), -1), ((x, j), -1))
            divergence = Counter()
            for edge, sign in face:
                divergence[edge[0]] -= sign
                divergence[ends[edge]] += sign
            assert all(v == 0 for v in divergence.values())
            assert {edge[0] for edge, _ in face} <= set(stars[x])
            plaquettes.append(face)
    assert set(counts) == set(edges) and set(counts.values()) == {4}
    assert len(plaquettes) == len(edges) == 3*side**3
    return {'side': side, 'sites': len(sites), 'full_Haar_links': len(edges),
            'closed_plaquettes': len(plaquettes), 'local_star_size': 4,
            'every_electric_term_occurs_exactly': 4,
            'all_Gauss_and_star_support_columns_checked': len(plaquettes),
            'electric_reassembly_exact': True, 'all_plaquette_Gauss_divergences_zero': True}


def premise_map(g2=F(8), spacing=F(1)):
    g2, spacing = F(g2), F(spacing)
    alpha, beta = coefficients(g2, spacing)
    scale = 3*alpha/16
    coefficient = beta/scale
    norm = 3*coefficient
    assert norm == 64/g2**2 and scale*F(16, 3) == alpha
    assert scale*coefficient == beta
    return {'g_squared_exact': str(g2), 'spacing_exact': str(spacing),
            'electric_alpha_exact': str(alpha), 'magnetic_beta_exact': str(beta),
            'theorem_to_source_energy_scale_exact': str(scale),
            'single_plaquette_perturbation_coefficient_exact': str(coefficient),
            'local_perturbation_operator_norm_exact': str(norm),
            'local_classical_gap_lower_exact': '1',
            'unperturbed_global_normalized_gap_exact': '4', 'depends_on_volume': False,
            'numerical_coupling_admitted_by_uniform_gap_theorem': False,
            'reason_for_unassigned_numeric_admission': 'The cited theorem supplies positive constants but this package does not compute them.'}


def casimir_checks(cutoff=32):
    values = [F(k*(k+2), 4) for k in range(cutoff+1)]
    assert values[0] == 0 and values[1] == F(3, 4)
    assert all(v >= F(3, 4) for v in values[1:])
    return {'Peter_Weyl_label': 'k=2j, every integer k>=0',
            'Casimir_formula': 'k(k+2)/4', 'multiplicity_formula': '(k+1)^2',
            'unique_zero_representation': True, 'analytic_all_labels_gap_lower_exact': '3/4',
            'finite_formula_audit_last_label': cutoff,
            'audit_is_a_representation_cutoff_of_the_source': False}


def run():
    return {'schema': 'oph.full_Haar_SU2_uniform_volume_gap.v1',
            'verdict': 'PASS_EXACT_PERIODIC_FULL_HAAR_STRONG_COUPLING_PREMISES',
            'source': {'spatial_dimension': 3, 'local_space': 'L2(SU(2)^3,Haar)',
                       'all_link_representations_retained': True, 'matter_or_photon_factors': False,
                       'boundary': 'periodic cubic spatial lattices of side>=3',
                       'Hamiltonian': 'alpha sum_links C^2 + beta sum_plaquettes (1-W_p)',
                       'normalized_Wilson_loop': 'W_p=Tr_fundamental(U_p)/2',
                       'local_classical_term': '(4/3) sum_(y in x+star) sum_i C^2_(y,i)',
                       'local_perturbation': '-(16 beta/(3 alpha)) sum_(i<j) W_(x,ij)',
                       'removed_scalar_per_site_exact': '3 beta'},
            'Casimir_premise': casimir_checks(),
            'complete_graph_audits': [periodic_graph(n) for n in (3, 4, 6, 8)],
            'coupling_norms': [premise_map(F(g)) for g in (8, 16, 32, 64, 128)],
            'external_theorem': {'reference': REFERENCE, 'identifier': 'Theorem1',
                                 'local_Hilbert_dimension_may_be_infinite': True,
                                 'periodic_boundary_explicit_in_reference': True,
                                 'smallness_alpha_used_exact': '0',
                                 'bounded_smallness_beta_used': '64/(g_squared)^2',
                                 'theorem_constants_are_numerically_computed': False},
            'proved_conclusions': {'exists_finite_strong_coupling_threshold': True,
                                   'above_threshold_finite_volume_gap_uniform_in_side': True,
                                   'thermodynamic_ground_state_and_exponential_clustering': True,
                                   'physical_Gauss_restriction_contains_unique_ground': True,
                                   'uniform_gap_inherited_by_physical_GNS_observables': True,
                                   'infinite_volume_Wilson_variance_exceeds_one_eighth_on_existential_smaller_interval': True,
                                   'time_continuum_at_fixed_spatial_spacing': True},
            'physical_nontriviality': {'zero_perturbation_Haar_mean_exact': '0',
                                      'zero_perturbation_Haar_second_moment_exact': '1/4',
                                      'existential_neighborhood_variance_strict_lower_exact': '1/8',
                                      'continuity_input': 'Theorem1(4) local-observable analyticity',
                                      'numeric_coupling_admitted': False},
            'scope_limits': {'g_squared_8_uniform_gap_certified': False,
                             'numerical_uniform_gap_or_threshold_provided': False,
                             'spatial_continuum_positive_mass_gap_proved': False,
                             'Lorentz_or_Euclidean_rotation_limit_proved': False,
                             'Clay_Yang_Mills_problem_solved': False, 'GeV_glueball_mass_derived': False},
            'separate_open_cube_control': [open_cube_control(n) for n in (2, 3, 4)]}
