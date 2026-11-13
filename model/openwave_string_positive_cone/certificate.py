"""Complete finite records and rational common-scale/positivity certificates."""
from fractions import Fraction as F
from itertools import product
from math import comb
import json
from pathlib import Path
import sympy as sy
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.openwave_colour_string.source import native_transfer
from .source import (local_colour_census, uniform_bound, shared_length_attachment,
                     incidence_certificates, three_quark_control, SELECTED_KAPPA)

HERE = Path(__file__).resolve().parent


def six_quark_negative_colour_control():
    """Actual fundamental/conjugate currents on the complete81-coordinate witness sector."""
    subsets, _, g, _ = representation()
    ids = [subsets.index((c, 3)) for c in range(3)]
    T = [a.extract(ids, ids) for a in list(g[:7]) + [(g[6] + 2 * g[7]) / sy.sqrt(3)]]
    eye = sy.eye(3)
    def lift(a, slot):
        return sy.kronecker_product(*[a if i == slot else eye for i in range(4)])
    current = [[lift(a, 0) for a in T],
               [lift(-a.T, 1) + lift(-a.T, 2) for a in T],
               [lift(a, 3) for a in T]]
    vector = sy.Matrix([int(i == a and b == j) + int(i == b and a == j)
                        for i, a, b, j in product(range(3), repeat=4)])
    norm = (vector.T * vector)[0]
    assert norm == 24
    C, E = [], []
    for x in range(3):
        C.append(sum((a * a * vector for a in current[x]), sy.zeros(81, 1)))
        if x < 2:
            prefix = [sum((current[y][a] for y in range(x + 1)), sy.zeros(81)) for a in range(8)]
            E.append(sum((a * a * vector for a in prefix), sy.zeros(81, 1)))
    assert all(sum((current[x][a] * vector for x in range(3)), sy.zeros(81, 1)) == sy.zeros(81, 1)
               for a in range(8))
    assert all(x == value * vector for x, value in zip(C, (sy.Rational(4, 3), sy.Rational(10, 3), sy.Rational(4, 3))))
    assert all(x == sy.Rational(4, 3) * vector for x in E)
    assert sum(E, sy.zeros(81, 1)) - sum(C, sy.zeros(81, 1)) / 2 == -vector / 3
    # The two middle conjugates are literal antisymmetric native colour pairs.
    # E maps the missing-colour frame to the wedge basis(01,02,12).
    wedge = sy.Matrix([[0, 0, 1], [0, -1, 0], [1, 0, 0]])
    pairs = [(0, 1), (0, 2), (1, 2)]
    for a in T:
        lifted = sy.zeros(3)
        for col, (i, j) in enumerate(pairs):
            for k in range(3):
                for pair, coefficient in (((k, j), a[k, i]), ((i, k), a[k, j])):
                    if pair[0] == pair[1]:
                        continue
                    sign = 1 if pair[0] < pair[1] else -1
                    lifted[pairs.index(tuple(sorted(pair))), col] += sign * coefficient
        assert wedge.T * lifted * wedge == -a.T
    return {'native_carrier_dimension': 81, 'native_particle_number': 6,
            'integer_vector_norm_squared': 24,
            'all_eight_native_Gauss_currents_annihilate_vector_exactly': True,
            'two_literal_antisymmetric_native_pairs_realize_middle_conjugates': True,
            'colour_remainder_eigenvalue': '-1/3',
            'full_source_is_not_claimed_negative': True}


def complete_span_records(length, kappa, a, unit):
    dimension, full = comb(4 * length + 2, 3), comb(12 * length, 3)
    records = []
    for span in range(length):
        rank = 20 * length if span == 0 else (length - span) * (
            comb(4 * (span + 1) + 2, 3) - 2 * comb(4 * span + 2, 3)
            + comb(4 * (span - 1) + 2, 3))
        records.append({'span_cells': span, 'effect_rank': rank,
                        'conditional_probability': str(F(rank, dimension)),
                        'original_probability': str(F(rank, full)),
                        'string_energy_source': str(F(4, 3) * kappa * span),
                        'string_energy_GeV': str(F(4, 3) * kappa * span * unit),
                        'separation_fm': str(span * a)})
    assert sum(row['effect_rank'] for row in records) == dimension
    failure = 1 - F(dimension, full)
    assert sum(F(row['original_probability']) for row in records) + failure == 1
    return {'cells': length, 'complete_three_quark_singlet_dimension': dimension,
            'original_three_quark_occupation_dimension': full,
            'singlet_preparation_success': str(F(dimension, full)),
            'retained_coloured_preparation_failure': str(failure), 'records': records,
            'outside_prepared_number_sector_record_probability': '0',
            'outside_prepared_number_sector_effect': 'I-P_(N_quark=3,N_antiquark=0)',
            'instrument': 'Native span projectors on the complete singlet, Gauss-dressed; the original coloured preparation complement is retained in the kinematic register.',
            'mean_work_on_uniform_input_exact': '0',
            'work_scope': 'The nonselective projective map is unital. General input/branch work is the actual Hamiltonian expectation change, not asserted zero.',
            'idle_reference_retained': True}


def calculate():
    census = local_colour_census()
    attachment = shared_length_attachment()
    bound = uniform_bound(SELECTED_KAPPA)
    assert bound['number_floor'] > 0 and bound['adjusted_local_colour_coefficient'] < 0
    selected = json.loads((HERE.parent / 'openwave_spin_gradient_nuclei/results.json').read_text())
    assert {k: str(v) for k, v in bound['parameters'].items()} == selected['parameters_rational']
    transfer = native_transfer()
    # The same Gauss identity is universal on a two-cell total-colour singlet;
    # the executed N3/N6 intertwiners are complete-source controls of that identity.
    local_records = [{'number': r['number'], 'highest_weight': r['highest_weight'],
                      'probability': str(F(r['multiplicity'] * r['irrep_dimension'], 2 ** 24))}
                     for r in census['irreducible_rows']]
    assert sum(F(r['probability']) for r in local_records) == 1
    exact = lambda v: str(v) if isinstance(v, F) else [exact(x) for x in v] if isinstance(v, tuple) else v
    return {'schema': 'oph.full_gauss_Fock_positive_string.v1',
            'verdict': 'PASS: complete local-Fock positivity admits the simultaneous conditional one-length attachment',
            'scope': 'Finite open chains, zero boundary colour flux, the selected one-generation colourful matter action; no pure Yang-Mills or 4D continuum mass-gap assertion.',
            'local_colour_census': census,
            'complete_local_central_records': {'rows': local_records,
                'reference': 'Dual local Fock reference; maximally entangled gauge-singlet preparation.',
                'instrument': 'Joint native number, quadratic and cubic colour-Casimir spectral projectors, tensored with the untouched reference.',
                'all_multiplicity_spaces_retained': True, 'total_probability_exact': '1'},
            'uniform_Gauss_incidence_checks': incidence_certificates(),
            'three_quark_colour_shortcut_and_six_quark_failure': three_quark_control(),
            'actual_native_six_quark_counterexample': six_quark_negative_colour_control(),
            'complete_positive_source': {
                'source_parameters_rational': {k: str(v) for k, v in bound['parameters'].items()},
                'selected_kappa': str(SELECTED_KAPPA),
                'adjusted_local_colour_coefficient': str(bound['adjusted_local_colour_coefficient']),
                'strict_positive_kappa_upper': str(bound['strict_positive_kappa_upper']),
                'strict_positive_kappa_upper_diagnostic': float(bound['strict_positive_kappa_upper']),
                'uniform_number_floor': str(bound['number_floor']),
                'uniform_number_floor_diagnostic': float(bound['number_floor']),
                'inequality': 'H >= [mbar-|delta|/2-3 max(0,kappa/4-U)] N',
                'domain': 'All physical zero-boundary-flux Gauss states, every local Fock occupation, every chain length; other spin/flavour/multiplicity and particle-number sectors are retained.',
                'coefficientwise_positivity_required': False,
                'unique_finite_empty_matter_vacuum': True,
                'three_interacting_generations_bound_claimed': False,
                'pure_gauge_mass_gap_claimed': False},
            'same_two_cell_source_transfer': {
                'identity': '(U-kappa/2)(C0^2+C1^2)+kappa C0^2=U(C0^2+C1^2) on every total-colour singlet',
                'complete_native_controls': transfer,
                'all_120_baryon_64_meson_2520_nuclear_singlet_operators_preserved': True,
                'all_native_colour_preserving_current_matrix_elements_preserved': True,
                'three_empirical_mass_calibrations_unchanged': selected['three_hadron_mass_calibrations'],
                'axial_interval': selected['native_axial_ground']['gA_interval'],
                'nuclear_state_energy_and_graph_quadrupole': {
                    k: v for k, v in selected['selected_nuclear_record'].items()
                    if k not in ('physical_length_calibration', 'scope')},
                'nuclear_record_units': 'The inherited complete native record is in graph units; common_length_attachment supplies its explicitly calibrated femtometre map.',
                'coloured_kinematic_complement_not_identified_by_this_identity': True},
            'common_length_attachment': {k: exact(v) for k, v in attachment.items()},
            'common_length_diagnostics': {
                'kappa': float(SELECTED_KAPPA), 'a_fm': float(attachment['a_fm']),
                'required_kappa_interval': list(map(float, attachment['kappa_required_interval'])),
                'quadrupole_fm2_interval': list(map(float, attachment['quadrupole_fm2_interval'])),
                'tension_GeV_per_fm': float(attachment['tension_GeV_per_fm'])},
            'complete_native_span_records': [complete_span_records(L, SELECTED_KAPPA, attachment['a_fm'], attachment['electron_unit_GeV'])
                                             for L in (2, 3, 4, 5, 16)],
            'remaining_boundaries': ['The quadrupole and illustrative tension select an explicit retrospective physical length/coupling attachment.',
                'Two-cell identities preserve all five calibrations; larger-chain spectra and axial responses need their own source tests.',
                'A one-dimensional matter-supported static string is not the four-dimensional pure-gauge Yang-Mills mass-gap problem.']}
