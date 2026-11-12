"""Constrained complete-source minima and retained screening alternatives."""
from functools import lru_cache
import sympy as s

from model.coherent_defect_wave.chiral_joint.interacting_charged_response.transfer import data as charged_data
from model.coherent_defect_wave.chiral_joint.interacting_charge_localization import screening
from model.coherent_defect_wave.chiral_joint.interacting_charge_localization.geometry import graph
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import REFERENCE_COST, FAMILIES
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.model import target, weight_basis
from model.coherent_defect_wave.chiral_joint.local_reference import physical
from model.coherent_defect_wave.chiral_joint.local_reference.reduction import quark_module, lepton_operators
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
from model.coherent_defect_wave.physical_attachment.chiral import representation


def encoded_interval(values):
    return [str(s.cancel(x)) for x in values]


def all_weak_keys(number, spin, hyper):
    qw, lw, sw = physical.labels()
    return tuple((i, mask, h) for i, (_, qy) in enumerate(qw)
                 for mask, (n, _, ly, magnetic) in enumerate(lw) if (n, magnetic) == (number, spin)
                 for h, (_, sy) in enumerate(sw) if qy+ly+sy == hyper)


@lru_cache(maxsize=None)
def matter_irrep_certificate(family):
    """Attach the abstract string endpoint to the actual local CAR/scalar block."""
    data = charged_data(family); root_interval = data['local']['rows'][0]['energy_interval']
    x = s.Symbol('x')
    # A common isolated polynomial root selects an actual native block.
    selected = None
    for spin in (0, 2):
        block = target(family, spin)
        polynomial = s.Poly(block['H'].charpoly(x).as_expr(), x).sqf_part()
        if polynomial.count_roots(*root_interval) == 1:
            selected = block, spin, polynomial; break
    assert selected is not None
    block, spin, polynomial = selected
    number, hyper = (4, -2) if family == 'ec' else (2, 1)
    keys = all_weak_keys(number, spin, hyper); positions = {key:i for i,key in enumerate(keys)}
    high = s.zeros(len(keys), block['basis'].cols)
    for i, key in enumerate(block['keys']):
        high[positions[key], :] = block['basis'][i, :]
    q = quark_module(); l = lepton_operators(); scalar = higgs_memory()[1]
    currents = [physical.tensor_map(keys, keys, quark=q['currents'][j],
                lepton=l['currents'][j], scalar=scalar[8+j]) for j in range(4)]
    d = len(FAMILIES[family]); k = high.cols
    # The dual native doublet has highest weight in its second coordinate.
    endpoint = high if family == 'ec' else s.Matrix.hstack(-(currents[0]-s.I*currents[1])*high, high)
    native = representation()[2]
    for j, current in enumerate(currents):
        dual = -native[8+j].extract(FAMILIES[family], FAMILIES[family]).T
        assert current*endpoint == endpoint*s.kronecker_product(dual, s.eye(k))
    gram = physical.metric_on(keys)
    assert endpoint.H*gram*endpoint == s.kronecker_product(s.eye(d), block['gram'])
    _, diagonal = block['gram'].LDLdecomposition(hermitian=False)
    assert all(value.is_positive is True for value in diagonal.diagonal())
    assert physical.hamiltonian_on(keys)*endpoint == endpoint*s.kronecker_product(s.eye(d), block['H'])
    return {'family': family, 'native_local_weight_rows': len(keys),
            'endpoint_multiplicity_dimension': k, 'native_irrep_dimension': d,
            'chosen_twice_spin': spin, 'exact_Gram_positive': True,
            'actual_local_H_intertwining_exact': True,
            'weak_and_hypercharge_dual_current_intertwining_exact': True,
            'colour_singlet_maximal_quark_module': True,
            'selected_local_energy_polynomial': str(polynomial.as_expr()),
            'selected_local_root_interval': encoded_interval(root_interval),
            'isolated_selected_root_count': 1}


@lru_cache(maxsize=None)
def family_certificate(family):
    data = charged_data(family)
    e0lo, e0hi = data['local']['e_interval']; v = data['v']; mu = REFERENCE_COST
    trial_lo, trial_hi = data['local']['rows'][0]['energy_interval']
    complete_gap = data['sector']['complete_compensated_local_gap_lower']
    trial_gap_hi = trial_hi-e0lo
    delta = e0hi+s.Rational(1171, 100)
    assert complete_gap > 0 and delta+complete_gap < 1
    assert data['neutral']['coupling'] == s.Rational(1, 2**34)
    assert data['neutral']['tau'] == 360 and v == s.Rational(360, 2**34)
    exact_minimum = False
    if data['sector']['global_compensated_ground_identified']:
        x = s.Symbol('x')
        minimum = min(data['sector']['blocks'], key=lambda row:row['ground_interval'][0])
        polynomial = s.Poly(s.sympify(minimum['complete_characteristic_polynomial'], locals={'x':x}), x).sqf_part()
        common = s.gcd(polynomial, data['local']['q'])
        low = max(minimum['ground_interval'][0], trial_lo)
        high = min(minimum['ground_interval'][1], trial_hi)
        assert low < high and common.count_roots(low, high) == 1
        assert delta+trial_gap_hi < 1
        exact_minimum = True
    rows = []
    for cell, ell in sorted(graph()[1].items()):
        absolute = (8*e0lo+ell+complete_gap+mu, 7*e0hi+trial_hi+ell+mu)
        source_excitation = (ell+complete_gap-v, ell+trial_gap_hi+v)
        binding = ((s.Integer(ell),)*2 if exact_minimum else
                   ((s.Integer(0),)*2 if ell == 0 else
                    (ell+complete_gap-trial_gap_hi, ell+trial_gap_hi-complete_gap)))
        assert absolute[0] <= absolute[1]
        assert ell == 0 or binding[0] > 0
        rows.append({'cell': cell, 'distance': ell,
            'constrained_total_energy_minimum_interval': encoded_interval(absolute),
            'source_excitation_above_actual_neutral_ground_interval': encoded_interval(source_excitation),
            'reference_cost_subtracted_explicitly': str(mu),
            'binding_cost_relative_to_origin_interval': encoded_interval(binding),
            'binding_cost_decimal_interval': [float(x) for x in binding],
            'exact_constrained_minimum_formula': '7E0+E_L+mu+distance' if exact_minimum else None})
    return {'family': family, 'complete_compensated_local_gap_lower': str(complete_gap),
            'actual_native_local_trial_gap_upper': str(trial_gap_hi),
            'extra_flux_sector_margin_above_lower_bound': str(1-delta-complete_gap),
            'actual_lowest_local_energy_identified': data['sector']['global_compensated_ground_identified'],
            'shortest_string_saturates_global_constrained_minimum': exact_minimum,
            'rows': rows, 'actual_local_endpoint': matter_irrep_certificate(family),
            'Pi_y_T_Pi_y_exactly_zero': True,
            'reference_subtraction_cancels_from_binding_difference': True,
            'actual_neutral_ground_shift_cancels_from_binding_difference': True,
            'all_nonorigin_constrained_minima_strictly_above_origin': True}


def certificate():
    rows = [family_certificate(family) for family in ('ec', 'L')]
    screened = screening.certificate()
    assert all(row['single_unscreened_charge_projectors_all_annihilate_state'] for row in screened['records'])
    return {'source': {'Hamiltonian': 'same complete N2 PW119 H0+2^-34 T+4 N_ref',
            'native_fermion_modes': 240, 'scalar_states_per_cell': 5,
            'links': 24, 'states_per_link': 119, 'reference_modes': 30,
            'onsite_Yukawa_scalar_and_link_coefficients': '1',
            'reference_cell': 0, 'reference_occupation': 1,
            'neutral_energy_reference': 'same actual unique interacting physical ground',
            'source_clock': 'exp[-it(H0+2^-34 T+4 N_ref)]'},
            'constrained_minima': rows, 'screened_cloud_counterexample': screened,
            'scope': {'same_source_constrained_binding_certified': True,
                      'stationary_charged_eigenvalue_or_free_dispersion_computed': False,
                      'Pi_y_invariant_under_actual_H': False,
                      'screened_cloud_minimum_or_string_breaking_computed': False,
                      'finite_string_cost_is_observed_Coulomb_potential': False,
                      'empirical_mass_or_length_used': False,
                      'reference_binding_removed_by_energy_zero_subtraction': False}}
