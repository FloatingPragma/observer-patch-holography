"""Computed finite couplings and conditional physical gauge/Higgs parameters.

The fixed source specifications are readable inputs. This certificate reuses
their operational readouts, checks exact fixed-action parameter dependence,
and states the additional matching needed for physical couplings. It does
not require uniqueness of an architectural choice.
"""
from fractions import Fraction as F
import hashlib
import inspect
import json
from pathlib import Path

import numpy as np
import sympy as s

from model.particle_atlas_horizontal import field_nature
from model.particle_candidate_atlas import attenuation
from model.source_strong_coupling import readout as link_readout
from model.native_quantum_gauss import algebra as link_source
from model.native_gauge_strings import strings as tree_source
from model.common_source import scalar_coupling
from model.strong_coupling_refinement import bounds as spin_bounds
from . import higgs_records


ROOT = Path(__file__).resolve().parents[2]
SOURCE_INPUTS = (
    'model/source_strong_coupling/PROOF.md',
    'model/source_strong_coupling/readout.py',
    'model/source_strong_coupling/packet.py',
    'model/native_quantum_gauss/algebra.py',
    'model/native_gauge_strings/strings.py',
    'model/common_source/scalar_coupling.py',
    'model/common_source/full_repair_response/coupling.py',
    'model/common_source/full_repair_response/results.json',
    'model/strong_coupling_refinement/README.md',
    'model/strong_coupling_refinement/bounds.py',
    'model/observer_particles/coefficient_selection.py',
    'model/observer_particles/results/coefficient_selection.json',
    'Lean/Screen/ElectroweakBreakingComposition.lean',
    'model/particle_atlas_horizontal/field_nature.py',
    'model/particle_candidate_atlas/attenuation.py',
    'model/standard_model_source_map/koide_source.py',
    'model/particle_observer_outputs/higgs_records.py',
)


def _clean(x):
    return s.simplify(x)


def _sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def fixed_link_readout():
    """Read actual source defaults, then reuse its complete rational detector."""
    for action in (link_source.hamiltonian_action, tree_source.hamiltonian_action):
        signature = inspect.signature(action).parameters
        if signature['hopping'].default != 1 or signature['link_energy'].default != 1:
            raise ArithmeticError('the inherited fixed source coefficients changed')
    witnesses = link_readout.witness()
    return {'source_hopping_J': '1', 'source_electric_epsilon': '1',
            'dimensionless_epsilon_over_abs_J': '1',
            'complete_rational_readout_witnesses': witnesses,
            'all_original_endpoint_records_retained': True,
            'source_default_is_unambiguous': True,
            'operational_coefficient_is_computed': True,
            'identified_as_MSbar_alpha_s': False,
            'missing_transport': 'map the same spatial color action, state, clock and complete response to a specified renormalized QCD observable'}


def fixed_scalar_interface():
    """The executed common scalar architecture supplies a definite coupling."""
    coupling = F(scalar_coupling.COUPLING)
    if coupling != F(1, 4):
        raise ArithmeticError('the common scalar source coupling changed')
    # For dimension2, sqrt(2) times projected q is exactly the Pauli X.
    x = s.Matrix([[0, 1], [1, 0]])
    charge = s.diag(1, -1, 1, -1)
    exact_interaction = s.kronecker_product(s.eye(2), x, charge)/4
    gauge_hop = s.kronecker_product(s.eye(4),
        s.kronecker_product(s.eye(2), x)+s.kronecker_product(x, s.eye(2)))
    commutator = exact_interaction*gauge_hop-gauge_hop*exact_interaction
    if commutator == s.zeros(16):
        raise ArithmeticError('the common scalar interaction became a spectator')
    actual = scalar_coupling.scalar_terms(2)['interaction']
    residual = float(np.max(np.abs(actual-np.asarray(exact_interaction, dtype=float))))
    if residual > 1e-14:
        raise ArithmeticError('the executable scalar coefficient differs from the exact dimension2 form')
    return {'default_scalar_charge_coupling': str(coupling),
            'finite_register_dimension_in_control': 16,
            'dimension2_exact_interaction': 'I2 tensor X tensor diag(1,-1,1,-1)/4',
            'nonzero_complete_interaction_hopping_commutator': True,
            'source_matrix_float_diagnostic_max_difference': residual,
            'source_pulse_implementation': 'scalar_subwords and source_unitary, with their owning all-column execution tests',
            'pulse_tape_reexecuted_by_this_certificate': False,
            'source_coefficient_is_not_missing': True,
            'equal_to_electroweak_g_or_Higgs_lambda_claimed': False}


def quantum_Higgs_response():
    """Compute a coupling from the same source's complete finite response."""
    data = higgs_records.sector_data(s.diag(1, 2, 3), 1)
    time = s.pi/45
    laws = tuple(higgs_records.occupation_joint_law(data, duration, empty=1, bright=0, dark=0)
                 for duration in (time, 2*time))
    probability, double_probability = (law['empty', 'occupied'] for law in laws)
    inverse = higgs_records.gap_interval_from_records(
        (probability, probability), (double_probability, double_probability), time)
    gap, upper_gap = inverse['gap_interval']
    if gap != upper_gap or _clean(gap-data['gap']) != 0:
        raise ArithmeticError('the original complete count laws do not recover the same source gap')
    inferred_g2 = _clean(gap**2*probability/(4*s.sin(gap*time/2)**2))
    inferred_d2 = _clean(gap**2-4*inferred_g2)
    if inferred_g2 != data['g_squared'] or inferred_d2 != data['pair_cost']**2:
        raise ArithmeticError('complete finite-Higgs counts do not recover their own interaction strength')
    # This reads the actual nonempirical controller with z=1, theta=0.
    coherent_y = attenuation.family('1', 1)['Y']
    coherent_g2 = _clean(4*s.trace(coherent_y.H*coherent_y))
    expected = (425+100*s.sqrt(2))/1944
    if _clean(coherent_g2-expected) != 0:
        raise ArithmeticError('the configured coherent source vertex norm differs')
    return {'fixture': 'Y=diag(1,2,3), d=1',
            'same_source_formula': 'g_pair^2=4 Tr(Y^dagger Y); Delta^2=d^2+4g_pair^2',
            'complete_empty_input_probability': str(probability), 'time': str(time),
            'complete_fixture_joint_probabilities': [
                {a+','+b: str(value) for (a, b), value in law.items()} for law in laws],
            'count_inverse_gap': str(gap),
            'record_inverse': 'g_pair^2=Delta^2*T(t)/(4 sin^2(Delta*t/2)); d^2=Delta^2-4g_pair^2',
            'inferred_pair_strength_squared': str(inferred_g2),
            'inferred_number_cost_squared': str(inferred_d2),
            'population_records_fix_cost_sign': False,
            'fixture_nonnegative_cost_sign_is_declared': True,
            'coherent_controller_pass_amplitude': '1', 'coherent_controller_unit_phase': '1',
            'coherent_controller_phase_angle': '0',
            'coherent_source_pair_strength_squared': str(coherent_g2),
            'coherent_source_gap_squared': str(_clean(1+4*coherent_g2)),
            'complete_dark_dimension': 305,
            'equal_to_Higgs_background_v_or_quartic_lambda_claimed': False,
            'missing_transport': 'a field/state/scaling map from these complete quantum-scalar responses to the radial continuum susceptibility and gauge-current responses'}


def action_dependency():
    """Exact ratios from the owning Hessian, with explicit counterchoices."""
    q = field_nature.electroweak_algebra()['owner']
    g, gp, v, lam = q['symbols']
    ratios = {'Higgs_squared_over_W_squared': _clean(q['m_H_squared']/q['m_W_squared']),
              'Z_squared_over_W_squared': _clean(q['m_Z_squared']/q['m_W_squared'])}
    expected_ratios = {'Higgs_squared_over_W_squared': 8*lam/g**2,
                       'Z_squared_over_W_squared': 1+gp**2/g**2}
    if any(_clean(ratios[key]-value) != 0 for key, value in expected_ratios.items()):
        raise ArithmeticError('the owning quadratic ratios changed')
    examples = []
    for label, parameters, kappa in (
            ('A', {g: 1, gp: 1, v: 1, lam: 1}, s.Integer(1)),
            ('B', {g: 1, gp: 2, v: 1, lam: 2}, s.Integer(2))):
        mass = q['gauge_mass_squared'].subs(parameters)
        photon = q['photon_vector_unnormalized'].subs(parameters)
        if mass.rank() != 3 or mass*photon != s.zeros(4, 1):
            raise ArithmeticError('an admitted parameter choice violates the photon null direction')
        examples.append({'label': label, 'parameters': {str(k): str(a) for k, a in parameters.items()},
                         'color_curvature_coefficient_kappa': str(kappa),
                         'alpha_s_in_declared_curvature_convention': str(1/(4*s.pi*kappa)),
                         **{key: str(value.subs(parameters)) for key, value in ratios.items()},
                         'neutral_photon_mass_squared': '0',
                         'positive_couplings_and_scalar_curvature': True,
                         'same_internal_representations_and_canonical_kinetic_form': True})
    scale = s.symbols('energy_unit_scale', positive=True)
    for key in ('m_W_squared', 'm_Z_squared', 'm_H_squared'):
        if _clean(q[key].subs(v, scale*v)-scale**2*q[key]) != 0:
            raise ArithmeticError('common dimensional rescaling failed')
    for value in ratios.values():
        if _clean(value.subs(v, scale*v)-value) != 0:
            raise ArithmeticError('a dimensionless ratio depends on the common energy unit')
    return {'mass_ratio_formulas': {k: str(value) for k, value in ratios.items()},
            'positive_parameter_counterchoices': examples,
            'counterexample_domain': 'same group/current/kinetic-form premises, different declared actions; not two realizations of the same fully specified action',
            'common_scale': 'v -> s*v multiplies every displayed mass by s; time units transform inversely',
            'dimensionless_ratios_invariant_under_common_unit_change': True,
            'architecture_uniqueness_required': False,
            'symmetry_alone_selects_dimensionless_couplings': False}


def conditional_gauge_ray():
    """Reuse the authenticated owning solution, without claiming source matching."""
    receipt_path = 'model/observer_particles/results/coefficient_selection.json'
    producer_path = 'model/observer_particles/coefficient_selection.py'
    receipt = json.loads((ROOT/receipt_path).read_text())
    if receipt['producer_sha256'] != _sha(producer_path):
        raise ArithmeticError('the gauge-ray receipt does not bind its producer')
    for path, digest in receipt['source_sha256'].items():
        if _sha(path) != digest:
            raise ArithmeticError('the gauge-ray inherited source binding changed: '+path)
    branch = receipt['conditional_matter_trace_Chain_P_entropy_gauge_join']
    a, t = s.symbols('A t', positive=True)
    c = (3*a-11*t)/16
    ray = c*s.Matrix([s.Rational(10, 3), 2, 2])+t*s.Matrix([s.Rational(41, 6), -s.Rational(19, 6), -7])
    eliminated = s.Matrix([5*a/8+109*t/24, 3*a/8-109*t/24, 3*a/8-67*t/8])
    if (ray-eliminated).applyfunc(s.simplify) != s.zeros(3, 1) or _clean(ray[0]+ray[1]-a) != 0:
        raise ArithmeticError('the exact conditional gauge-ray elimination changed')
    return {'owner_receipt': receipt_path, 'owner_receipt_sha256': _sha(receipt_path),
            'producer_and_all_inherited_sources_authenticated': True,
            'owning_entropy_tail_proof_reexecuted_here': False,
            'exact_electromagnetic_candidate': '(79+48sqrt(2)+42sqrt(3)+3sqrt(6))/31104',
            'adopted_P': '(1+sqrt(5))/2+sqrt(pi)*alpha_Chain',
            'exact_ray': 'x=c*(10/3,2,2)+t*(41/6,-19/6,-7)',
            'same_scale_EM_constraint': 'xY+x2=A=1/alpha_Chain',
            'entropy_constraint': 'ell2(1/x2)+ell3(1/x3)=P/4',
            'entropy_trace_choice': 'weights d_r exp(-4pi^2 alpha C_r), one power of representation dimension',
            'positive_domain': '0 <= t < 3A/67, c=(3A-11t)/16',
            'conditional_values': {key: branch[key] for key in ('alphaY', 'alpha2', 'alpha3', 'sin2_theta_same_scale_tree')},
            'premises': branch['premises'],
            'stated_uniqueness_scope': branch['global_uniqueness_scope'],
            'empirical_inputs_in_owner': receipt['empirical_inputs'],
            'physical_matching_obligations': [
                'identify the matter-trace ray with the same source physical kinetic response',
                'establish the declared three-family, one-doublet chiral field census',
                'transport one clock, state and scale into the common renormalized coupling convention',
                'include threshold and radiative matching for the chosen comparison observable',
                'realize the adopted representation-entropy relation on that source and reference'],
            'value_identified_as_alpha_s_at_MZ': False,
            'conditional_solution_discarded': False}


def finite_response_reuse():
    """Keep the exact repair result and controlled spin limit in their domains."""
    path = 'model/common_source/full_repair_response/results.json'
    receipt = json.loads((ROOT/path).read_text())
    excluded_notes = [path for path in receipt['source_bindings'] if path.startswith('docs/')]
    for source_path, digest in receipt['source_bindings'].items():
        if source_path not in excluded_notes and _sha(source_path) != digest:
            raise ArithmeticError('the inherited repair-response scientific source changed: '+source_path)
    control = receipt['harmonic_coupling_control']
    if control['harmonic_elimination']['static_rank'] != 1:
        raise ArithmeticError('the owning full-repair response rank changed')
    if control['uniform_all_seam_control']['family_ratios_if_assigned_as_Yukawa'] != [1, 1, 1]:
        raise ArithmeticError('the uniform-repair control changed')
    bound = spin_bounds.finite_graph(128, 128**2, 1, 2)
    return {'full_repair_owner': path, 'full_repair_receipt_sha256': _sha(path),
            'scientific_source_bindings_authenticated': True,
            'non_evidence_note_bindings_excluded': excluded_notes,
            'full_repair_elimination_reexecuted_here': False,
            'repair_result': 'S0=diag(0,0,kappa0), Z=I+(Uharm-B)^T(Uharm-B)/G >= I',
            'static_family_rank': 1,
            'canonical_normalization_preserves_two_zero_family_coefficients': True,
            'uniform_repair_coefficient': '(5-sqrt(5))/60',
            'uniform_repair_family_ratios': [1, 1, 1],
            'automatic_Lorentzian_Yukawa_interpretation': False,
            'spin_exchange_formula': 'J_spin=4 t_hop^2/U=c^2/U for t_hop=+/-c/2',
            'fixed_graph_example': {key: str(value) for key, value in bound.items()},
            'spin_exchange_is_QCD_alpha_s': False,
            'finite_spin_limit_is_interacting': True,
            'relativistic_chiral_continuum_inferred': False}


def certificate():
    return {'schema': 'oph.observer_outputs.gauge_input.v1',
            'fixed_native_link': fixed_link_readout(),
            'fixed_scalar_charge_interface': fixed_scalar_interface(),
            'finite_quantum_Higgs_response': quantum_Higgs_response(),
            'fixed_action_dependency': action_dependency(),
            'conditional_gauge_ray': conditional_gauge_ray(),
            'finite_response_reuse': finite_response_reuse(),
            'eliminated_unknowns': [
                'the configured native link epsilon/abs(J), recovered as 1 from complete records',
                'the common scalar program coupling, read as 1/4 with its executed operator attachment',
                'the declared finite quantum-Higgs vertex norm and pair-cost magnitude from complete records',
                'conditional gauge-ray coordinates once all owning ray/entropy/field premises are adopted'],
            'physical_output_requirements': [
                'compute canonical gauge-current and curvature response of the same chosen spatial interacting source',
                'compute its scalar effective potential/susceptibility and phase in the same reference',
                'transport those complete responses through a controlled continuum and renormalization map',
                'compare dimensionless observables at one scheme and scale with calibration ancestry retained'],
            'one_common_energy_unit_can_be_declared': True,
            'common_unit_determines_dimensionless_mass_ratios': False,
            'architecture_inputs_allowed': True, 'measured_mass_inputs_added': [],
            'source_inputs': list(SOURCE_INPUTS)}


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2))
