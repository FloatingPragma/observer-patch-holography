"""Native neutral composite records and a declared sterile spectator lift."""
import hashlib
import json
from functools import lru_cache
import sympy as s

from model.coherent_defect_wave.chiral_joint.higgs_composites import model as car
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.resolvent import field_value
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.defect_particle_classes.model import patterns
from .source import ROOT, ground_pair, encoded_value
from .vertices import action, scalar_matrix

RECEIPT = 'model/coherent_defect_wave/chiral_joint/ground_composite_response/results.json'


@lru_cache(maxsize=1)
def calculate():
    h, _, _ = action()
    phi = higgs_memory()[0]
    # A declared neutral scalar orientation is only an algebraic contact slice.
    neutral = s.zeros(30)
    for i in range(30):
        derivative = car.commutator(car.monomial(((i, False),)), h)
        for j in range(30):
            k = scalar_matrix(car.anticommutator(derivative, car.monomial(((j, False),))))
            neutral[i, j] = k[0, 2]+k[2, 0]
    assert neutral.T == -neutral and neutral.rank() == 28
    kernel = tuple(i for i in range(30) if neutral[:, i] == s.zeros(30, 1))
    assert len(kernel) == 2 and kernel[1]-kernel[0] == 15
    native_charge = representation()[2][10]+representation()[2][11]
    assert all(native_charge[i % 15, i % 15] == 0 for i in kernel)
    composite_records = []
    for spin in (0, 1):
        j = car.neutral_composite(spin)
        assert j
        for current in range(12):
            assert not car.commutator(car.gauge_current(current), j)
        derivative = car.commutator(h, j)
        cubic = [(key, value) for key, value in derivative.items() if len(key[0]) == 3]
        assert cubic
        word, row, col = cubic[0][0]
        composite_records.append({'spin': spin, 'global_gauge_commutators_zero': 12,
            'nonzero_cubic_EOM_terms': len(cubic),
            'cubic_witness': {'word': [[i, create] for i, create in word],
                              'scalar_row': row, 'scalar_col': col,
                              'coefficient': str(cubic[0][1])}})
    receipt = json.loads((ROOT/RECEIPT).read_text())
    for path, digest in receipt['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    spectrum = receipt['actual_ground_composite_spectrum']
    setting = spectrum['recorded_settings'][0]
    x = s.Symbol('x')
    p = s.Poly(setting['exact_source_minimal_polynomial'], x)
    field = ground_pair()['field']
    numerator = field.zero
    for i, trace in enumerate(setting['exact_cross_moments_target_rows_source_columns'][0]):
        coefficient = field_value(s.Poly(sum(p.nth(j)*x**(j-i-1)
                                       for j in range(i+1, p.degree()+1)), x))
        numerator += coefficient*field.convert(s.Rational(trace))
    z_composite = numerator/field_value(p.diff())
    # P_g has rank two; addition equals removal. These factors cancel in Z_J.
    success_per_setting = z_composite/32
    from fractions import Fraction
    interval = encoded_value(success_per_setting)['interval']
    supplied = setting['joint_success_probability_interval']
    assert Fraction(supplied[0]) <= Fraction(interval[0]) <= Fraction(interval[1]) <= Fraction(supplied[1])
    # The empty five-level occupation is the available gauge-singlet weight.
    empty = next(row for row in patterns() if row['mask'] == 0)
    assert empty['k'] == 0 and empty['Y'] == 0 and empty['colour'] == empty['weak'] == '1'
    lift_mode = lambda i: (i % 15)+16*(i//15)
    lifted_h = {((tuple((lift_mode(i), creation) for i, creation in word)), row, col): value
                for (word, row, col), value in h.items()}
    sterile = (15, 31)
    for index in sterile:
        c = car.monomial(((index, False),))
        assert not car.commutator(lifted_h, c)
        assert car.anticommutator(c, car.adjoint(c)) == car.monomial()
        for other in range(32):
            assert not car.anticommutator(car.commutator(c, lifted_h), car.monomial(((other, False),)))
    return {'minimal_native_neutral_orientation_contact': {
                'scalar_orientation': 'phi_1=phi_1^dagger=1, phi_0=0, only for this algebraic slice',
                'antisymmetric_contact_rank': 28, 'null_native_modes': list(kernel),
                'neutral_scalar_background_selected': False,
                'physical_massless_neutrino_inferred': False},
            'actual_gauge_singlet_neutral_composite': {
                'operator': 'J_spin=phi_0 c_(spin,11)-phi_1 c_(spin,10)',
                'full_operator_controls': composite_records,
                'source_receipt': RECEIPT,
                'source_receipt_sha256': hashlib.sha256((ROOT/RECEIPT).read_bytes()).hexdigest(),
                'complete_addition_removal_spectral_normalization': encoded_value(z_composite),
                'each_setting_original_success_probability': encoded_value(success_per_setting),
                'complete_records': spectrum,
                'literal_free_linear_neutral_EOM': False,
                'physical_neutrino_species_assigned': False},
            'declared_full32_sterile_spectator_extension': {
                'extension': 'R16=R15 direct_sum gauge singlet, H32=H30 tensor I_sterile_Fock',
                'singlet_weight_antecedent': 'empty five-level occupation, mask 0',
                'spin_major_sterile_modes': list(sterile),
                'ground': '(P_ground30/2) tensor (I_sterile_Fock/4)',
                'physical_ground_rank': 8,
                'scalar_normalization_unchanged_by_tensor_factor': True,
                'sterile_contact_coefficients': '0 for every pair and scalar orientation',
                'complete_records_per_sterile_mode': [
                    {'operation': 'addition', 'source_energy_transfer': '0', 'weight': '1/2'},
                    {'operation': 'removal', 'source_energy_transfer': '0', 'weight': '1/2'},
                    {'operation': 'other_energy', 'weight': '0'}],
                'direct_spectator_CAR_H_commutators_checked': 2,
                'sterile_pair_contact_checks': 64,
                'nonzero_neutrino_interaction_inserted': False,
                'source32_interacting_ground_copied_from_different_action': False,
                'massive_neutrino_or_oscillation_prediction': False},
            'additional_mass_mechanism': {
                'existing_optional_operator': 'model/observer_particles/neutrino_action_extension.py',
                'requires_independent_coefficient_and_scale': True,
                'optional_operator_used_in_this_calculation': False,
                'measured_neutrino_inputs': []}}
