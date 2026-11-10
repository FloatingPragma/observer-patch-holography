"""Same-native-H elementary/Higgs cross residue and complete interference records."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys

import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from model.native_abelian_field_join import certificate as native
from model.coherent_defect_wave.chiral_joint.higgs_composites import model as poly
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_joint.peter_weyl_links.clebsch import irreps, native_charges
from model.particle_common_contract.provenance import import_closure

NATIVE = 'model/native_abelian_field_join/results.json'
NATIVE_SHA = 'f6aadcd79eec889de0c3243925c950e482bd456ea5c22fcab6a8406978472327'
HIGGS = 'model/coherent_defect_wave/chiral_joint/higgs_composites/results.json'
HIGGS_SHA = '6bac059d0b1b2c86d6f7915783c3964c6ebfbc877031d636fe814b18cb412441'
CHARGED = 'model/coherent_defect_wave/chiral_joint/interacting_charged_response/results.json'
PREPARATION = 'model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json'
GEOMETRY = 'model/coherent_defect_wave/chiral_joint/atomic_spin_geometry/results.json'
RESEARCH = HERE/'research_reference.json'
RESEARCH_SHA = '59902a655b4d314a32681445559858e076df1be29b522eb0398339b77cf4ff43'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(pins):
    for name, expected in pins.items():
        p = Path(name)
        if p.is_absolute() or '..' in p.parts or digest(ROOT/p) != expected:
            raise ValueError('source custody mismatch: '+name)


@lru_cache(None)
def source():
    if digest(ROOT/NATIVE) != NATIVE_SHA or digest(ROOT/HIGGS) != HIGGS_SHA:
        raise ValueError('the frozen native source or full-domain Higgs operator receipt changed')
    n = json.loads((ROOT/NATIVE).read_text()); h = json.loads((ROOT/HIGGS).read_text())
    validate(n['input_sha256']); validate(h['input_sha256'])
    if native.certificate() != n: raise ValueError('fresh actual native source replay differs')
    # Only full-domain algebra/compiler fields are consumed. No two-cell
    # periodic reference, residue or eigenvalue is transferred.
    if poly.tensor_certificate() != h['exact_operator_identities']:
        raise ValueError('actual full-domain current/spin/EOM identities differ')
    if poly.compiler_certificate() != h['compiler_ancestry']:
        raise ValueError('actual complete Yukawa compiler differs')
    charged = json.loads((ROOT/CHARGED).read_text())['same_H_interacting_channels']['ec']
    geometry = json.loads((ROOT/GEOMETRY).read_text())['full_interacting_transfer']
    prep = json.loads((ROOT/PREPARATION).read_text())['preparation']
    assert F(charged['source']['hopping_coefficient']) == F(1, 2**34)
    assert F(geometry['full_hopping_norm_upper']) == 360
    assert charged['source']['onsite_Yukawa_scalar_link_coefficients'] == '1'
    return n, h, charged, prep


@lru_cache(None)
def operator_controls():
    _, h_saved, _, _ = source()
    parts, ys = poly.quantum_yukawa()
    HY = poly.add(*parts.values())
    raw_weights = [1 if j%15 in (10, 11) else -1 if j%15 == 9 else 0 for j in range(30)]
    raw = poly.one_body(s.diag(*raw_weights))
    assert not poly.commutator(raw, HY)
    # Actual native coordinates never connect different raw-L multiplets.
    native_blocks = []
    for q in native_charges():
        indices = irreps()[q]['native_indices']
        values = {raw_weights[j] for j in indices}
        assert len(values) == 1
        native_blocks.append({'integer_charge_6Y': q, 'native_indices': list(indices), 'raw_L_weight': values.pop()})
    gauge_checks = spin_checks = 0
    probes = []
    for alpha in (0, 1):
        A = poly.monomial(((15*alpha+9, True),))
        B = poly.screened_eom('charged_lepton', 9, alpha)
        assert poly.commutator(HY, A) == poly.scale(B, ys['charged_lepton'])
        for probe in (A, B):
            assert poly.commutator(raw, probe) == poly.scale(probe, -1)
            assert all(len(word)%2 == 1 for word, _, _ in probe)
            for j in range(12):
                assert poly.commutator(poly.gauge_current(j), probe) == (probe if j == 11 else {})
                gauge_checks += 1
            for axis in range(3):
                expected = poly.add(*[poly.scale(
                    poly.monomial(((15*beta+9, True),)) if probe is A else poly.screened_eom('charged_lepton', 9, beta),
                    poly.PAULI[axis][beta, alpha]/2) for beta in (0, 1)])
                assert poly.commutator(poly.spin_current(axis), probe) == expected
                spin_checks += 1
        probes.append({'spin': alpha, 'elementary': f'c_dagger[{15*alpha+9}]',
                       'composite': '-sum_beta epsilon[alpha,beta] (phi0_dagger c[15beta+10]+phi1_dagger c[15beta+11])',
                       'full_Yukawa_equation_exact': True})
    return {'all_local_occupation_sectors': list(range(31)), 'native_scalar_dimension': 5,
            'actual_gauge_component_checks': gauge_checks, 'actual_spin_component_checks': spin_checks,
            'full_Yukawa_raw_lepton_commutator_zero': True,
            'actual_PW_hopping_blocks_preserve_raw_L': native_blocks,
            'source_creation_probe_quantum_numbers': {'color': 'singlet', 'weak': 'singlet', 'Y': '1',
                'integer_6Y': '6', 'raw_lepton_number_change': '-1', 'fermion_parity': 'odd',
                'origin_proper_cubic_type': 'G1 relative to the trivial actual ground character'},
            'conjugate_probe_quantum_numbers': {'color': 'singlet', 'weak': 'singlet', 'Y': '-1',
                'integer_6Y': '-6', 'raw_lepton_number_change': '1', 'fermion_parity': 'odd'},
            'probes': probes, 'composite_full_operator_norm_upper': '2',
            'elementary_full_operator_norm': '1',
            'scalar_expectation_or_canonical_composite_CAR_assumed': False,
            'full_domain_operator_parent_checks': h_saved['exact_operator_identities']}


def down(value, digits=40):
    value = F(value); scale = 10**digits
    return str(F(value.numerator*scale//value.denominator, scale))


def up(value, digits=40):
    value = F(value); scale = 10**digits
    return str(F(-(-value.numerator*scale//value.denominator), scale))


def cross_bound(energy_low, weight_low, hopping_coefficient, hopping_norm):
    L, w, lam, tau = map(F, (energy_low, weight_low, hopping_coefficient, hopping_norm))
    if not (L > 0 and 0 < w <= 1 and lam >= 0 and tau >= 0):
        raise ValueError('positive spectral window/residue and nonnegative source perturbation required')
    # ||P A g||<=||A||=1 and ||[T,A]||<=2||T|| are full-source bounds.
    return L*w-2*lam*tau


def bounds():
    n, _, charged, prep = source()
    low, high = map(F, n['intervals']['native_selected_E_interval'])
    wlo, whi = map(F, n['whole_bin_coupling']['bin_weight_unnormalized_interval'])
    old = charged['full_interacting_energy_readout']
    assert [str(wlo), str(whi)] == old['conditional_probabilities']['line']
    lam = F(charged['source']['hopping_coefficient']); tau = F(360)
    assert lam*tau == F(old['full_hopping_perturbation_norm'])
    gamma = cross_bound(low, wlo, lam, tau)
    assert 0 < gamma < 1 and whi < 1
    trace = F(prep['finite_energy_agreement_filter']['conditional_one_source_trace_distance_upper'])
    success = F(prep['finite_energy_agreement_filter']['success_probability_lower'])
    conditional_contrast = gamma/2-2*trace
    joint_contrast = gamma/4-trace
    assert joint_contrast > 0 and 0 < success <= 1
    absolute_source_bin = [F(x)-4 for x in old['selected_total_absolute_energy_bin']]
    return {'actual_source_H': 'H0+2^-34 T on the complete native eight-cell PW119 source',
            'actual_state': 'same unique neutral full physical ground g',
            'actual_clock': 'exp(-it h); no two-cell or free-macro clock substitution',
            'absolute_native_source_spectral_bin': list(map(str, absolute_source_bin)),
            'excitation_energy_enclosure_for_that_exact_bin': list(map(str, (low, high))),
            'elementary_unnormalized_bin_weight_interval': list(map(str, (wlo, whi))),
            'full_hopping_norm_upper': str(tau), 'source_hopping_coefficient': str(lam),
            'full_commutator_norm_upper': str(2*tau),
            'real_elementary_composite_cross_residue_lower_exact': str(gamma),
            'composite_unnormalized_bin_weight_lower_exact': str(gamma**2/whi),
            'ideal_conditional_interference_hit_contrast_lower': down(gamma/2),
            'ideal_joint_setting_hit_contrast_lower': down(gamma/4),
            'actual_prepared_conditional_interference_hit_contrast_lower': down(conditional_contrast),
            'actual_prepared_joint_setting_hit_contrast_lower': down(joint_contrast),
            'unconditional_full_preparation_setting_hit_contrast_lower': down(success*joint_contrast),
            'native_preparation_success_lower': down(success),
            'native_preparation_failure_upper': up(1-success),
            'native_preparation_conditional_trace_error_upper': up(trace),
            'complete_detector_law': {
                'setting_probability_each': '1/2', 'pump_M_pm': '(A +/- B/2)/2 with physical even port attachment',
                'success_inside': '||P_J M_pm g||^2', 'success_outside': '||(I-P_J)M_pm g||^2',
                'failure': '1-||M_pm g||^2', 'exact_completeness_each_setting': '1',
                'retained_apparatus_spin_doublet_dimension': 4,
                'retained_energy_flag_dimension': 2, 'retained_setting_dimension': 2,
                'regular_native_cubic_reference_dimension': 48,
                'full_setting_spin_occupation_energy_record_count': 16,
                'potentially_nonzero_records_from_declared_empty_apparatus': 6,
                'all_other_apparatus_records_retained_with_exact_zero': True,
                'complete_preparation_filter_pointer_dimension': int(prep['finite_energy_agreement_filter']['explicit_integer_M']),
                'preparation_partners_spin_references_and_all_filter_labels_retained': True,
                'spectral_energy_bin_readout_is_exact_finite_service_not_QPE_implemented': True,
                'autonomous_native_local_controller_claimed': False},
            'positive_cross_selector': {
                'definition': 'energy in the frozen bin maximizing Re <A_up g,P_E B_up g>; tie goes to smaller energy',
                'existence_and_positive_residue_proved': True,
                'same_energy_elementary_composite_G1_doublet': True,
                'source_only_energy_interval': list(map(str, (low, high))),
                'is_the_maximum_elementary_residue_E_star_selector': False,
                'maximum_elementary_E_star_common_composite_overlap_proved': False,
                'selector_identified_by_numerical_full_native_diagonalization': False,
                'physical_SM_weak_chirality_interactions_or_GeV_identity_proved': False}}


def input_hashes():
    n, h, _, _ = source(); seeds = set(n['input_sha256'])|set(h['input_sha256'])
    seeds.update((NATIVE, HIGGS, CHARGED, PREPARATION, GEOMETRY,
                  'model/particle_common_contract/provenance.py', 'model/requirements.txt'))
    seeds.update(str(path.relative_to(ROOT)) for path in HERE.iterdir() if path.suffix in ('.py', '.md'))
    seeds.add(str(RESEARCH.relative_to(ROOT)))
    paths = import_closure(seeds)
    return {name: digest(ROOT/name) for name in sorted(paths)}


def certificate():
    result = {'schema': 'oph.native_higgs_spectral_join.v1', 'input_sha256': input_hashes(),
            'local_sha256': {p.name: digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py', '.md')},
            'exact_full_source_operator_controls': operator_controls(), 'same_H_cross_residue_and_records': bounds(),
            'excluded_transfer': 'no two-cell periodic spectral reference, no full-source canonical composite CAR, no free-block/affine Hamiltonian substitution'}
    if digest(RESEARCH) != RESEARCH_SHA:
        raise ValueError('the immutable independent Higgs common-energy calculation changed')
    ignored = {'schema', 'input_sha256', 'local_sha256'}
    scientific = lambda value: {key: item for key, item in value.items() if key not in ignored}
    if scientific(result) != scientific(json.loads(RESEARCH.read_text())):
        raise ValueError('public native Higgs science differs from frozen research')
    result['independent_research_comparison'] = {
        'reference_sha256': RESEARCH_SHA, 'all_scientific_fields_equal': True,
        'ignored_only_provenance_fields': sorted(ignored), 'fresh_calculation_was_replayed': True}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true'); mode.add_argument('--check', action='store_true')
    args = parser.parse_args(); path = HERE/'results.json'
    if args.check:
        saved = json.loads(path.read_text()); validate(saved['input_sha256'])
        for name, expected in saved['local_sha256'].items():
            if digest(HERE/name) != expected: raise ValueError('local custody mismatch: '+name)
    actual = certificate()
    if args.write: path.write_text(json.dumps(actual, indent=2, sort_keys=True)+'\n')
    elif actual != saved: raise ValueError('same-native-H Higgs receipt differs')
    print('PASS: full-native current/EOM ancestry, positive common-energy cross residue and complete interference records')


if __name__ == '__main__':
    main()
