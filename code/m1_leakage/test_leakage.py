"""Full-channel, exact-basis, hostile-input and artifact-custody controls."""

import copy
import fnmatch
import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest
from scipy.linalg import expm

from . import build, correction, correction_check, interface_check, interfaces
from . import noise, noise_check, scaling, scaling_check, composition_check, verify
from .format import pack, unpack


@pytest.fixture(scope='module')
def evidence():
    return build.candidate()


def test_receipt_and_recursive_parent_custody():
    verify.verify()


def test_rebuild_replays_all_science(evidence):
    verify.verify_evidence(evidence)
    stored = verify.load(verify.HERE/'receipt.json')['evidence']
    assert evidence['scaling'] == stored['scaling']
    assert evidence['correction']['basis_census'] == stored['correction']['basis_census']


@pytest.mark.parametrize('d', [2, 4])
def test_native_reduction_is_tp_on_code_complement_and_coherence(d):
    groups = [[pack(k) for k in group] for group in interfaces.reduction(d)]
    interface_check.check_groups(groups, d)
    # A mathematical dilation retains private branch information, but does
    # not grant the source a reversible operation on the full input algebra.
    ks = sum(interfaces.reduction(d), [])
    v = np.concatenate(ks, axis=0)
    np.testing.assert_allclose(v.conj().T@v, np.eye(6), atol=1e-14)
    coherence = np.zeros((6, 6))
    coherence[0, d] = 1
    assert np.linalg.norm(sum(k@coherence@k.conj().T for k in ks)) == 0


@pytest.mark.parametrize('d', [2, 4])
@pytest.mark.parametrize('mutation', ['omit_failure', 'coherent_merge', 'wrong_blank', 'code_phase'])
def test_physical_reduction_mutants_fail_without_hashes(d, mutation):
    groups = interfaces.reduction(d)
    if mutation == 'omit_failure':
        groups[1] = [np.zeros_like(k) for k in groups[1]]
    elif mutation == 'coherent_merge':
        groups[1][0] += groups[1][1]
        groups[1][1] *= 0
    elif mutation == 'wrong_blank':
        for k in groups[1]:
            k[1] = k[0]
            k[0] = 0
    else:
        # This preserves all computational input probabilities, but changes
        # an input coherence. Every public label is still present.
        groups[0][0][1, 1] *= -1
    if mutation in ('wrong_blank', 'code_phase'):
        np.testing.assert_allclose(sum(k.conj().T@k for group in groups for k in group), np.eye(6))
    with pytest.raises(ValueError):
        interface_check.check_groups([[pack(k) for k in group] for group in groups], d)


def test_pair_channel_mutation_is_tp_but_changes_the_other_carrier(evidence):
    row = copy.deepcopy(evidence['interfaces']['pair_gates'][0])
    # After a first-carrier leak, illegally flip the surviving second qubit.
    flip = np.kron(np.eye(2), np.array([[0, 1], [1, 0]]))
    branch = row['branches'][2]
    branch['maps'] = [pack(flip@unpack(k, (4, 36))) for k in branch['maps']]
    all_maps = [unpack(k, (4, 36)) for b in row['branches'] for k in b['maps']]
    np.testing.assert_allclose(sum(k.conj().T@k for k in all_maps), np.eye(36), atol=1e-14)
    with pytest.raises(ValueError, match='complete two-carrier'):
        interface_check.check_pair(row, 'cz')


@pytest.mark.parametrize('value', [True, 2., 3, None])
def test_native_code_dimension_is_strict(value):
    with pytest.raises(ValueError):
        interfaces.reduction(value)


@pytest.mark.parametrize('mutation', ['shape', 'duplicate', 'order', 'nan', 'inf', 'boolean', 'zero', 'index', 'extra'])
def test_sparse_matrix_rejects_trash(mutation):
    row = pack(np.eye(2))
    if mutation == 'shape':
        row['shape'] = [2., 2]
    elif mutation == 'duplicate':
        row['entries'].append(row['entries'][0])
    elif mutation == 'order':
        row['entries'].reverse()
    elif mutation == 'nan':
        row['entries'][0][2] = float('nan')
    elif mutation == 'inf':
        row['entries'][0][2] = float('inf')
    elif mutation == 'boolean':
        row['entries'][0][2] = True
    elif mutation == 'zero':
        row['entries'][0][2] = 0.
    elif mutation == 'index':
        row['entries'][0][0] = -1
    else:
        row['extra'] = 1
    with pytest.raises(ValueError):
        unpack(row, (2, 2))


def test_noise_on_a_noncommuting_pulse_is_not_an_endpoint_surrogate():
    _, t, h, b, jumps = noise.specification('driven_return')
    actual = expm(t*noise.generator(h+b, jumps))
    endpoint = expm(t*noise.generator(b, jumps))@expm(t*noise.generator(h, []))
    assert np.linalg.norm(actual-endpoint) > .01
    independent = noise_check.evolution(h+b, jumps, t)
    np.testing.assert_allclose(actual, independent, atol=2e-13)


@pytest.mark.parametrize('mutation', ['endpoint', 'no_return', 'no_leak', 'wrong_flag', 'false_bound'])
def test_dissipative_channel_mutants_fail_independent_replay(evidence, mutation):
    rows = copy.deepcopy(evidence['noise']['channels'])
    row = rows[2]
    d, t, h, b, jumps = noise.specification('driven_return')
    if mutation == 'endpoint':
        wrong = expm(t*noise.generator(b, jumps))@expm(t*noise.generator(h, []))
        row['channel_choi'] = pack(noise.choi(wrong, 6))
    elif mutation == 'no_return':
        wrong = expm(t*noise.generator(h+b, jumps[:1]))
        row['channel_choi'] = pack(noise.choi(wrong, 6))
    elif mutation == 'no_leak':
        wrong = expm(t*noise.generator(h+b, jumps[1:]))
        row['channel_choi'] = pack(noise.choi(wrong, 6))
    elif mutation == 'wrong_flag':
        j = unpack(row['flagged_choi'], (8, 8))
        # Relabel success and failure; still CPTP, but changes public meaning.
        swap = np.kron(np.eye(2), np.roll(np.eye(4), 2, axis=0))
        row['flagged_choi'] = pack(swap@j@swap.T)
    else:
        row['diamond_upper'] = 0.
    if mutation in ('endpoint', 'no_return', 'no_leak'):
        noise_check.cptp(unpack(row['channel_choi'], (36, 36)), 6, 6)
    with pytest.raises(ValueError):
        noise_check.verify_channels(rows)


def test_coherent_public_outcomes_are_detected_even_when_channel_is_tp(evidence):
    rows = copy.deepcopy(evidence['noise']['channels'])
    d, t, h, b, _ = noise.specification('coherent_leak')
    u = expm(-1j*t*(h+b))
    v = np.zeros((4, 2), complex)
    v[:2] = u[:2, :2]
    v[2] = u[4, :2]
    np.testing.assert_allclose(v.conj().T@v, np.eye(2), atol=1e-14)
    vector = v.T.reshape(-1)/np.sqrt(2)
    coherent = np.outer(vector, vector.conj())
    noise_check.cptp(coherent, 2, 4)
    rows[4]['flagged_choi'] = pack(coherent)
    with pytest.raises(ValueError, match='unconditional flagged'):
        noise_check.verify_channels(rows)


@pytest.mark.parametrize('gamma', [.01, .2, .7])
@pytest.mark.parametrize('weight', [1e-5, .1, .8])
def test_amplitude_damping_has_no_positive_identity_component(gamma, weight):
    k0 = np.diag([1., math.sqrt(1-gamma)])
    k1 = np.array([[0., math.sqrt(gamma)], [0., 0.]])
    j = sum(np.outer(k.T.reshape(-1), k.T.reshape(-1)) for k in (k0, k1))
    identity = np.array([1., 0, 0, 1.])
    witness = np.array([math.sqrt(1-gamma), 0, 0, -1.])
    actual = witness@(j-weight*np.outer(identity, identity))@witness
    assert actual < 0
    assert actual == pytest.approx(-weight*(1-math.sqrt(1-gamma))**2, abs=3e-16)
    # The natural common-ancilla dilation has an amplitude of square-root
    # order even though the channel error has first-order damping strength.
    v = np.vstack((k0, k1))
    ideal = np.vstack((np.eye(2), np.zeros((2, 2))))
    amplitude_squared = np.linalg.norm(v-ideal, 2)**2
    assert amplitude_squared == pytest.approx(2-2*math.sqrt(1-gamma))
    assert gamma <= amplitude_squared <= 2*gamma


def test_common_bath_cannot_be_replaced_by_fresh_baths(evidence):
    row = copy.deepcopy(evidence['noise']['bath'])
    row['shared_choi'] = copy.deepcopy(row['reset_choi'])
    with pytest.raises(ValueError, match='same bath retained'):
        noise_check.verify_bath(row)


def test_exact_decoder_certificate_and_two_site_failure(evidence):
    row = evidence['correction']
    correction_check.verify(row)
    assert row['basis_census']['branch_identities'] == 26880
    assert row['basis_census']['scalar_numerator_histogram'] == {'4': 126, '0': 26712, '-4': 42}
    bell = np.array([1, 0, 0, 1])/np.sqrt(2)
    assert abs(bell@unpack(row['two_site_phase_choi'], (4, 4))@bell) < 1e-14


def test_wrong_decoder_fails_the_linear_basis_identity():
    code, maps = correction_check.integer_code()
    maps[1, 1] *= -1
    with pytest.raises(ValueError, match='preserves the logical qubit'):
        correction_check.exact_basis_check(code, maps)


@pytest.mark.parametrize('site', [0, 3, 6])
def test_arbitrary_coherent_full_interface_error_with_a_reference(site):
    rng = np.random.default_rng(821+site)
    arbitrary, _ = np.linalg.qr(rng.normal(size=(6, 2))+1j*rng.normal(size=(6, 2)))
    code, maps = correction_check.integer_code()
    code, maps = code/np.sqrt(8), maps/np.sqrt(8)
    choi = np.zeros((4, 4), complex)
    for _, k in correction_check.reduced_operators([arbitrary]):
        for branch in maps@correction_check.local_action(k, site, code):
            vector = branch.T.reshape(-1)/np.sqrt(2)
            choi += np.outer(vector, vector.conj())
    bell = np.array([1, 0, 0, 1])/np.sqrt(2)
    np.testing.assert_allclose(choi, np.outer(bell, bell), atol=1e-13)


@pytest.mark.parametrize('mutation', ['census', 'site', 'drop_branch', 'false_identity', 'accounting', 'omit_case'])
def test_correction_and_accounting_mutants_fail(evidence, mutation):
    row = copy.deepcopy(evidence['correction'])
    if mutation == 'census':
        row['basis_census']['scalar_numerator_histogram']['0'] += 1
    elif mutation == 'site':
        row['cases'][0]['site'] = True
    elif mutation == 'drop_branch':
        row['cases'][0]['noise_operators'].pop()
    elif mutation == 'false_identity':
        row['two_site_phase_choi'] = row['cases'][0]['recovered_choi']
    elif mutation == 'accounting':
        row['cases'][0]['corrected_accounting'] = row['cases'][0]['product_abort_accounting']
    else:
        row['cases'].pop()
    with pytest.raises(ValueError):
        correction_check.verify(row)


def test_rejecting_any_raw_leak_discards_a_correctable_experiment(evidence):
    row = evidence['correction']['cases'][0]
    assert row['raw_carrier_acceptance'] == 0
    assert row['product_abort_accounting'] == 1
    assert row['corrected_accounting'] < .5


@pytest.mark.parametrize('mutation', ['amplitude', 'old_levels', 'decreasing_error', 'insufficient_levels', 'lifetime', 'boolean'])
def test_false_asymptotic_claims_are_rejected(evidence, mutation):
    row = copy.deepcopy(evidence['scaling'])
    if mutation == 'amplitude':
        row['fixed_rate'][0]['amplitude_decay_denominator'] = 1
    elif mutation == 'old_levels':
        row['fixed_rate'][0]['selected_level'] = 4
    elif mutation == 'decreasing_error':
        row['growing'][0]['contraction_bits'] += 1
    elif mutation == 'insufficient_levels':
        row['growing'][-1]['levels'] -= 1
    elif mutation == 'lifetime':
        row['resources']['physical_storage_volume'] = row['resources']['diagnostic_bits']
    else:
        row['accounting_q_degree'] = True
    with pytest.raises(ValueError):
        scaling_check.verify(row)


@pytest.mark.parametrize('q,bits', [(True, 2), (1, 2), (2., 2), (2, True), (2, 0), (2, float('inf'))])
def test_level_chooser_rejects_trash(q, bits):
    with pytest.raises(ValueError):
        scaling.levels(q, bits)


def test_fixed_error_requires_growing_redundancy(evidence):
    rows = [r for r in evidence['scaling']['growing'] if r['contraction_bits'] == 3]
    assert rows[-1]['levels'] > rows[0]['levels']
    assert all((1 << r['suppression_bits']) >= (2*r['q'])**16 for r in rows)
    errors = [r['half_trace_distance'] for r in evidence['noise']['calibration']]
    np.testing.assert_allclose(errors, [errors[0]]*4, atol=1e-14)
    assert errors[0] > .01


@pytest.mark.parametrize('section', ['interfaces', 'noise', 'correction', 'scaling', 'composition'])
def test_evidence_sections_are_required(evidence, section):
    row = copy.deepcopy(evidence)
    del row[section]
    with pytest.raises(ValueError):
        verify.verify_evidence(row)


def test_strict_json_custody_and_nonfinite_inputs(tmp_path):
    for raw in ('{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}'):
        path = tmp_path/'bad.json'
        path.write_text(raw, encoding='utf-8')
        with pytest.raises(ValueError):
            verify.load(path)
    packet = verify.load(verify.HERE/'receipt.json')
    for section in ('sources', 'claims'):
        bad = copy.deepcopy(packet)
        bad[section][next(iter(bad[section]))] = '0'*64
        path = tmp_path/(section+'.json')
        path.write_text(json.dumps(bad), encoding='utf-8')
        with pytest.raises(ValueError, match='custody'):
            verify.verify(path)


def test_optimized_replay_with_all_producers_disabled():
    script = '''
import sys
class Block:
    def find_spec(self, fullname, path=None, target=None):
        if fullname in ('m1_leakage.interfaces', 'm1_leakage.noise', 'm1_leakage.correction', 'm1_leakage.scaling', 'm1_leakage.composition', 'm1_leakage.build'):
            raise RuntimeError('producer called')
sys.meta_path.insert(0, Block())
from m1_leakage import verify
packet = verify.load(verify.HERE/'receipt.json')
verify.verify_evidence(packet['evidence'])
packet['evidence']['scaling']['fixed_rate'][0]['selected_level'] = 1
try:
    verify.verify_scaling(packet['evidence']['scaling'])
except ValueError:
    pass
else:
    raise RuntimeError('optimized verification accepted false levels')
'''
    result = subprocess.run([sys.executable, '-O', '-c', script], cwd=verify.ROOT,
                            env=dict(os.environ, PYTHONPATH=str(verify.ROOT/'code')),
                            text=True, capture_output=True, timeout=120)
    assert result.returncode == 0, result.stdout+result.stderr


@pytest.mark.parametrize('mutation', ['stochastic_tail', 'multiplicity', 'boolean', 'amplitude', 'old_levels',
                                     'old_volume', 'discard_lifetimes', 'synthesis', 'false_scope', 'nonfinite'])
def test_combined_noisy_archive_claims_reject_false_shortcuts(evidence, mutation):
    row = copy.deepcopy(evidence['composition'])
    if mutation == 'stochastic_tail':
        row['tails'][1]['coherent_bad_majorant'] = [1, 64]
    elif mutation == 'multiplicity':
        row['tails'][1]['coefficient_by_fault_count'][-1] = 2
    elif mutation == 'boolean':
        row['tails'][1]['coefficient_by_fault_count'][-1] = True
    elif mutation == 'amplitude':
        row['amplitude_decay_denominator'] = 1
    elif mutation == 'old_levels':
        row['quantum_level'] = 5
    elif mutation == 'old_volume':
        row['active_locations'] = [4, 6]
    elif mutation == 'discard_lifetimes':
        row['resources']['physical_storage_volume'] = [5, 16]
    elif mutation == 'synthesis':
        row['synthesis_per_gate_decay'] = 16
    elif mutation == 'false_scope':
        row['model'] = 'arbitrary noise at constant per-operation error'
    else:
        row['archive_tail'][0]['log_tail_bound'] = float('inf')
    with pytest.raises(ValueError):
        composition_check.verify(row)


def test_combined_tail_is_asymptotic_and_keeps_the_amplitude_loss(evidence):
    row = evidence['composition']
    composition_check.verify(row)
    # Normalized illustrative constants give no useful early bound. They
    # must not be advertised as a physical onset or a machine threshold.
    assert row['archive_tail'][0]['log_tail_bound'] > 0
    final = row['archive_tail'][-1]
    log_q = final['log2_q']*math.log(2)
    assert final['log_tail_bound'] < -30*log_q
    assert final['log_tail_bound'] > -log_q**2/(2*math.log(2))


def test_every_merged_parent_claim_is_replayed():
    assert verify.parent.CLAIMS == ('OPH-SOURCE-NOISY-CONTROL-HISTORIES', 'OPH-SOURCE-LIVE-RECORD-PROTECTION')
    assert 'code/m1_noisy_records/receipt.json' in verify.SOURCES
    assert 'code/m1_fixed_noise/receipt.json' in verify.SOURCES


@pytest.mark.parametrize('kind', ['erasure', 'leak_damping', 'coherent_leak'])
def test_corrected_leakage_preserves_actual_adaptive_histories_and_aborts(kind):
    from m1_noisy_records import circuits, independent
    code, maps = correction.integer_maps()
    code, maps = code/np.sqrt(8), maps/np.sqrt(8)
    ks = correction.corrected(maps, code, correction.noise_operators(kind), 3,
                              sum(interfaces.reduction(2), []))
    angles = (.37, -.61, .43)
    v = circuits.execute(circuits.tape(*angles))
    output = np.zeros((128, 128), complex)
    for k in ks:
        if np.linalg.norm(k) > 1e-14:
            output += independent.choi_from_isometry(v@np.kron(np.eye(2), k))
    np.testing.assert_allclose(output, independent.reference(angles), atol=2e-13)
    # This is a finite composition of the exact local decoder and actual
    # causal instrument, not a simulation of the six-level computation code.
    # In output-input ordering, every named history retains nonzero weight.
    for r0 in (0, 1):
        for r1 in (0, 1):
            begin = 4*(4*r0+8*r1+16*r0*r1)
            assert np.trace(output[begin:begin+16, begin:begin+16]).real > .01


@pytest.mark.parametrize('d', [2, 4])
def test_private_kraus_basis_changes_preserve_public_instruments(d):
    groups = interfaces.reduction(d)
    count = 6-d
    basis = np.exp(2j*np.pi*np.outer(np.arange(count), np.arange(count))/count)/math.sqrt(count)
    groups[1] = [sum(basis[i, j]*groups[1][j] for j in range(count)) for i in range(count)]
    interface_check.check_groups([[pack(k) for k in group] for group in groups], d)


def test_cross_public_kraus_rotation_fails_even_with_identical_unflagged_channel():
    original = interfaces.reduction(2)
    groups = copy.deepcopy(original)
    a, b = groups[0][0].copy(), groups[1][0].copy()
    groups[0][0], groups[1][0] = (a+b)/math.sqrt(2), (a-b)/math.sqrt(2)
    for i in range(6):
        for j in range(6):
            np.testing.assert_allclose(interface_check.image(sum(groups, []), i, j),
                                       interface_check.image(sum(original, []), i, j), atol=1e-14)
    with pytest.raises(ValueError, match='complete flagged reduction'):
        interface_check.check_groups([[pack(k) for k in group] for group in groups], 2)


@pytest.mark.parametrize('mutation', ['overlap', 'boolean', 'missing_region', 'wrong_threshold',
                                     'naive_sum', 'positive_mixture', 'changed_operator', 'missing_case'])
def test_coherent_union_countercontrols_reject_invalid_shortcuts(evidence, mutation):
    rows = copy.deepcopy(evidence['composition']['coherent_unions'])
    row = rows[0]
    if mutation == 'overlap':
        row['regions'][1] = [1, 2]
    elif mutation == 'boolean':
        row['regions'][0][0] = False
    elif mutation == 'missing_region':
        row['regions'].pop()
    elif mutation == 'wrong_threshold':
        row['required_faults'][0] = 1
    elif mutation == 'naive_sum':
        row['product_majorant'] = row['naive_sum']
    elif mutation == 'positive_mixture':
        row['joint_bad_norm'] = row['naive_sum']
    elif mutation == 'changed_operator':
        row['bad_operator'] = pack(np.eye(2)*row['local_bad_norm'])
    else:
        rows.pop()
    with pytest.raises(ValueError):
        composition_check.verify_unions(rows)


def test_joint_coherent_fault_norm_needs_intersection_terms(evidence):
    rows = evidence['composition']['coherent_unions']
    composition_check.verify_unions(rows)
    assert all(row['joint_bad_norm'] > row['naive_sum'] for row in rows)
    assert all(row['joint_bad_norm'] <= row['product_majorant'] for row in rows)


def test_disjoint_sparse_rectangles_require_four_faults_at_two_levels():
    # Finite control of the combinatorial induction, not an execution of the
    # imported computation-code gadgets. Three children at each of two levels.
    bad_weights = []
    for mask in range(1 << 9):
        bad_children = sum(((mask >> (3*i)) & 7).bit_count() >= 2 for i in range(3))
        if bad_children >= 2:
            bad_weights.append(mask.bit_count())
    assert min(bad_weights) == 4
    assert max(bad_weights) == 9
