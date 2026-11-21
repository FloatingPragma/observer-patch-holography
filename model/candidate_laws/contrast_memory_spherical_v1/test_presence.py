"""Source, archive, numerical and receipt checks for the two-reading instrument."""
from copy import deepcopy
from fractions import Fraction
from itertools import product

import numpy as np
import pytest
from scipy.linalg import expm

from model.executable_model import carrier
from . import presence as p, presence_certificate as cert, reader, sectors


def test_public_code_and_reversible_reduction():
    # Every dirty eleven-bit word has one decoded status, including uncommitted
    # or malformed public words. Each parent XOR is reversible for dirty targets.
    for word in range(reader.APP):
        decision, triplet = (word >> 6) & 3, word >> 8
        expected = {(1, 7): 1, (1, 0): 0, (2, 0): 2, (3, 0): 4}.get((decision, triplet), 8)
        assert p.status_bits(word >> 6) == expected
    workspace = np.array(list(product(range(16), repeat=3)))
    dirty = workspace.copy()
    p.apply_layer(workspace, [(0, 1, 2)])
    np.testing.assert_array_equal(workspace[:, 2], dirty[:, 2] ^ (dirty[:, 0] | dirty[:, 1]))
    p.apply_layer(workspace, [(0, 1, 2)])
    np.testing.assert_array_equal(workspace, dirty)
    for leaves in (1, 2, 3, 29, 30):
        layers = p.tree_layers(leaves)
        assert sum(map(len, layers)) == leaves-1
        assert len(layers) == (leaves-1).bit_length()


def test_all_source_columns_and_mixed_public_failures():
    sector = sectors.MultisetSector(carrier.cell(), 2)
    flags = p.flags_with_fuel(sector)
    assert flags.shape == (79, 30)
    for placement in ('nominal', 'first-seam mismatch', 'first-seam unavailable', 'all common flip'):
        controls = p.controls_for(placement, 30)
        coarse = p.compile_records(flags, controls, False)
        fine = p.compile_records(flags, controls, True)
        for i, basis_state in enumerate([(0,)*12] + sector.states):
            codes = tuple(reader.gate(int(basis_state[a] == basis_state[b] == 1), 0,
                                      **reader.CONTROLS[control]) >> 6
                          for (a, b), control in zip(carrier.CELL_SEAMS, controls))
            label = (int(29 in codes) | 2*int(2 in codes) | 4*int(3 in codes))
            assert coarse[i] == (label, (0,)*30)
            assert fine[i] == (label, codes)
        if placement == 'nominal':
            assert sorted(map(len, p.partition(coarse).values())) == [30, 49]
            assert sorted(map(len, p.partition(fine).values())) == [1]*30+[49]
        if placement == 'all common flip':
            assert all(label == 1 for label, _ in coarse)
    mixed = ['mismatch', 'unavailable'] + ['nominal']*28
    assert {label for label, _ in p.compile_records(flags, mixed, False)} == {6, 7}


def test_certain_states_and_reference_coherence():
    sector = sectors.MultisetSector(carrier.cell(), 2)
    flags = p.flags_with_fuel(sector)
    active = np.flatnonzero(flags.any(axis=1))
    inactive = np.flatnonzero(~flags.any(axis=1))
    nominal = ['nominal']*30
    for keep, expected in [(False, 1.), (True, 1/30)]:
        groups = p.partition(p.compile_records(flags, nominal, keep))
        # Gram matrix of the actual final apparatus vectors. It multiplies the
        # source matrix elementwise after forgetting the resolved seam labels.
        gram = np.zeros((79, 79))
        for indices in groups.values():
            gram[np.ix_(indices, indices)] = 1
        np.testing.assert_array_equal(gram.diagonal(), np.ones(79))
        np.testing.assert_array_equal(gram[np.ix_(active, inactive)], 0)
        np.testing.assert_array_equal(gram[np.ix_(inactive, inactive)], 1)
        assert gram[np.ix_(active, active)].mean() == expected
        # On sum_e |e>|e> / sqrt(30), the difference to the ideal pure output
        # lives on the correlated 30-dimensional subspace, with this matrix.
        delta = (1 - gram[np.ix_(active, active)])/30
        assert abs(np.sum(abs(np.linalg.eigvalsh(delta)))/2 - (1-expected)) < 1e-14


@pytest.fixture(scope='module')
def campaign():
    return p.campaign()


def test_complete_frozen_replay(campaign):
    p.scientific_gates(campaign)
    p.verify(campaign, p.load_json(p.HERE/'presence_results.json'))
    assert [(r['kappa'], r['first_reader']) for r in campaign['cases']] == [
        (8, 'nominal'), (8, 'first-seam mismatch'), (8, 'first-seam unavailable'),
        (8, 'all common flip'), (0, 'nominal'), (-8, 'nominal')]


@pytest.mark.parametrize('kappa,placement', [
    (8, 'nominal'), (8, 'first-seam mismatch'), (8, 'first-seam unavailable'),
    (8, 'all common flip'), (0, 'nominal'), (-8, 'nominal')])
def test_full_dense_history_and_retained_fuel(campaign, kappa, placement):
    sector = sectors.MultisetSector(carrier.cell(), 2)
    flags = p.flags_with_fuel(sector)
    # A separate dense calculation in the unshifted energy convention. The
    # explicit masks use basis occupations, with fuel first and doubled loads intact.
    source = np.zeros((79, 79))
    source[1:, 1:] = sector.hamiltonian(kappa)
    emitting = np.flatnonzero(flags[:, 0]).item()
    source[0, 0] = source[emitting, emitting]
    pulse = source.copy()
    pulse[emitting, 0] = pulse[0, emitting] = 64
    free = expm(-4j*source)
    state = free @ expm(-1j*np.pi*pulse/128)[:, 0]
    diagonal = np.r_[source[0, 0], kappa*sector.C]
    state *= np.exp(-1j*diagonal/16384)
    masks = [~flags.any(axis=1), *flags.T]
    row = next(r for r in campaign['cases'] if (r['kappa'], r['first_reader']) == (kappa, placement))
    assert abs(abs(state[0])**2 - row['residual_fuel']) < 1e-15

    last = [masks[address.index(29)+1] if code else masks[0] for code, address in campaign['final_labels']]
    # Independently decode both retained alphabets, using the documented bit
    # operations for all basis inputs, including fuel and the doubled sector.
    words = np.where(flags, 29, 1)
    if placement == 'first-seam mismatch':
        words[:, 0] = 2
    elif placement == 'first-seam unavailable':
        words[:, 0] = 3
    elif placement == 'all common flip':
        words = np.where(flags, 1, 29)
    codes = (np.any(words == 29, axis=1).astype(int)
             + 2*np.any(words == 2, axis=1) + 4*np.any(words == 3, axis=1))
    for detailed, labels, actual in [
            (False, [(code, ()) for code in row['aggregate_codes']], row['aggregate_joint']),
            (True, row['resolved_first_labels'], row['resolved_joint'])]:
        for (code, address), law in zip(labels, actual, strict=True):
            mask = codes == code
            if detailed:
                mask &= np.all(words == address, axis=1)
            after = free @ (mask*state)
            expected = [np.vdot(after[mask], after[mask]).real for mask in last]
            np.testing.assert_allclose(law, expected, rtol=1e-9, atol=2e-14)


def test_integer_generators_match_the_primary_emitter_and_free_law():
    free, pulse, seams, emitting = cert.generators()
    pairs = sectors.PairSector(carrier.cell())
    h = pairs.hamiltonian(8).toarray()
    initial = pairs.index(0, 1)
    assert emitting == initial+1
    np.testing.assert_array_equal(seams, pairs.seam_states+1)
    expected = np.zeros((67, 67))
    expected[1:, 1:] = 2*(h - h[initial, initial]*np.eye(66))
    for rows in (free, pulse):
        if rows is pulse:
            expected[0, emitting] = expected[emitting, 0] = 2*sectors.G_EMIT
        literal = np.zeros_like(expected)
        for i, entries in enumerate(rows):
            for j, value in entries.items():
                literal[i, j] = value
        np.testing.assert_array_equal(literal, expected)


def test_rational_witness_exceeds_full_finite_duration_error(campaign):
    result = cert.event_certificate()
    assert result['always_on_tv_lower'] > Fraction(817, 10000)
    assert result['each_probability_error'] < Fraction(3, 10**7)
    # Identify the witness by its copied words, independently of label sorting.
    final = [(code, tuple(address)) for code, address in campaign['final_labels']]
    emitting_label = (1, (29,) + (1,)*29)
    column = final.index(emitting_label)
    row = campaign['cases'][0]
    coarse = row['aggregate_joint'][row['aggregate_codes'].index(1)][column]
    fine = sum(law[column] for (code, _), law in zip(row['resolved_first_labels'], row['resolved_joint'])
               if code == 1)
    for name, value in [('aggregate_event', coarse), ('resolved_event', fine)]:
        assert abs(Fraction(value)-result[name]) <= result['each_probability_error']


def test_gate_rejects_normalized_but_false_joint_laws(campaign):
    stored = p.load_json(p.HERE/'presence_results.json')
    for mutation in ('relabel', 'alter', 'prune', 'negative', 'bool', 'calibration'):
        bad = deepcopy(stored)
        row = bad['cases'][0]
        if mutation == 'relabel':
            row['resolved_first_labels'][0][1][0] ^= 1
        elif mutation == 'alter':
            row['aggregate_joint'][0].reverse()
        elif mutation == 'prune':
            row['resolved_joint'].pop()
        elif mutation == 'negative':
            row['aggregate_joint'][0][0] = -1e-16
        elif mutation == 'bool':
            row['kappa'] = True
        else:
            row['public_joint_tv'] = 0.
        with pytest.raises((ValueError, IndexError)):
            p.verify(campaign, bad)


def test_replay_checks_scientific_gates_on_the_supplied_receipt(campaign):
    stored = p.load_json(p.HERE/'presence_results.json')
    bad = deepcopy(stored)
    row = bad['cases'][0]
    for key in ('aggregate_joint', 'resolved_joint'):
        row[key] = (np.array(row[key])*(1+5e-9)).tolist()
    with pytest.raises(ValueError, match='incomplete joint law'):
        p.verify(campaign, bad)
    bad = deepcopy(stored)
    bad['cases'][0]['public_joint_tv'] *= 1+5e-9
    with pytest.raises(ValueError, match='inconsistent public_joint_tv'):
        p.verify(campaign, bad)
    bad = deepcopy(stored)
    bad['cases'].pop()
    with pytest.raises(ValueError, match='incomplete campaign'):
        p.scientific_gates(bad)


def test_dirty_word_must_define_a_hermitian_pulse(monkeypatch):
    native_stages = reader.stages
    def corrupt(flag, **kwargs):
        stages = native_stages(flag, **kwargs)
        stages[0][-1] = stages[0][-2]
        return stages
    monkeypatch.setattr(reader, 'stages', corrupt)
    with pytest.raises(ValueError, match='involution'):
        p.compile_records(np.array([[False], [True]]), ['nominal'], False)
