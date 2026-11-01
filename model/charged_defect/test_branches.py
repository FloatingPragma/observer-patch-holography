"""Independent occupation/CAR and exact-algebra checks of the charged lift."""
from itertools import product

import numpy as np
import pytest

from ..exchange_sector import statistics
from . import branches as b


def car_move(mask, origin, destination):
    """Apply c_destination^dagger c_origin on an ordered 24-mode wedge."""
    assert mask >> origin & 1
    sign = (-1)**(mask & ((1 << origin)-1)).bit_count()
    mask ^= 1 << origin
    assert not mask >> destination & 1
    sign *= (-1)**(mask & ((1 << destination)-1)).bit_count()
    return mask | (1 << destination), sign


@pytest.mark.parametrize('law', ('C', 'E'))
def test_every_internal_column_against_independent_spinful_CAR(law):
    count = 0
    for seam in statistics.SEAMS:
        left, right = seam
        for coin in (0, 1):
            oracle = [[] for _ in range(4 if law == 'C' else 8)]
            for x, pair in enumerate(b.PAIRS):
                l, r = int(left in pair), int(right in pair)
                current = ((l-r)-(2*coin-1)*((l+r) % 2))//2
                raw = {0: 0, 1: 1, -1: 2}[current]
                if law == 'E':
                    raw += 4*((l+r) % 2)
                for s, t in product(range(2), repeat=2):
                    mask = (1 << (2*pair[0]+s)) | (1 << (2*pair[1]+t))
                    sign = 1
                    if current:
                        origin, destination = (left, right) if current > 0 else (right, left)
                        spin = s if pair[0] == origin else t
                        mask, sign = car_move(mask, 2*origin+spin, 2*destination+spin)
                    modes = tuple(j for j in range(24) if mask >> j & 1)
                    output_pair = tuple(j//2 for j in modes)
                    output = 4*b.INDEX[output_pair]+2*(modes[0] % 2)+modes[1] % 2
                    oracle[raw].append((4*x+2*s+t, output, sign))
            for raw, expected in enumerate(oracle):
                actual = b.branch(law, seam, coin, raw)
                assert list(zip(actual.inputs, actual.outputs, actual.signs, strict=True)) == expected
                count += len(expected)
    assert count == 15840


@pytest.mark.parametrize('law', ('C', 'E'))
def test_complete_singlet_triplet_intertwiner_and_full_commutant(law):
    cert = b.decomposition_certificate(law)
    assert cert['all_30_seams_complete'] and cert['all_30_seams_unital']
    assert cert['triplet_noiseless_factor_dimension'] == 3
    proof = b.commutant_certificate(law)
    rows = statistics.branches(law, -1)
    signatures = [tuple(i in row.inputs for row in rows) for i in range(66)]
    assert len(set(signatures)) == 66
    parity = {0: 1}
    for x, y, number, sign in proof['transition_tree']:
        assert x in parity and y not in parity
        row = rows[number]
        assert (x, y, sign) in list(zip(row.inputs, row.outputs, row.signs)) or (
            y, x, sign) in list(zip(row.inputs, row.outputs, row.signs))
        parity[y] = parity[x]*sign
    x, y, number, sign = proof['odd_chord']
    row = rows[number]
    assert (x, y, sign) in list(zip(row.inputs, row.outputs, row.signs)) or (
        y, x, sign) in list(zip(row.inputs, row.outputs, row.signs))
    assert parity[x]*sign != parity[y]
    assert len(parity) == 66
    assert proof['commutant_complex_dimension'] == 10


@pytest.mark.parametrize('law', ('C', 'E'))
def test_arbitrary_triplet_reference_state_survives_complete_record_instrument(law):
    # Include source/internal/reference correlations.  Only the COMPLETE
    # output, not each postselected normalized branch, must preserve the
    # internal/reference marginal for such correlated inputs.
    rng = np.random.default_rng(1139)
    state = rng.normal(size=(66, 3, 2))+1j*rng.normal(size=(66, 3, 2))
    state /= np.linalg.norm(state)
    triplet = b.COUPLING[:, 1:]/np.sqrt(np.array([1., 2., 1.]))
    embedded = np.einsum('aj,xjr->xar', triplet, state).reshape(264, 2)
    expected = np.einsum('xar,xbs->arbs', state, state.conj()).reshape(6, 6)
    actual = np.zeros((6, 6), complex)
    probability = 0.
    for coin in (0, 1):
        for raw in range(4 if law == 'C' else 8):
            out = b.branch(law, (0, 5), coin, raw).matrix()@embedded/np.sqrt(2)
            internal = np.einsum('aj,xar->xjr', triplet.conj(), out.reshape(66, 4, 2))
            probability += float(np.vdot(out, out).real)
            actual += np.einsum('xar,xbs->arbs', internal, internal.conj()).reshape(6, 6)
    np.testing.assert_allclose(actual, expected, atol=5e-16, rtol=0)
    assert probability == pytest.approx(1, abs=5e-16)


def test_singlet_triplet_coherence_is_not_a_four_state_passive_memory():
    singlet, triplet = b.COUPLING[:, 0], b.COUPLING[:, 1]
    cross = np.outer(singlet, triplet)
    assert np.any(cross@b.SWAP != b.SWAP@cross)
    # The negative exchange affects the two sectors with opposite signs.
    np.testing.assert_array_equal(-b.SWAP@singlet, singlet)
    np.testing.assert_array_equal(-b.SWAP@triplet, -triplet)


@pytest.mark.parametrize('args', [('X', (0, 5), 0, 0), ('C', (0, 5), True, 0),
                                 ('C', (0, 5), 0, 4), ('E', (5, 0), 0, 0)])
def test_invalid_native_branch_labels_reject(args):
    with pytest.raises(ValueError):
        b.branch(*args)
