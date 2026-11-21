"""Full two-constituent coordinate check of the one-rail clock and its local reader."""
import math

import numpy as np
import pytest

from . import formed_clock_readout as measured


def _operators(size):
    # Basis index 8*z+s. Only the third spatial axis has multiple positions.
    equal = np.full(8, 1 / math.sqrt(8))
    bright = np.outer(equal, equal)
    reflection = 2 * bright - np.eye(8)
    phi = math.pi / 16
    contact = np.exp(-.5j * phi) * bright - np.exp(.5j * phi) * (np.eye(8) - bright)
    to = np.empty(8 * size, dtype=int)
    for z in range(size):
        for slot in range(8):
            shift = 1 if slot == 4 else -1 if slot == 5 else 0
            to[8 * z + slot] = 8 * ((z + shift) % size) + (slot ^ 1)
    return equal, reflection, contact, np.argsort(to)


def _coin_pair(wave, ordinary, contact):
    # Apply the independent ordinary coins, then replace all same-site blocks.
    out = ordinary @ wave @ ordinary.T
    for z in range(len(wave) // 8):
        slots = np.arange(8 * z, 8 * z + 8)
        out[np.ix_(slots, slots)] = wave[np.ix_(slots, slots)]
        out[slots, slots] = contact @ wave[slots, slots]
    return out


def _root_rotation(words, size, wire_length, equal, angle):
    # Within the vacuum/one-marker history, a root couples the head-on pair
    # to the equal matched pair with that site's wire bit set and age zero.
    vacuum = words[0, -1]
    for site in range(size):
        ids = np.arange(8 * site, 8 * site + 8)
        root = (1 << (site * wire_length), 0)
        x = vacuum[ids[5], ids[4]]
        if x == 0 and root not in words:
            continue
        marked = words.setdefault(root, np.zeros_like(vacuum))
        y = equal @ marked[ids, ids]
        vacuum[ids[5], ids[4]] = math.cos(angle) * x - math.sin(angle) * y
        marked[ids, ids] += (math.sin(angle) * x + (math.cos(angle) - 1) * y) * equal


def _read_train(size, far, angle, *, reverse=False, competing=False):
    equal, coin, contact, route = _operators(size)
    ordinary = np.kron(np.eye(size), coin)
    a = np.zeros((size, 8))
    b = np.zeros_like(a)
    for distance, power in ((1, 1 - far), (2, far)):
        a[0 if competing and distance == 2 else -distance] = math.sqrt(power) * coin[:, 4]
        b[distance] = math.sqrt(power) * coin[:, 5]
    vacuum = np.outer(b, a).astype(complex) if reverse else np.outer(a, b).astype(complex)
    words = {(0, -1): vacuum}
    wire_length = 257
    _root_rotation(words, size, wire_length, equal, angle)
    early = sum(float(np.vdot(v, v).real) for (mask, _), v in words.items() if mask)
    # The early outcome terminates that trial. Advance both endpoints of
    # every unmarked amplitude to the next vertex before reading its wire.
    words = {(0, -1): _coin_pair(vacuum, ordinary, contact)[np.ix_(route, route)]}
    _root_rotation(words, size, wire_length, equal, angle)
    unmarked = float(np.vdot(words[0, -1], words[0, -1]).real)
    elsewhere = sum(float(np.vdot(v, v).real) for (mask, _), v in words.items() if mask not in (0, 1))
    pair = words.get((1, 0), np.zeros_like(vacuum))
    weight = float(np.vdot(pair, pair).real)
    np.testing.assert_allclose(early + elsewhere + unmarked + weight, 1, atol=2e-13, rtol=0)
    mask, age = 1, 0
    # Each physical wire exceeds the run, so shifting its one occupied bit
    # never brings it back to the root. All eight slot pairs retain both tags.
    local_words = []
    ids = np.arange(8 * size)
    for appointment in range(1, 256):
        theta = math.pi / 2 if age >= 128 else math.pi * (1 - math.cos(math.pi * (age + 1) / 128)) / 4
        coined = _coin_pair(pair, ordinary, contact)
        routed = coined[np.ix_(route, route)]
        # Spectral projectors of the literal two-endpoint routing involution.
        pair = ((coined + routed) + np.exp(-2j * theta) * (coined - routed)) / 2
        mask = ((mask << 1) % (1 << wire_length)) | (mask >> (wire_length - 1))
        age = (age + 1) % 8192
        if appointment < 128:
            continue
        assert mask == 1 << appointment and age == appointment
        local = pair[ids, ids].reshape(size, 8)
        np.testing.assert_allclose(np.vdot(local, local).real, weight, atol=2e-12, rtol=0)
        # At one named address the odd mode is a two-tag swap interferometer.
        inc = np.sum(abs(local)**2, axis=1)
        odd = abs(local[:, 4] - local[:, 5])**2 / 2
        local_words.append(np.stack((weight - inc, inc - odd, odd), axis=-1))
    return np.array((early, elsewhere, unmarked)), np.asarray(local_words)


@pytest.mark.parametrize('size', (5, 6))
def test_two_tag_wire_words_equal_every_local_clock_probability(size):
    for far in (0., .125, .25, .5):
        refused, local = _read_train(size, far, math.pi / 4)
        q = (1 - far)**2 / 2
        np.testing.assert_allclose(refused, (0, 0, 1 - q), atol=2e-13, rtol=0)
        probes = {'ordinary': local[:, 1], 'reflected': local[:, -1], 'translated': local.mean(axis=1)}
        for menu, probe in probes.items():
            expected = np.column_stack((np.tile(refused, (128, 1)), probe))
            np.testing.assert_allclose(measured.native_words(size, far, menu), expected, atol=2e-12, rtol=0)
        off_refused, off_local = _read_train(size, far, 0.)
        off = np.column_stack((np.tile(off_refused, (128, 1)), off_local[:, 1]))
        np.testing.assert_allclose(measured.native_words(size, far, 'vertex_off'), off, atol=2e-12, rtol=0)
        for menu in ('reversed_collision', 'competing_birth'):
            r, site = _read_train(size, far, math.pi / 4, reverse=menu == 'reversed_collision',
                                  competing=menu == 'competing_birth')
            booked = (0, 0, 1) if menu == 'reversed_collision' else (0, far**2 / 2, 1 - q - far**2 / 2)
            np.testing.assert_allclose(r, booked, atol=2e-13, rtol=0)
            np.testing.assert_allclose(measured.native_words(size, far, menu),
                                       np.column_stack((np.tile(r, (128, 1)), site[:, 1])), atol=2e-12, rtol=0)
        np.testing.assert_allclose(local[:, 1, 2], q * measured.predicted_intensity(size), atol=2e-12, rtol=0)
