"""Independent blank-seam reconstruction, time comparator and integer shape checks."""
import gzip
from itertools import combinations
import json

import numpy as np
from scipy import linalg

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import cell
from .apparatus import ARMS, NATIVE_LINE, STARTS, STEP
from .measure import ARCHIVE


def archive():
    return json.loads(gzip.decompress(ARCHIVE.read_bytes()))


def independently_assemble():
    graph = cell()
    words = [sum(1 << p for p in word) for word in combinations(range(12), 6)]
    index = {mask: i for i, mask in enumerate(words)}
    kinetic, contrast = law.local_term_parts()
    h99 = kinetic + 8*contrast
    np.testing.assert_allclose(h99[9:, :9], 0, atol=2e-13)
    matrix = np.zeros((len(words), len(words)))
    for i, mask in enumerate(words):
        matrix[i, i] = 6*8  # six degree-one private load penalties
        for a, b in graph.seams:
            x, y = (mask >> a) & 1, (mask >> b) & 1
            src = law.SRC.index((x, y))
            matrix[i, i] += h99[src, src]
            if x != y:
                j = index[mask ^ (1 << a) ^ (1 << b)]
                matrix[j, i] += h99[law.SRC.index((y, x)), src]
    return words, matrix


def test_shape_bounds_and_polarity():
    record = archive()
    words, matrix = independently_assemble()
    faces = [sum(1 << v for v in face) for face in cell().faces()]
    both = np.array([any(a & b == 0 and a | b == word for a, b in combinations(faces, 2))
                     for word in words])
    cap = np.diag(matrix) < 134
    assert cap.sum() == 12 and both.sum() == 100 and not any(cap & both)
    controls = record["sharp_standards"]
    assert controls["axis_max"]["cell_q"][0] > 2.399
    assert controls["axis_min"]["cell_q"][0] < -2.399
    assert controls["two_faces"]["cell_q"][1] > 2.399
    assert max(abs(np.array(controls["cap"]["cell_q"]))) < 1e-13
    assert controls["cap"]["centroid_q"][0] < -3.49
    for (a, b), expected in zip([(cap, ~cap), (both, ~(cap | both))],
                                ((10, 1), (18, 6)), strict=True):
        block = matrix[np.ix_(a, b)] < -0.1
        assert (block.sum(axis=1).max(), block.sum(axis=0).max()) == expected
    np.testing.assert_allclose(matrix[np.ix_(cap, both)], 0)
    antipodes = []
    graph = cell()
    for i in range(12):
        near = {i, *graph.neighbours[i]}
        reach = near | {b for a in near for b in graph.neighbours[a]}
        j, = set(range(12))-reach
        antipodes.append(j)
    index = {word: i for i, word in enumerate(words)}
    partner = [index[4095 ^ sum(1 << antipodes[k] for k in range(12) if mask >> k & 1)]
               for mask in words]
    np.testing.assert_array_equal(cap[partner], cap)
    np.testing.assert_allclose(matrix, matrix[np.ix_(partner, partner)], atol=3e-12)
    cap_mask = sum(1 << k for k in record["source_words"]["cap"])
    assert partner[index[cap_mask]] == index[cap_mask]


def test_native_joint_counts_against_binary_delay_tree():
    record = archive()
    words, matrix = independently_assemble()
    values, vectors = linalg.eigh(matrix)
    gate = np.ones(len(words), dtype=complex)
    for stage in range(13):
        gate *= (1+np.exp(-1j*(values-NATIVE_LINE)*STEP*(1 << stage)))/2
    assert 1 << 13 == ARMS
    for source, word in record["source_words"].items():
        idx = words.index(sum(1 << p for p in word))
        for t in (0, 8):
            state = vectors @ (gate*np.exp(-1j*values*t)*vectors[idx])
            probs = abs(state)**2
            bank = record["laws"]["native"]["compact"][source][f"filter_{NATIVE_LINE}_{t}"]
            n = np.array(bank["counts"])
            expected = STARTS*np.append(probs, 1-sum(probs))
            assert np.all(abs(n-expected) <= 8*np.sqrt(np.maximum(expected, 0))+12)
    assert sum(values < 144) == 12 and values[12] > 148


def test_archive_totals_and_changed_laws():
    record = archive()
    for section in record["laws"].values():
        for runs in section["compact"].values():
            for result in runs.values():
                counts = result["counts"]
                assert len(counts) == 925 and sum(counts) == STARTS
                assert counts[-1] == result["dark"] == STARTS-result["hits"]
        for result in section["separated"].values():
            assert len(result["counts"]) == 220**2+1
            assert sum(result["counts"]) == STARTS
            assert result["counts"][-1] == result["dark"]
    zero = record["laws"]["exchange_only"]
    assert max(abs(v) for v in zero["floors"].values()) < 1e-10
    negative = record["laws"]["reversed"]["floors"]
    assert negative["6"] > 2*negative["3"] > 6*negative["1"]
    native = record["laws"]["native"]
    assert native["floors"]["6"] < min(native["fragment_floors"].values())-45
    high = native["high_energy_aperture_probability_bound"]
    bright = native["compact"]["cap"][f"filter_{NATIVE_LINE}_0"]["hits"]/STARTS
    separated = 1-native["separated"][str(NATIVE_LINE)]["dark"]/STARTS
    assert bright > 0.99 and bright > 100*max(high, separated, 0.0001)
    assert all(board["proper_split_lower"] == 176 for board in record["boards"])
    assert all(board["free_connected_probability_lower"] >= 7/8 for board in record["boards"])
    assert record["signed_and_timing_calibration"]["signed_echo"]["1"][3] == STARTS
