"""Measure three ordinary exchange energies and resolve the retained component price."""

from fractions import Fraction as Q
import itertools
import math

import numpy as np
from scipy.linalg import expm

from . import sector_spectrum as spectrum
from . import sector_thermal as thermal


def calibrate(seed=621779):
    words, _, gid, blocks, ops = thermal.apparatus()
    _, _, four_h, h0 = thermal.pair.stock.kernels()
    H = four_h / 4
    levels = np.array([0, 1, 8, 9, 32, 33])
    # These exact native identities imply the integer support for each
    # local energy measurement. The 64-tick Fourier clock resolves it.
    D = four_h - h0
    identity = np.eye(len(D), dtype=np.int64)
    assert np.array_equal(h0 @ h0, 4 * h0) and np.array_equal(h0 @ D, D @ h0)
    assert np.all(D @ (D - 32 * identity) @ (D - 128 * identity) == 0)
    tick = expm((-2j * np.pi / 64) * H)
    unit = np.eye(len(H), dtype=complex)
    amplitudes = []
    for k in range(64):
        amplitudes.append(unit)
        unit = unit @ tick
    P = np.fft.ifft(amplitudes, axis=0)
    distribution = (abs(P) ** 2).sum(1).T
    assert distribution[:, np.setdiff1d(np.arange(64), levels)].max() < 1e-21
    law = distribution[:, levels]
    np.testing.assert_allclose(law.sum(1), 1, atol=1e-12)
    reads = ops[0].diagonal(axis1=1, axis2=2).T
    edges = thermal.pair.stock.carrier()[3]
    codes = np.column_stack(
        [
            9 * words[:, 4 + j] + 3 * words[:, left] + words[:, right]
            for j, (left, right) in enumerate(edges)
        ]
    )
    assert np.array_equal(four_h.diagonal()[codes], reads)
    active = np.unique(codes)
    rng = np.random.default_rng(seed)
    shots = 2**22
    err = 33 * math.sqrt(math.log(2 * len(active) / 0.001) / (2 * shots)) + 1e-10
    assert err < 1 / 8
    measured = np.zeros(len(H), int)
    clocks = {}
    for idx in active:
        counts = rng.multinomial(shots, law[idx] / law[idx].sum())
        mean = counts @ levels / shots
        value = round(4 * mean)
        assert abs(4 * mean - value) + 4 * err < 0.5 and value == four_h[idx, idx]
        clocks[str(idx)] = counts.tolist()
        measured[idx] = value
    assert np.array_equal(measured[codes], reads)
    shapes, _ = spectrum.spectra()
    spectra = [
        (shape, Q(int(reads[g].sum()), 4 * len(g)))
        for shape, g in zip(shapes, blocks, strict=True)
    ]
    keys = sorted(set(spectra))
    spectral_id = np.array([keys.index(spectra[i]) for i in gid])
    rows = []
    for triple in np.unique(reads, axis=0):
        ix = np.all(reads == triple, axis=1)
        labels = np.unique(spectral_id[ix])
        assert len(labels) == 1
        rows.append(
            dict(
                four_times_contact_means=triple.tolist(),
                spectral_id=int(labels[0]),
                word_count=int(ix.sum()),
            )
        )
    pairs = {}
    for cols in itertools.combinations(range(3), 2):
        mapping = {}
        for r in rows:
            key = tuple(r["four_times_contact_means"][i] for i in cols)
            mapping.setdefault(key, set()).add(r["spectral_id"])
        pairs[str(cols)] = sum(len(v) > 1 for v in mapping.values())
    return dict(
        clock_levels=levels.tolist(),
        shots_per_local_input=shots,
        clock_risk=0.001,
        absolute_clock_error=err,
        local_energy_counts=clocks,
        word_to_spectral_class=spectral_id.tolist(),
        contact_triples=rows,
        ambiguous_two_contact_readings=pairs,
        spectral_classes=[
            dict(
                centered_factors=[
                    dict(polynomial=list(map(str, p)), multiplicity_per_word=str(n))
                    for p, n in k[0]
                ],
                mean_energy=str(k[1]),
            )
            for k in keys
        ],
    )
