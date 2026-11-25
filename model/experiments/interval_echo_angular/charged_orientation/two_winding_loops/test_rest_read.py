import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from . import measure, rest_read


def test_known_time_distributions_agree_with_direct_evolution():
    _, ix, terms = measure.build(True)
    h = sum(terms) / 8
    e, v = np.linalg.eigh(h)
    width = 2.0
    weighted = rest_read.wait_kernel(e, width)
    # Independent Gauss quadrature directly on positions and material head words.
    nodes, weights = np.polynomial.legendre.leggauss(32)
    integrated = sum(
        w * measure.choi(measure.kraus(expm(-1j * (16 + width * (1 + x) / 2) * h), ix))
        for x, w in zip(nodes, weights / 2, strict=True)
    )
    by_pole = np.array([measure.kraus(np.outer(u, u), ix) for u in v.T])
    amplitudes = by_pole.transpose(0, 1, 3, 2).reshape(len(e), 2, 4)
    averaged = np.einsum("ij,ila,jlb->ab", weighted, amplitudes, amplitudes.conj()) / 2
    np.testing.assert_allclose(integrated, averaged, atol=1e-12)


def test_local_banks_resolve_channel_signs_and_single_poles():
    out = rest_read.analyze()
    held = json.loads(Path(rest_read.__file__).with_name("rest_read.json").read_text())
    assert out["simple"] and out["positive_energies"] == 8
    assert 1.33320 < out["largest_to_smallest_positive"] < 4 / 3
    assert out["coherent_positive_pt_min"] < -0.007
    assert out["perfect_resolution_pt_min"] > 0.05
    assert min(out["arrival_weight_by_pole"]) > 0
    for j, (row, archived) in enumerate(
        zip(out["finite_waits"], held["finite_waits"], strict=True)
    ):
        assert row["counts"] == archived["counts"]
        assert all(sum(bank) == 2**24 for setting in row["counts"] for bank in setting)
        radius = row["simultaneous_radius"]
        assert row["pt_minimum"] < -radius if j < 2 else row["pt_minimum"] > radius
