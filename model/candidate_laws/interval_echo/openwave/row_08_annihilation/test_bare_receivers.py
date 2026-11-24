"""Independent source conservation, detector and propagation checks."""

import gzip
import json

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from . import bare_receivers as bare, experiment as ex, walls


def test_local_parity_readout_and_native_count_conservation():
    kinetic, contrast = ex.LAW.local_term_parts()
    labels, table = bare.analyzers()
    for k in (0, 8):
        full = kinetic + k * contrast
        h = full[:9, :9]
        np.testing.assert_allclose(full[9:, :9], 0, atol=2e-13)
        np.testing.assert_allclose(h, ex.LAW.seam_term(k), atol=2e-13)
        parity = np.asarray([v for a, b, v in labels])
        values = [k * (a - b) ** 2 + (v == -1) for a, b, v in labels]
        np.testing.assert_allclose(table.T @ h @ table, np.diag(values), atol=2e-13)
        np.testing.assert_allclose(
            expm(-1j * np.pi * h), (table * parity) @ table.T, atol=3e-13
        )
        for load in (0, 1, 2):
            histogram = np.array(
                [int(i == load) + int(j == load) for i, j in ex.LAW.SRC]
            )
            np.testing.assert_allclose(
                (histogram[:, None] - histogram[None, :]) * h, 0, atol=2e-13
            )


def test_extinction_with_remote_brightness_exists_in_larger_charge_sector():
    full_codes, _ = bare.basis(True)
    narrow_codes, _ = bare.basis(False)
    ring = walls.loads(walls.reference(0))
    allowed_by_total = (ring + [0, 2, 1, 1]) @ bare.POWERS
    assert allowed_by_total in full_codes
    assert allowed_by_total not in narrow_codes
    assert (ring + [1] * 4) @ bare.POWERS in narrow_codes
    assert len(narrow_codes) == 34650
    assert len(full_codes) == 73789


def test_wavefunction_envelope_against_independent_sparse_exponential():
    codes, words = bare.basis(False)
    h = bare.generator(codes, words, 8, False)
    start = np.zeros((len(codes), 1))
    start[
        np.searchsorted(codes, (walls.loads((0, 1, 1, 0)) + [1] * 4) @ bare.POWERS), 0
    ] = 1
    actual = expm_multiply(-0.5j * h, start)
    estimated, _ = bare.flights(h, start, (0.5,))
    # The Chebyshev propagation drops its common scalar phase.
    np.testing.assert_allclose(
        abs(np.vdot(actual, estimated[0])), 1, atol=2e-12, rtol=0
    )


def test_saved_records_keep_each_failure_and_active_comparison():
    data = json.loads(
        gzip.decompress((ex.HERE / "bare_receivers.json.gz").read_bytes())
    )
    for case in data["cases"]:
        for setting in case["readings"]:
            for row in setting["preparations"]:
                count = np.asarray(row["count"])
                assert count.shape == (5, 9, 9)
                assert count.sum() == data["shots"]
                if case["pair_mixing"]:
                    assert row["ordered_and_any_energy"] > 0
                else:
                    assert row["ordered_and_any_energy"] == 0
                for count in row["outside_material_counts"]:
                    assert 0 <= count <= data["shots"]
    assert (
        data["cases"][0]["readings"][-1]["preparations"][1]["outside_material_counts"][
            0
        ]
        > 0
    )
