"""Cross-check physical orientation-space records and the saved counts-only bank."""

import hashlib

import numpy as np
from scipy.linalg import expm

from .decode_example import decode
from .experiment import HERE, H, K, dual, exact_check, native_checks
from .packets import occupations, synthesize


def test_local_columns_and_even_sector():
    assert native_checks()["max_amplitude_residual"] < 1e-12
    for n in (5, 6):
        assert exact_check(n)["even_string_trace_error"] < 1e-12


def test_single_wall_record_and_weight():
    n = 7
    words = np.arange(2**n)
    z = 1 - 2 * ((words[:, None] >> np.arange(n)) & 1)
    ends = np.column_stack([np.ones(len(z)), z, -np.ones(len(z))])
    walls = (1 + ends[:, :-1] * ends[:, 1:]) / 2
    mat = np.diag(K * (np.sum(z[:, :-1] * z[:, 1:], axis=1) + z[:, 0] - z[:, -1]))
    for i in range(n):
        mat[words, words ^ (1 << i)] -= H
    x = np.arange(n + 1)
    f = np.exp(-0.1 * (x - 3) ** 2 + 0.8j * x)
    f = synthesize(f / np.linalg.norm(f))
    psi = np.zeros(2**n, dtype=complex)
    for i in x:
        target = np.zeros(n + 1)
        target[i] = 1
        psi[np.flatnonzero(np.all(walls == target, axis=1))[0]] = f[i]
    u, e, v = dual(n, K)
    for t in (0.0, 0.81, 2.37):
        np.testing.assert_allclose(
            occupations(f, u, e, v, t),
            abs(expm(-1j * t * mat) @ psi) ** 2 @ walls,
            atol=1e-12,
        )
    _, physical = np.linalg.eigh(mat)
    actual = np.linalg.norm(physical[:, : n + 1].conj().T @ psi) ** 2
    vacmix = (u - v) / 2
    bound = 1 - ((vacmix**2).sum() - 2 * np.linalg.norm(vacmix.T @ f) ** 2) / 2
    assert actual >= bound - 1e-13 and bound > 0.99


def test_count_bank_checksums_and_independent_decode():
    import json

    receipt = json.loads((HERE / "results.json").read_text())["boards"][0]["observer"]
    bank = HERE / "example_counts.npz"
    with np.load(bank) as rr:
        for key in rr:
            if key in receipt["count_sha256"]:
                assert (
                    hashlib.sha256(
                        np.asarray(rr[key], dtype="<u8").tobytes()
                    ).hexdigest()
                    == receipt["count_sha256"][key]
                )
    records = decode(bank)
    for i, row in enumerate(records["representative_modes"]):
        assert abs(row["omega"] - receipt["modes"][i]["omega"]) < 1e-10
        # Full receipt uses also the higher modes' spatial scans. The example
        # only scans the first three pairs and independently reconstructs k.
        assert abs(row["k"] - receipt["modes"][i]["k"]) < 0.001
