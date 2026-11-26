import gzip
import json
from itertools import combinations

import numpy as np

from .. import measure as stock
from . import measure


def test_pair_distance_identity_and_complemented_support():
    xyz = measure.vertices()
    w, _ = stock.sector(2, 2)
    centre = (w != 0) @ xyz / 4
    moment = (w * ((xyz[None, :, :] - centre[:, None, :]) ** 2).sum(axis=2)).sum(axis=1)
    square = {}
    for polarity in (1, -1):
        pair = np.array([np.flatnonzero(row == polarity) for row in w])
        square[polarity] = ((xyz[pair[:, 0]] - xyz[pair[:, 1]]) ** 2).sum(axis=1)
    np.testing.assert_allclose(2 * moment, square[1] - square[-1], atol=2e-14)
    assert len(w) == len(list(combinations(range(12), 4))) * 6


def test_current_meter_matches_a_direct_two_word_commutator():
    # Thirty copies of a seam allow the physical readout's native scaling.
    coords = np.array([-0.31, 0.23])
    rho = np.array([[0.6, 0.27], [0.27, 0.4]])
    swap = np.array([[0.0, 1.0], [1.0, 0.0]])
    h = 15 * (np.eye(2) - swap)
    x = np.diag(coords)
    current = 1j * (h @ x - x @ h)
    bank = measure.inertia(
        np.random.default_rng(731), rho, np.tile([1, 0], (30, 1)), coords, 1
    )
    for u in measure.KICKS:
        kick = np.diag(np.exp(1j * u * coords))
        prepared = kick @ rho @ kick.conj().T
        for key, levels, operator in (
            ("current", bank["current_levels"], current),
            ("energy", [0, 30], h),
        ):
            est, err = measure.read(bank["banks"][str(u)][key], levels)
            assert abs(est - np.trace(prepared @ operator).real) < err + 2e-12


def test_archived_counts_and_calibrated_profiles():
    out = json.loads(gzip.decompress((measure.HERE / "counts.json.gz").read_bytes()))
    for row in out["panels"].values():
        for axis in row["axes"]:
            for raw in axis["banks"].values():
                for bank in raw.values():
                    assert sum(bank["counts"]) + bank["refused"] == measure.STARTS
        for prof in (
            row["radial"],
            *([row["conjugated"]] if "conjugated" in row else []),
        ):
            for banks in prof["banks"].values():
                for bank in banks.values():
                    assert sum(bank["counts"]) + bank["refused"] == measure.STARTS
        for u in (0.25, 0.5):
            estimates = measure.recoil_decode(row["axes"], u)
            # A read is a binary64 dot product of L integer counts below 2**53
            # with levels, divided by their sum. Its error is below
            # (L + 1) * 2**-53 * max|level|: 1e-14 for the energy levels {0, 30}
            # and 7e-15 for the ten current levels below 6. The even work weights
            # energy reads by 4 / u**2 <= 64, so each decoded value errs by less
            # than 1e-12.
            np.testing.assert_allclose(
                estimates, row["recoil"][str(u)], rtol=0, atol=1e-12
            )
            w, dw, j, dj = estimates
            assert abs(w - j) <= dw + dj


def test_exchange_control_against_uniform_stock():
    out = json.loads(gzip.decompress((measure.HERE / "counts.json.gz").read_bytes()))
    v = measure.vertices()
    a, b = stock.EDGES[0]
    l2 = np.sum((v[a] - v[b]) ** 2)
    for name, d in [("neutral", 4), ("positive", 3)]:
        truth = 30 * 2 * d * (12 - d) / (12 * 11) * l2 / (6 * d * d)
        for u in (0.25, 0.5):
            _, _, seen, err = out["panels"][name + "/exchange"]["recoil"][str(u)]
            assert abs(truth - seen) < err


def test_fixed_inventory_counts_and_blind_predictions():
    from . import storage

    raw = json.loads(gzip.decompress((measure.HERE / "counts.json.gz").read_bytes()))
    data = json.loads(gzip.decompress((measure.HERE / "storage.json.gz").read_bytes()))
    for tag, arm in data["settings"].items():
        for name, joint in arm["joint"].items():
            for axis in joint["axes"]:
                for row in axis.values():
                    assert sum(row["counts"]) == row["admissions"][0]
                    assert sum(row["admissions"]) == storage.JOINT_STARTS
            for u in (0.25, 0.5):
                seen, err = storage.decode(joint["axes"], u)
                first = "neutral" if name == "neutral_filled" else "positive"
                part = 0.75 if name.endswith("equal_load") else 1.0
                _, _, origin, os = raw["panels"][f"{first}/{tag}"]["recoil"][str(u)]
                _, _, vac, vs = arm["vacancy_recoil"][str(u)]
                forecast = part**2 * origin + (1 - part) ** 2 * vac
                assert abs(seen - forecast) < err + part**2 * os + (1 - part) ** 2 * vs
