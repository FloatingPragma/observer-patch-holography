"""Check aperture selection, the exact all-open certificate and stored observer banks."""

import hashlib
import json
from itertools import combinations

import numpy as np
from scipy.linalg import eigh
import sympy as sp

from model.candidate_laws.interval_echo import currents
from model.executable_model import cell
from .certificate import PRIME, ROOT, witness
from .measure import HERE, SHOTS, endpoint, fixtures, intervals, turn_preserves


def test_modular_witness_and_apertures():
    x = sp.Symbol("t")
    f = sp.Poly(x**4 + 34 * x**3 - 128 * x**2 - 6904 * x + 41936, x, domain=sp.ZZ)
    assert f.is_irreducible and f.eval(ROOT) % PRIME == 0
    assert f.diff().eval(ROOT) % PRIME != 0
    setups, p, q = fixtures()
    for i, (board, dest) in enumerate(setups):
        assert turn_preserves(board, 0, endpoint(p)["turn"]) == (i == 0)
        assert turn_preserves(board, dest, endpoint(q)["turn"]) == (i == 0)
        for k, distortion in ((8, 98), (0, 85)):
            determinant, unequal_axes = witness(board, dest, p, q, k)
            if i == 0:
                assert not determinant and not unequal_axes
            else:
                assert determinant[0] == (18, 69)
                assert unequal_axes[0] == (18, distortion)


def test_guarded_flight_against_raw_swaps():
    # Independently build the complete 276-word generator by acting on sharp word sets.
    setups, p, q = fixtures()
    board, _ = setups[0]
    words = list(combinations(range(board.n), 2))
    index = {w: i for i, w in enumerate(words)}
    local = list(combinations(range(12), 2))
    src = [index[x] for x in local]
    dst = [index[(a + 12, b + 12)] for a, b in local]
    h = np.zeros((len(words), len(words)))
    for i, (a, b) in enumerate(words):
        for x, y in board.seams:
            if {a, b} == {x, y}:
                h[i, i] -= 2 * currents.KAPPA
            elif (x in (a, b)) != (y in (a, b)):
                w = tuple(sorted({a, b} ^ {x, y}))
                h[i, i] += 0.5
                h[index[w], i] -= 0.5
    energy, vectors = eigh(h)
    receipt = json.loads((HERE / "records.json").read_text())
    home = endpoint(p)
    seam_ids = [local.index(word) for word in cell().seams]
    for entry in receipt["arms"][0]["samples"]:
        transfer = (vectors[dst] * np.exp(-1j * entry["exposure"] * energy)) @ vectors[
            src
        ].T
        assert np.linalg.norm(transfer @ home["frame"][:, 1:]) < 2e-13
        bright = transfer @ home["frame"][:, 0]
        np.testing.assert_allclose(
            np.vdot(bright, bright).real,
            entry["bright_to_whole_receiver_squared"],
            atol=1e-14,
            rtol=2e-9,
        )
        action = transfer @ home["stars"][:, 1]
        bank = next(
            r
            for r in entry["phase_records"]
            if r["source_kind"] == "cross_axis" and not r["receiver_tag"]
        )
        np.testing.assert_allclose(
            abs(action[seam_ids]) ** 2, bank["law"][0][:30], atol=1e-14, rtol=2e-9
        )


def test_acquisition_and_significance():
    receipt = json.loads((HERE / "records.json").read_text())
    assert (
        receipt["protocol_sha256"]
        == hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest()
    )
    for arm in receipt["arms"]:
        for entry in arm["samples"]:
            for record in entry["phase_records"]:
                counts = np.array(record["counts"])
                assert np.all(counts.sum(axis=1) == SHOTS)
                np.testing.assert_allclose(
                    np.sum(record["law"], axis=1), 1, atol=1e-10, rtol=0
                )
                lo0, hi0 = intervals(counts[0, :30])
                lo1, hi1 = intervals(counts[1, :30])
                assert (
                    np.flatnonzero((lo0 > hi1) | (lo1 > hi0)).tolist()
                    == record["resolved_seams"]
                )
                if (
                    arm["independent_turns"][0]
                    or record["source_kind"] == "on_axis"
                    or entry["exposure"] == 0
                ):
                    assert not record["resolved_seams"]
                else:
                    assert record["resolved_seams"]
