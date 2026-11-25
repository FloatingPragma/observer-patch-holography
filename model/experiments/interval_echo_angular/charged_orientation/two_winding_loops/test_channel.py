import json

import numpy as np
from scipy.linalg import expm

from . import measure as m
from .certificate import witness


def test_explicit_gate_and_trace_and_native_local_records():
    archive = json.loads((m.HERE / "results.json").read_text())
    for odd, label in [(False, "even"), (True, "odd")]:
        words, ix, terms = m.build(odd)
        step = expm(-1j * sum(terms) / 64)
        propagated = np.linalg.matrix_power(step, 128)
        # Construct the output/reference amplitude with a third explicit
        # index for the discarded source memory. Trace only that third index.
        state = np.zeros((2, 2, 2), complex)
        for s0 in range(2):
            for b in range(2):
                for a in range(2):
                    state[s0, b, a] = propagated[
                        ix[3, 7 * a + 56 * b], ix[0, 7 * s0]
                    ] / np.sqrt(2)
        direct = np.einsum("rbe,sce->rbsc", state, state.conj()).reshape(4, 4)
        np.testing.assert_allclose(m.choi(m.kraus(propagated, ix)), direct, atol=1e-12)
        pt = np.linalg.eigvalsh(m.pt(direct))
        saved = archive[label]
        np.testing.assert_allclose(saved["records"][-1]["min_pt"], pt[0], atol=1e-10)
        bank = saved["tomography"]
        assert np.all(np.asarray(bank["counts"]).sum(axis=2) == m.SHOTS)
        assert abs(bank["measured_min"] - pt[0]) < bank["simultaneous99_radius"]
        assert (bank["measured_min"] + bank["simultaneous99_radius"] < 0) == odd
        assert saved["simple_spectrum"] == odd
        if not odd:
            flip = m.source_reflection(words, ix)
            assert [flip[ix[0, k]] for k in (0, 7)] == [ix[0, k] for k in (7, 0)]
            for a in (0, 7):
                for b in (0, 7):
                    assert flip[ix[3, a + 8 * b]] == ix[3, (7 - a) + 8 * b]


def test_exact_witness_and_phase_shutters():
    saved = json.loads((m.HERE / "results.json").read_text())
    assert witness() == saved["exact_witness"]
    assert json.loads(json.dumps(m.fixture())) == saved["apparatus"]
    for phase in (0, 0.74, -np.pi / 2, np.pi):
        np.testing.assert_allclose(m.phase_lane(phase), np.exp(1j * phase), atol=1e-12)
    assert m.winding.flight.meter() == [0, 1, 2]
