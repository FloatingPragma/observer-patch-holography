"""Cross-check the physical weak strand, source populations and finite end clock."""

from hashlib import sha256
import json
import numpy as np
from scipy.linalg import eigh, expm
from . import weak_strand as ws, strand_clock as clock
from .measure import HERE
from ..joined_decay.measure import core_pulse
from ..passive_decay.measure import hamiltonians


def record(module):
    return json.loads((HERE / f"{module}.json").read_text())


def test_periodic_obstruction_and_wordwise_ordinary_energy():
    data = record("weak_strand")
    for name, digest in data["protocol_sha256"].items():
        assert digest == sha256((HERE / name).read_bytes()).hexdigest()
    for n in (6, 8):
        b, d, heads, _ = ws.ring_layout(n, True, False)
        a, z = d[3]
        used = set(sum(d[:5], ())) | {
            p for edge in ws.route((3,), crossing=True)[2] for p in edge
        }
        assert set(b.neighbours[a]) & set(b.neighbours[z]) <= used
        assert len(set(b.neighbours[z]) - used) == 1
        assert all(p // 12 != 6 for p in d[-1])
        h, _ = ws.constructed(d, b, heads, "native")
        # Direct one-word toggles independent of sparse flip construction.
        words = np.arange(2 ** len(d))
        for p in (0, 43, len(words) - 1):
            v = h[:, p].toarray().ravel()
            assert np.all(v[p ^ (1 << np.arange(len(d)))] == -0.5)
            assert np.count_nonzero(v[:p]) + np.count_nonzero(v[p + 1 :]) == len(d)


def test_source_at_the_primary_clock_by_independent_dense_evolution():
    b, d, heads, _ = ws.ring_layout(6, True, False)
    h, _ = ws.constructed(d, b, heads, "native")
    energies, v = eigh(h.toarray())
    _, f = eigh(hamiltonians(4)[0])
    wave = (core_pulse(f) @ v[:, 0].reshape(4, -1)).ravel()
    state = v @ (np.exp(-32j * energies) * (v.T @ wave))
    medium = ws.strand(6).toarray()
    e, w = eigh(medium)
    lines = (e - e[0] > 1e-10) & (e - e[0] < 0.934101)
    projector = w[:, lines] @ w[:, lines].T
    expected = record("weak_strand")["settings"]["6"]["source"]["trials"]["native"][
        "core_pulse"
    ][-1]
    for i in range(4):
        remnant = f[:, i] @ state.reshape(4, -1)
        p = np.linalg.norm(remnant.reshape(2, -1) @ projector) ** 2
        np.testing.assert_allclose(expected["core_and_band"][i], p, atol=3e-12)
    source = record("strand_clock")["settings"]["6"]
    assert (
        source["pulse_and_idle_cut_work"]["native"][0]
        < source["pulse_and_idle_cut_work"]["native"][1]
    )


def test_fourier_energy_instrument_against_time_controlled_local_operation(monkeypatch):
    monkeypatch.setattr(clock, "LABELS", 32)
    monkeypatch.setattr(clock, "PERIOD", 8)
    n = 3
    h = ws.strand(n).toarray()
    effects, v, gaps, band, false = clock.instrument(h, 1, n, 0.934101)
    rng = np.random.default_rng(43214)
    psi = rng.normal(size=len(v)) + 1j * rng.normal(size=len(v))
    psi /= np.linalg.norm(psi)
    z = np.diag(1 - 2 * ws.fields(n)[0][:, 1])
    duration = 2 * np.pi * (clock.LABELS - 1) / clock.PERIOD
    amplitudes = []
    for label in range(clock.LABELS):
        t = 2 * np.pi * label / clock.PERIOD
        evolved = expm(-1j * (duration - t) * h) @ z @ expm(-1j * t * h) @ psi
        amplitudes.append(v[:, 0] @ evolved)
    outputs = abs(np.fft.ifft(amplitudes)) ** 2
    frequency = np.fft.fftfreq(clock.LABELS) * clock.PERIOD
    accepted = (frequency > 0) & (frequency < 0.934101)
    av = v.T @ psi
    np.testing.assert_allclose(
        outputs[accepted].sum(), np.vdot(av, effects[2] @ av).real, atol=5e-13
    )
    np.testing.assert_allclose(
        1 - outputs.sum(), np.vdot(av, effects[0] @ av).real, atol=5e-13
    )
    assert band.any()
    np.testing.assert_allclose(
        false, eigh(effects[2][np.ix_(~band, ~band)], eigvals_only=True)[-1], atol=1e-13
    )


def test_clock_core_coincidences_and_material_charge_standards():
    data = record("strand_clock")
    assert (
        data["protocol_sha256"]
        == sha256((HERE / clock.DESIGN).read_bytes()).hexdigest()
    )
    linear = record("strand_lines")
    assert (
        linear["protocol_sha256"]
        == sha256((HERE / "STRAND_LINES.md").read_bytes()).hexdigest()
    )
    for n in (6, 8):
        tested = data["settings"][str(n)]
        end = tested["results"][str(n - 1)]
        assert all(s["binary"] == [3, 0, 0, 1] for s in tested["inventory"])
        true = eigh(ws.strand(n).toarray(), eigvals_only=True)
        for r in linear["compared"][str(n)]:
            assert abs(true[r["level"]] - true[0] - r["gap"]) < r["radius"] < 0.00005
        for setting, banks in end["rows"].items():
            for rec in banks.values():
                assert abs(np.array(rec["core_by_bin"]).sum() - 1) < 1e-10

        def lower(setting, arm):
            return np.array(end["rows"][setting][arm]["core_by_bin"])[1:3, 2].sum()

        assert lower("native", "pulse") > 0.011
        assert lower("source_cut", "pulse") < 0.0005
        assert lower("native", "idle") < 1e-9
        assert end["false_upper"] < 0.0001
        reset = end["rows"]["native"]["pulse"]["reset_core_soft"]
        assert sum(reset[1:3]) < 0.000009
