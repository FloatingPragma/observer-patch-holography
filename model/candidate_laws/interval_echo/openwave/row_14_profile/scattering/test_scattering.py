"""Independent contact, coherent-source, profile, control and complete-bank checks."""

import json
from itertools import product

import numpy as np
import pytest
from scipy.linalg import expm

from . import measure as meter
from model.candidate_laws.interval_echo import candidate

RECORD = json.loads((meter.HERE / "results.json").read_text())


@pytest.mark.parametrize("kappa", [0, 8])
def test_full_seam_address_amplitude_on_coherent_matter(kappa):
    law = candidate.native_law
    hx, hc = law.local_term_parts()
    full_gate = expm(-1j * meter.DURATION * (hx + kappa * hc))
    np.testing.assert_allclose(full_gate[9:, :9], 0, atol=2e-14)
    # An arbitrary pure source, retaining coherence between its nine words.
    words = list(product(range(3), repeat=2))
    rng = np.random.default_rng(7031)
    psi = rng.normal(size=9) + 1j * rng.normal(size=9)
    psi /= np.linalg.norm(psi)
    wave = np.zeros((2, 9, 3, 3), complex)
    for address in (0, 1):
        for j, word in enumerate(words):
            for active, spectator in ((2, 0), (0, 2)):
                initial = law.SRC.index((word[address], active))
                for i, (target_out, probe_out) in enumerate(law.SRC):
                    out = list(word)
                    out[address] = target_out
                    # Equal-amplitude address and probe preparations.
                    wave[address, words.index(tuple(out)), probe_out, spectator] += (
                        full_gate[i, initial] * psi[j] / 2
                    )
    np.testing.assert_allclose(np.linalg.norm(wave), 1, atol=1e-14)
    odd = (wave[:, :, 2, 0] - wave[:, :, 0, 2]) / np.sqrt(2)
    matched = odd.sum(axis=0) / np.sqrt(2)
    charges = (np.array(words) - 1).mean(axis=1)
    np.testing.assert_allclose(matched, meter.gain(kappa) * charges * psi, atol=1e-14)
    mean = np.dot(abs(psi) ** 2, charges)
    power = np.dot(abs(psi) ** 2, charges**2)
    p = [np.linalg.norm((psi + sign * matched) / 2) ** 2 for sign in (1, -1)]
    np.testing.assert_allclose(
        p, meter.detector(mean, power, meter.gain(kappa), True)[:2]
    )
    incoherent = (abs(odd) ** 2).sum() / 2
    expected = np.dot(abs(psi) ** 2, ((np.array(words) - 1) ** 2).sum(axis=1)) / 4
    assert incoherent == pytest.approx(abs(meter.gain(kappa)) ** 2 * expected)


def test_direct_port_profiles_predict_the_independent_phase_detector():
    direct = json.loads((meter.HERE.parent / "records.json").read_text())
    for name, source_name, setting in (
        ("neutral", "neutral", "native"),
        ("positive", "positive", "native"),
        ("neutral_reference_shut", "neutral", "no_reference"),
    ):
        source = direct["panels"][source_name][setting]["profile"]
        panel = RECORD["panels"][name]
        for key, rows in panel["settings"].items():
            center, aperture = key.split("/")
            indices = list(meter.APERTURES[aperture])
            charges = np.array(source[center]["true_signed_means"])
            expected = charges[indices].sum()
            size = sum(np.isin(meter.source.DIST[0], indices))
            for law, k in (("native", 8), ("exchange", 0)):
                row = rows[f"{law}_fringe"]
                p = np.reshape(row["probability"][:-1], (2, 3)).sum(axis=0)
                signed = (
                    size * (p[0] - p[1]) / (panel["acceptance"] * meter.gain(k).real)
                )
                assert signed == pytest.approx(expected, abs=2e-10)
                lo, hi = row["signed_shell_interval"]
                assert lo <= expected <= hi


def test_neutral_forward_dark_recording_paths_and_compact_scattering():
    neutral = RECORD["panels"]["neutral"]["settings"]
    all_ports = neutral["intrinsic/all"]
    for key in ("native_elastic", "exchange_elastic"):
        assert sum(all_ports[key]["probability"][:-1:2]) == 0
        assert sum(all_ports[key]["counts"][:-1:2]) == 0
    acceptance = RECORD["panels"]["neutral"]["acceptance"]
    recorded = all_ports["recorded_address"]["probability"]
    assert sum(recorded[:-1:2]) / acceptance == pytest.approx(
        abs(meter.gain(8)) ** 2 / 36
    )
    positive = RECORD["panels"]["positive"]
    matched = positive["settings"]["intrinsic/all"]["native_elastic"]["probability"]
    assert sum(matched[:-1:2]) / positive["acceptance"] == pytest.approx(
        abs(meter.gain(8)) ** 2 / 144
    )
    center = neutral["intrinsic/center"]
    for key, bounds in (
        ("native_elastic", (0.51, 0.52)),
        ("exchange_elastic", (0.0003, 0.0004)),
    ):
        p = center[key]["probability"][2:4]
        assert bounds[0] < p[0] / sum(p) < bounds[1]
    p = center["native_elastic"]["probability"]
    assert 0.0044 < sum(p[:-1:2]) / acceptance < 0.0046
    for aperture in ("center", "outer"):
        row = neutral[f"intrinsic/{aperture}"]["native_fringe"]
        size = 1 if aperture == "center" else 6
        lo, hi = meter.decode(row, RECORD["standards"]["native"]["2"], size, 1)
        assert lo < 0 < hi


def test_all_counts_and_custody():
    for level, board in enumerate(("20", "80", "320")):
        assert meter.routes(level) == RECORD["routes"][board]
    banks = [b for s in RECORD["standards"].values() for b in s.values()]
    for p in RECORD["panels"].values():
        for s in p["settings"].values():
            banks.extend(s.values())
    for bank in banks:
        assert sum(bank["counts"]) == RECORD["starts"]
        assert sum(bank["grid_weights"]) == RECORD["grid"]
        law = np.array(bank["probability"])
        assert min(law) >= 0 and abs(sum(law) - 1) < 2e-11
        assert (
            np.abs(np.array(bank["grid_weights"]) / RECORD["grid"] - law).sum() < 3e-11
        )
    with pytest.raises(ValueError):
        meter.verify({"counts": [0, 2]}, {"counts": [1, 1]})


def test_two_probe_population_intervals_on_compact_events():
    from .readout import center_composition

    banks = RECORD["panels"]["neutral"]["settings"]["intrinsic/center"]
    standard = RECORD["standards"]["native"]["2"]
    inferred = center_composition(
        banks["native_fringe"], banks["native_elastic"], standard, 1
    )
    # Independent local pointer populations on the earlier connected bank.
    archive = json.loads((meter.HERE.parent / "records.json").read_text())
    profile = archive["panels"]["neutral"]["native"]["profile"]["intrinsic"]
    labels = np.array(profile["keys"])
    connected = labels[:, 4] == 1
    counts = np.array(profile["counts"])
    for label, charge in (("vacancy", -1), ("excess", 1)):
        prob = (
            counts[connected & (labels[:, 0] == charge)].sum() / counts[connected].sum()
        )
        # The independent original cohort has only about 23,500 connected
        # records. Its own sampling radius applies when comparing the two banks.
        err = np.sqrt(np.log(400) / (2 * counts[connected].sum()))
        lo, hi = inferred[label]
        assert 0.25 < lo < hi < 0.27
        assert lo < prob + err and prob - err < hi
