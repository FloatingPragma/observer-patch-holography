import json
from functools import lru_cache

import numpy as np
from scipy.linalg import eigh

from . import calibration, identity
from .material import SOURCE, X, Z, formation, generator, mount
from .measure import (
    ACQUISITION,
    HERE,
    N_CLOCK,
    PERIOD,
    SHOTS,
    AcquisitionLaw,
    Reader,
    local,
    run,
)
from .source import inlet


@lru_cache()
def data():
    return json.loads((HERE / "records.json").read_text())


def test_native_charge_calibration_and_spectrum():
    board, d, heads, stock = mount()
    assert stock["total_load"] == 237 and inlet()["words_before"] == [2, 1, 1, 0]
    assert inlet()["words_after"] == [1, 1, 1, 1]
    assert d[1][0] // 12 != d[2][0] // 12
    h = generator(board, d, heads)
    cut = generator(board, d, heads, (1, 2))
    np.testing.assert_array_equal(
        cut, local(SOURCE, 0) - local(X, 1) / 2 - local(X, 2) / 2
    )
    assert np.count_nonzero(np.diag(h - cut)) == 4
    np.testing.assert_allclose(formation(board, d, heads), SOURCE + 104.5 * np.eye(2))
    result = calibration.acquire()
    for arm in result["distal_charge_calibration"]:
        assert arm["charge_bins"]["01"] == arm["charge_bins"]["10"]
        assert arm["charge_bins"]["00"][7] == arm["charge_bins"]["11"][1] == SHOTS
    np.testing.assert_allclose(
        result["source_energy_above_filled"], [104.381966, 106.618034], atol=0.00268
    )


def test_independent_real_time_fourier_instrument():
    b, d, heads, _ = mount()
    h = generator(b, d, heads)
    e, v = eigh(h)
    psi = np.arange(1, 9) + 1j * np.arange(8)
    psi = psi / np.linalg.norm(psi)
    probe = v.T @ local(Z, 1) @ v
    phases = np.exp(-2j * np.pi * np.outer(np.arange(N_CLOCK) / PERIOD, e))
    # Enumerate the physical alternatives before taking the Fourier transform.
    paths = phases[::-1] * ((phases * psi) @ probe.T)
    fft = np.fft.ifft(paths, axis=0)
    low = fft[N_CLOCK // (4 * PERIOD) : 2 * N_CLOCK // PERIOD + 1]
    no, yes = Reader(e, v, 1, 10).records(np.outer(psi, psi.conj()))
    np.testing.assert_allclose(yes, low.T @ low.conj(), atol=2e-12)
    np.testing.assert_allclose(no + yes, fft.T @ fft.conj(), atol=2e-12)


def test_all_count_banks_and_bounds():
    result = run(AcquisitionLaw(ACQUISITION))
    banks = {r["name"]: r for r in result["banks"]}
    for r in result["banks"]:
        old = next(a for a in data()["banks"] if a["name"] == r["name"])
        assert r["flight"] == old["flight"]
        for x, y in zip(r["signed"], old["signed"], strict=True):
            assert x["arm_loss_core"] == y["arm_loss_core"]
            assert sum(x["initial_vacancy_arm_loss_final_vacancy"]) == SHOTS
    active = banks["native_pulse"]
    assert active["flight"][5]["core_x_x"][3] / SHOTS > 0.36
    for order, bank in enumerate(active["signed"]):
        assert (
            bank["arm_loss_core"][6] / SHOTS
            > bank["false_double_loss_lower_core"] + 0.15
        )
        assert (
            bank["initial_vacancy_arm_loss_final_vacancy"][14] / SHOTS
            > bank["false_double_loss_primary_vacancy"] + 0.14
        )
        for cut in ("cut_left_pulse", "cut_right_pulse"):
            assert banks[cut]["signed"][order]["arm_loss_core"][6] == 0
        assert (
            banks["parked_pulse"]["signed"][order][
                "initial_vacancy_arm_loss_final_vacancy"
            ][14]
            == 0
        )


def test_held_species_and_lifetime_controls():
    r = identity.run()
    old = json.loads((HERE / "identity.json").read_text())
    assert r["lifetime"] == old["lifetime"]
    for s in r["charge_only"]:
        i = 1 - s["start"]
        n = s["all_core_left_right"]
        assert n[0] + n[4] == SHOTS
        assert abs(n[i * 4] / SHOTS - s["predicted_conversion"]) < 0.00051
    active = {
        x["wait"]: sum(x["core_left_right"][4:]) for x in r["lifetime"] if x["active"]
    }
    assert (active[32] - active[16]) / SHOTS > 0.31
