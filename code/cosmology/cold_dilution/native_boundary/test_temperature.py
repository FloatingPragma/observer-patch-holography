"""Direct comparison of the archived peak arithmetic with observed bandpowers."""

import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).parent


def test_measured_bandpowers_and_control():
    record = json.loads((HERE / "results.json").read_text())
    raw = HERE / "planck_tt.dat"
    assert hashlib.sha256(raw.read_bytes()).hexdigest() == record["planck_tt"]["sha256"]
    ell, powers = np.loadtxt(raw, usecols=(0, 1), unpack=True)
    ix = [np.argmin(abs(ell - x)) for x in (225, 525, 825)]
    observed = record["measured"]
    np.testing.assert_allclose(observed["measured_peaks"], powers[ix])
    np.testing.assert_allclose(
        observed["third_over_second"], powers[ix[2]] / powers[ix[1]]
    )
    control = record["variants"]["ordinary_cold_control"]
    layer = record["variants"]["maximal_one_layer"]
    np.testing.assert_allclose(
        layer["omega_c"], record["control_background"]["ombh2"] * 17 / 260
    )
    for panel in (control, layer):
        np.testing.assert_allclose(
            panel["at_observed_crests_aligned"][0], powers[ix[0]]
        )
    assert abs(control["third_over_second"] - observed["third_over_second"]) < 0.002
    assert layer["at_observed_crests_aligned"][2] < powers[ix[2]] / 2
