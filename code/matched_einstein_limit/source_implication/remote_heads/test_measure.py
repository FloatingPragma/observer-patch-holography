import json
from pathlib import Path

import numpy as np

from . import measure


def test_rank_and_control():
    hs = measure.generators()
    record = json.loads(Path(measure.__file__).with_name("results.json").read_text())
    for modulus, residue in record["determinant_residues"].items():
        rows = measure.derivative_rows(hs, 5, 40, int(modulus))
        minor = ((n, row) for n, row in rows if n in record["orders"])
        orders, det = measure.modular_minor(minor, int(modulus))
        assert orders == record["orders"] and det == residue
    assert measure.exchange_control()["rank"] == 5


def test_signed_count_replay():
    record = json.loads(Path(measure.__file__).with_name("results.json").read_text())
    measured = measure.counted_survey(measure.generators())
    collar = measured["collar_charge"]
    np.testing.assert_allclose(collar["unpair_slope"], 0.2311061985, rtol=0, atol=1e-9)
    assert collar["unread_flight_charge_change"] > 9e-5
    assert abs(collar["unread_flight_parked_energy_change"]) < 2e-13
    low = measured["low_dose"]
    assert low["estimated_min_eigenvalue"] < 0
    high = measured["high_dose"]
    archived = record["finite_read"]["high_dose"]
    assert high["smallest_estimated_eigenvalue"] > high["hilbert_schmidt_radius"]
    for label, sign in (("ordinary", 0.29756316), ("pairing", 0.0)):
        arm = high["held_contrasts"][label]
        assert arm["modular_nats"] > arm["modular_halfwidth"]
        assert abs(arm["parked_energy"] - sign) < arm["energy_halfwidth"]
        np.testing.assert_allclose(
            arm["modular_nats"],
            archived["held_contrasts"][label]["modular_nats"],
            rtol=0,
            atol=2e-5,
        )
