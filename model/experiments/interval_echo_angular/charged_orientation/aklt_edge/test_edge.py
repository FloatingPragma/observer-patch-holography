"""The same detector separates a native charge orbital from a projective edge."""

import json
import numpy as np

from .measure import HERE, SHOTS, run


def check_bank(bank, n):
    counts = np.array(bank["counts"])
    expected = np.array([bank[k] for k in ("bright", "dark", "other")])
    assert counts.sum() == n
    np.testing.assert_allclose(expected.sum(), 1, atol=1e-13)
    # Binomial marginal variances, including nearly impossible outcomes. Exact
    # seeded count equality across numerical libraries has no physical meaning.
    assert np.all(
        abs(counts - n * expected) <= 8 * np.sqrt(n * expected * (1 - expected)) + 2
    )


def test_orders_charges_and_receipts():
    fresh = run()
    saved = json.loads((HERE / "readings.json").read_text())
    for material, stored in zip(fresh["chains"], saved["chains"], strict=True):
        assert material["material"] == stored["material"]
        assert len(set(material["material"]["ports"])) == material["length"]
        for row in (material, stored):
            for key in ("registered_blank", "AKLT"):
                result = row[key]
                finite = result["finite_clock"]
                check_bank(
                    result["endpoint_bank"], SHOTS - result["preparation_refused"]
                )
                check_bank(finite["endpoint_bank"], SHOTS - finite["refused"])
                assert finite["excited_fraction"] < 1e-10
                if key == "AKLT":
                    assert (
                        result["ground_dimension"] == 4
                        and abs(result["edge_commutator"] + 1) < 1e-12
                    )
                    assert (
                        finite["endpoint_bank"]["dark"] > 0.22
                        and finite["endpoint_bank"]["bright"] < 1e-12
                    )
                else:
                    assert result["ground_dimension"] == 3
                    assert finite["endpoint_bank"]["dark"] < 1e-12
                    assert abs(finite["endpoint_bank"]["bright"] - 1 / 3) < 1e-12
            for q in ("-1", "1"):
                line = row["registered_blank"][q]
                assert 8 < line["lowest_energy"] < line["second_energy"] < 8.6
                assert line["same_end_arrival_rank"] == 1
                assert 0 <= line["clock_excited_fraction"] < 1e-4
                n = SHOTS - line["clock_refusals"]
                assert line["clock_arrivals"] + line["clock_other"] == n
                assert abs(
                    line["clock_arrivals"] / n - line["clock_left_probability"]
                ) < 4 / np.sqrt(n)
        for key in ("registered_blank", "AKLT"):
            for prob in ("bright", "dark", "other"):
                np.testing.assert_allclose(
                    material[key]["finite_clock"]["endpoint_bank"][prob],
                    stored[key]["finite_clock"]["endpoint_bank"][prob],
                    atol=1e-12,
                )
