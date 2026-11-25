"""One retained material event plus a physical message must produce a diagonal joint record."""

import json
import numpy as np
from . import shared_report as report


def test_material_report_channels_and_event_counts():
    fresh = report.run()
    saved = json.loads((report.HERE / "shared_report.json").read_text())
    for result in (saved, fresh):
        geo = result["material"]
        assert (
            geo["alice"] != geo["bob"]
            and len(
                set(geo["retained"])
                & {p for route in geo["transfer_paths"] for p in route}
            )
            == 0
        )
        assert len({tuple(w) for w in result["address_words"]}) == 6
        assert all(sum(w) == 5 for w in result["address_words"])
        assert result["conditional_channel"] == [
            dict(
                field_pointer=i,
                bob=i,
                alice=i,
                opposite_face=2 * ((-(i // 2)) % 3) + i % 2,
                no_send=report.REFUSED,
            )
            for i in range(6)
        ]
        assert len(result["rows"]) == 4
        for row in result["rows"]:
            latent = np.array(row["field_event_counts"])
            p = np.array(row["field_event_probabilities"])
            ci = np.array(row["field_event99"])
            assert np.all(ci[:, 0] < p) and np.all(p < ci[:, 1])
            correct, wrong, stopped = (
                row["settings"][k] for k in ("alice", "opposite_face", "no_send")
            )
            np.testing.assert_array_equal(
                np.array(correct["joint"])[:6, :6], np.diag(latent)
            )
            np.testing.assert_allclose(
                np.array(correct["predicted_joint"])[:6, :6], np.diag(p), atol=1e-14
            )
            assert (
                correct["agreements"] == report.SHOTS
                and correct["transmission_refusals"] == correct["disagreements"] == 0
            )
            assert (
                wrong["disagreements"] == sum(latent[2:])
                and wrong["transmission_refusals"] == 0
            )
            wrong_proportion = sum(p[2:])
            assert (
                wrong["disagreement99"][0]
                < wrong_proportion
                < wrong["disagreement99"][1]
            )
            np.testing.assert_array_equal(np.array(stopped["joint"])[:6, -1], latent)
            assert stopped["transmission_refusals"] == report.SHOTS
            for bank in row["settings"].values():
                assert (
                    bank["agreements"]
                    + bank["disagreements"]
                    + bank["transmission_refusals"]
                    == report.SHOTS
                )
    for old, new in zip(saved["rows"], fresh["rows"], strict=True):
        np.testing.assert_allclose(
            old["field_event_probabilities"],
            new["field_event_probabilities"],
            atol=1e-13,
        )
