import json

import numpy as np

from . import measure as m


def test_native_full_energy_and_binding_at_fixed_particle_content():
    for a, b in m.MATERIALS:
        h, odd = m.source(a, b)
        np.testing.assert_allclose(odd + np.diag(h.diagonal() - 0.5), h, atol=5e-15)
        np.testing.assert_allclose(odd @ odd, odd, atol=5e-15)
        expected_binding = 16 * a * b + 1 - int(a != b) / 2
        assert h[1, 1] - h[0, 0] == expected_binding
        assert 0 <= min(h.diagonal()) - 0.5 <= max(h.diagonal()) - 0.5 <= m.SCALE - 1
        np.testing.assert_allclose(
            m.beat((a, b)), np.cos(expected_binding * m.SHUTTER / 2) ** 2, atol=1e-14
        )


def test_all_terminal_records_and_silenced_exchange():
    raw = json.loads((m.HERE / "records.json").read_text())
    derived = m.audit(raw)
    assert json.dumps(derived) == json.dumps(
        json.loads((m.HERE / "results.json").read_text())
    )
    for row in derived:
        forecast = row["forecast_from_ramsey"]
        assert all(
            0 < lo < forecast[0] < forecast[1] < hi
            for lo, hi in row["recoil_near_minus_far"]
        )
    for pair in m.MATERIALS:
        rows = {
            r["stock"]: r
            for r in raw["banks"]
            if r["dwell"] == 0 and tuple(r["pair"]) == pair
        }
        full = (
            sum(rows["odd"]["sharp_clicks"])
            + rows["odd"]["parity_clicks"]
            - sum(rows["even"]["sharp_clicks"])
            - rows["even"]["parity_clicks"]
        )
        silent = sum(rows["odd"]["sharp_clicks"]) - sum(rows["even"]["sharp_clicks"])
        assert abs(full / m.N - 1 / m.SCALE) < 2 * m.WIDTH
        assert abs(silent / m.N) < 2 * m.WIDTH
