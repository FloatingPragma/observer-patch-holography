import json

import numpy as np

from ..receipt_replay import check_binomial, pinned_receipt
from . import certificate, counts, measure


def test_five_pixel_first_arrivals_and_correlation():
    source = measure.front.optics.pairing.formation
    board = source.federation(1)
    route = source.region(board)
    rings = measure.screens(board, route)
    states = measure.aperture.survey.preparations()
    for q in (-1, 1):
        kernels = measure.effects(board, route, rings, q, None)
        measure.check_leading(kernels, states["pairing"])
        weights = measure.probabilities(kernels, states["equilibrium"])
        cov = weights[0, 0] - weights[0, 1] * weights[1, 0]
        j, _ = measure.image(weights)
        assert (
            cov > 0
            and np.linalg.norm(j - np.diag([50 * cov / weights.sum() ** 2, 0])) < 2e-15
        )
    census = measure.gateway_census()
    assert [c["bridges"] for c in census] == [30, 120, 480]


def test_exact_source_and_finite_dwell_replay():
    exact = certificate.run()
    assert exact == json.loads((measure.HERE / "certificate.json").read_text())
    record = json.loads((measure.HERE / "results.json").read_text())
    for q, data in exact.items():
        lo, hi = json.loads(data["ordinary_calibrated_ratio"])
        nominal = record["charges"][q]["finite"]["4.0"]
        assert lo < nominal["observed_to_predicted"]["unpairing"] < hi
        assert 0.53 < lo < hi < 0.58 if q == "-1" else -0.50 < lo < hi < -0.45
        for name in ("ordinary", "unpairing"):
            low, high = json.loads(data[name])
            assert low < nominal["area_slopes"][name] < high
        assert min(json.loads(data["signed_area"])) > 0
        assert (
            abs(
                record["charges"][q]["zero_contrast_at_four"]["area_slopes"][
                    "unpairing"
                ]
            )
            < 1e-15
        )


def test_image_area_from_held_finite_counts():
    laws = counts.detector_laws()
    fresh = counts.run(laws=laws)
    frozen = pinned_receipt(
        measure.HERE / "counts.json",
        "97d61082df75b5baf52bf31de048b0545fa2ebf3efd6a40d5824a77875058ff3",
    )
    assert {k: v for k, v in fresh.items() if k != "species"} == {
        k: v for k, v in frozen.items() if k != "species"
    }
    # Each of the six outcomes is a binomial marginal. Test every bin in
    # both independent replays, then reconstruct the signed intervals.
    log_level = np.log(2 * 2 * 2 * 5 * 5 * 6 / 0.01)
    for raw in (frozen, fresh):
        shots = raw["shots_per_pixel_per_preparation"]
        for q, panel in raw["species"].items():
            observed = panel["hit_pixels_then_miss"]
            for setting, rows in observed.items():
                for row, law in zip(rows, laws[q][setting].T, strict=True):
                    assert sum(row) == shots
                    for value, p in zip(row, [*law, 1 - law.sum()], strict=True):
                        check_binomial(value, shots, p, log_level)
            rebuilt = counts.area_intervals(observed, shots)
            for label, interval in rebuilt.items():
                # Centering translations preserve the determinant exactly.
                # The interval widths can depend on summation order in BLAS.
                np.testing.assert_allclose(
                    json.loads(interval),
                    json.loads(panel["signed_slopes_99pct"][label]),
                    rtol=1e-10,
                    atol=1e-13,
                )
            ordinary = json.loads(rebuilt["ordinary"])
            collar = json.loads(rebuilt["pairing"])
            assert min(collar) > 0
            assert min(ordinary) > 0 if q == "-1" else max(ordinary) < 0
