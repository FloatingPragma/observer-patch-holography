import json

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.experiments.interval_echo_einstein_station.collar_pairing.ranging import (
    measure as radar,
)
from . import chart, returns


def test_two_charts_and_held_fifth_position():
    output = chart.analyze()
    for panel in output.values():
        assert all(lo > 0 and hi < 1 for lo, hi in panel["weights"])
        assert (
            sum(lo for lo, _ in panel["weights"])
            < 1
            < sum(hi for _, hi in panel["weights"])
        )
        assert panel["affine_residual_length"][0] > 0.0014
    assert output["pooled"]["directional_moment_eigenvalues_between"] == [0.29, 0.40]
    lo, hi = output["pooled"]["fifth_position_prediction_error"]
    assert 0.0073433 < lo < hi < 0.0073435


def test_detector_counts_and_independent_ordinary_standard():
    archive = json.loads((chart.HERE / "results.json").read_text())
    records = returns.report(
        archive["geometry"], archive["probabilities"], archive["raw"]
    )
    for q in ("vacancy", "excess"):
        for quantity, interval in records[q]["counts"].items():
            a, b = records[q]["native"][quantity]
            assert interval[0] < a <= b < interval[1]
            np.testing.assert_allclose(
                interval, archive["readout"][q]["counts"][quantity], rtol=0, atol=1e-13
            )
    assert records["vacancy"]["counts"]["equivalent_exchange_energy"][0] > 2.83
    assert records["excess"]["counts"]["equivalent_exchange_energy"][1] < 0.164


def test_independent_native_exponential_at_the_five_local_emitters():
    archive = json.loads((chart.HERE / "results.json").read_text())
    ports = archive["raw"]["ports"]
    board = radar.formation.federation(1)
    marks = radar.formation.region(board)
    rho, tau, _ = radar.prepared_states(board, marks, 1)
    ix = [4 * p + j for p in ports for j in range(4)]
    for q, name in ((-1, "vacancy"), (1, "excess")):
        rows = {}
        for odd in (False, True):
            h = radar.generator(board, marks, q, 8, odd)
            src = np.zeros((h.shape[0], len(ix)))
            src[ix] = np.eye(len(ix))
            wave = expm_multiply(-4j * h, src, traceA=-4j * h.diagonal().sum())[
                ix
            ].reshape(5, 4, 5, 4)
            local = np.stack([wave[i, :, i, :] for i in range(5)])
            for s, coin in (("formed", rho), ("independent", tau)):
                pop = np.einsum("yoi,ij,yoj->y", local, coin, local.conj()).real
                if odd and s == "independent":
                    rows["ordinary"] = (rows[s] + pop) / 2
                elif not odd:
                    rows[s] = pop
        for s, data in rows.items():
            exact = np.mean(archive["probabilities"][name][s], axis=1)
            # The sparse exponential is an independent binary64 path. The
            # tight native enclosures are provided by the integer propagator.
            np.testing.assert_allclose(data, exact, rtol=0, atol=1e-10)
