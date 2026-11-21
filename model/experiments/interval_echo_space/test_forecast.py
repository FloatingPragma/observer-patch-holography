"""Independent checks on the anonymous survey and its held receiver forecasts."""
from itertools import combinations

import networkx as nx
import numpy as np

from .forecast import HOP, infer_operator
from .instrument import HERE, WORDS, slow_frame


def test_survey_geometry_and_slow_projector_are_equivariant_under_anonymous_labels():
    with np.load(HERE / "forecast_records.npz", allow_pickle=False) as records:
        survey = records["n80_survey"]
    rng = np.random.default_rng(478)
    permutation = np.concatenate([12 * i + rng.permutation(12) for i in range(10)])
    reordered = survey[np.ix_(np.r_[permutation, len(permutation)], permutation)]
    original, _, _ = infer_operator(survey)
    permuted, _, _ = infer_operator(reordered)
    np.testing.assert_array_equal(permuted, original[np.ix_(permutation, permutation)])
    frame, _ = slow_frame(survey[:12, :12])
    reordered_frame, _ = slow_frame(reordered[:12, :12])
    projector = frame @ frame.T
    np.testing.assert_allclose(reordered_frame @ reordered_frame.T, projector[np.ix_(permutation[:12], permutation[:12])],
                               atol=1e-12, rtol=0)


def test_three_layer_rooted_surveys_distinguish_boards_at_the_outer_aperture():
    with np.load(HERE / "forecast_records.npz", allow_pickle=False) as records:
        graphs = []
        for size in (20, 80, 320):
            _, degrees, adjacency = infer_operator(records[f"n{size}_survey"])
            g = nx.from_numpy_array(adjacency)
            for i in g:
                g.nodes[i]["type"] = (int(degrees[i]), int(i >= 12) + int(i >= 48))
            graphs.append(g)
    for ports, expected in ((12, [True, True, True]), (48, [False, False, True]), (120, [False, False, False])):
        same = [nx.vf2pp_is_isomorphic(a.subgraph(range(ports)), b.subgraph(range(ports)), node_label="type")
                for a, b in combinations(graphs, 2)]
        assert same == expected


def test_held_forecast_reconstructed_without_access_to_any_held_receiver_table():
    with np.load(HERE / "forecast_records.npz", allow_pickle=False) as records:
        survey = records["n20_survey"]
        h, _, _ = infer_operator(survey)
        # Independent Hermitian spectral propagation from the survey alone.
        spectrum, modes = np.linalg.eigh(h)
        transfers = (modes[12:48] * np.exp(-3j * spectrum)) @ modes[:12].T
        clicks = abs(transfers.reshape(3, 12, 12))**2
        ports = np.concatenate((clicks, 1 - clicks.sum(axis=1, keepdims=True)), axis=1)
        weights = 1 / 12 + WORDS @ records["n20_frame"].T / 24
        predicted = np.einsum("ws,ros->wro", weights, ports)
        np.testing.assert_allclose(predicted, records["n20_t3_forecast"], atol=1e-12, rtol=0)
        assert np.max(abs(predicted - records["n20_t3_held"])) < 6.32e-7
        # Doubling the calibrated rate resolves the same support but changes its forecast.
        wrong = (modes[12:48] * np.exp(-3j * spectrum / HOP)) @ modes[:12].T
        wrong = np.einsum("ws,ros->wro", weights, abs(wrong.reshape(3, 12, 12))**2)
        assert np.max(abs(wrong - records["n20_t3_held"][:, :, :12])) > 0.001
