"""Independent charge, propagation and fragmentation checks for the two full-sector audits."""

import json
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.executable_model import carrier
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls import (
    belt,
    ternary,
    upper,
)

HERE = Path(__file__).parent


def test_ternary_receipt_and_histogram_independent_colour_certificate():
    recorded = json.loads((HERE / "ternary.json").read_text())
    assert ternary.run() == recorded
    for case in recorded["belts"]:
        channels = case["reservoir_delta_holes_accepted"]
        for counts in channels.values():
            assert sum(counts.values()) == pytest.approx(1, abs=1e-9)
        if case["total"] == 13 and case["holes"] == 1:
            assert 0.099 < channels["1"]["1,0,1"] < 0.100
            assert 0.099 < channels["2"]["1,1,1"] < 0.100
    assert max(recorded["local_histogram"]["reservoir_swap_error"].values()) < 1e-10
    for total, size in ((13, 4740), (14, 2850)):
        blocks = [b for b in recorded["all_histograms"] if b["total"] == total]
        assert (
            sum(b["dimension"] * (1 if b["momentum"] in (0, 5) else 2) for b in blocks)
            == size
        )
        assert all(b["krylov_rank"] == b["dimension"] for b in blocks)
    cases = recorded["long_cores"]["cases"]
    for s in cases:
        assert (s["right_clock"] - s["left_clock"]) % 3 == 1
        assert s["holes"] == s["counted_vacancies"]
        charge = s["core_window_excess"] / s["unit"]
        assert charge == pytest.approx(-1 / 3 if s["shift"] == 1 else 2 / 3)
        strength = (
            s["phase_records"]["8.0"]["phase"]
            / recorded["long_cores"]["independent_unit"]["phase"]
        )
        assert abs(strength - charge) < 0.005


def test_ternary_propagation_against_direct_seam_action():
    case = json.loads((HERE / "ternary.json").read_text())["belts"][1]
    states = list(ternary.words(3, 1))
    addr = {s: i for i, s in enumerate(states)}
    edges = belt.ring(10).seams
    rows, cols, data = [], [], []
    local = law.seam_term(8)
    for j, s in enumerate(states):
        diagonal = 16 * sum(v * v for v in s)
        for a, b in edges:
            at = law.SRC.index((s[a], s[b]))
            diagonal += local[at, at]
            if s[a] != s[b]:
                t = list(s)
                t[a], t[b] = t[b], t[a]
                rows.append(addr[tuple(t)])
                cols.append(j)
                data.append(local[law.SRC.index((t[a], t[b])), at])
        rows.append(j)
        cols.append(j)
        data.append(diagonal)
    h = sparse.csr_matrix((data, (rows, cols)), shape=(len(states), len(states)))
    state = np.zeros(len(states))
    state[addr[tuple(case["sharp_source"])]] = 1
    final = expm_multiply(-0.25j * h, state)
    mask = [bool(ternary.addresses(s, 3)) for s in states]
    assert sum(abs(final[mask]) ** 2) == pytest.approx(
        case["joint_probabilities"]["8.0"][1], abs=2e-10
    )


def test_binary_compression_against_full_carrier_sector():
    words, cost, degree, edges = upper.product_space(upper.bins(), 2, [(4, 11)], 40, 2)
    direct = BinarySector(carrier.Carrier("linked belt test", 20, edges), 2)
    order = [direct.index[tuple(i for i in range(20) if w & (1 << i))] for w in words]
    sub = sparse.diags(degree / 2 - 16 * cost) - 0.5 * upper.hopping(words, edges)
    np.testing.assert_allclose(
        direct.hamiltonian(8)[np.ix_(order, order)], sub.toarray()
    )
    # A conservative estimate displaced below the top must fail the Gram check.
    with pytest.raises(np.linalg.LinAlgError):
        upper.checked_ceiling(
            sub.toarray(), np.linalg.eigvalsh(sub.toarray())[-1] - 0.01
        )


def test_fragment_witnesses_and_full_replay():
    saved = json.loads((HERE / "upper.json").read_text())
    replay = upper.run()
    delta = replay["ordering"]
    assert delta["absolute_positive_minus_neutral"] == pytest.approx(
        saved["ordering"]["absolute_positive_minus_neutral"]
    )
    assert (
        32.86
        < delta["absolute_positive_minus_neutral"][0]
        < delta["absolute_positive_minus_neutral"][1]
        < 33.25
    )
    assert delta["positive_minus_neutral_with_rest_unit"][1] < -14.75
    for ref, got in zip(saved["sectors"], replay["sectors"]):
        for key in (
            "trial_energy",
            "top_eigenstate_wall_probability_lower",
            "complement_ceiling",
            "pair_plus_single_upper",
        ):
            assert got[key] == pytest.approx(ref[key], abs=2e-8)
        np.testing.assert_equal(got["contact_read"], ref["contact_read"])
        edge = got["joined_ceiling"]["upper"]
        assert got["pair_plus_single_upper"] + 0.22 < got["trial_energy"] < edge
        assert edge + 0.23 < got["capped"]["pair_plus_single_numerical_edge"]
        assert got["high_energy_filter"]["all_times_wall_probability_lower"] > 0.57
        assert got["high_energy_filter"]["coherent_trial_acceptance_lower"] > 0.36
        assert (
            got["contact_read"]["energies"]["0.0"]["mean"]
            < got["exchange_only_three_separated_edge"]
        )
        assert got["contact_read"]["energies"]["8.0"]["mean"] == pytest.approx(
            got["trial_energy"], abs=1e-8
        )
