import gzip
import json

import numpy as np

from .. import experiment
from . import growth, single_write, unsorted


def test_repair_branches_and_coverage():
    observed = growth.instrument()
    data = json.loads(gzip.decompress((unsorted.HERE / "records.json.gz").read_bytes()))
    assert observed == data["law"]
    assert sum(data["first_pointer_counts"]) == 2**20
    waves = 0
    for row in data["boards"]:
        setup = growth.commission(row["level"], row["first_pointer"])
        assert setup["complete"] and setup["reads"] == row["setup"]["reads"]
        assert setup["vacancy_route_sha256"] == row["setup"]["vacancy_route_sha256"]
        assert all(
            obs["pointer"][-1] == row["first_pointer"]
            and len(obs["reset"]) == len(obs["pointer"]) - 1
            and all(b != row["first_pointer"] for b in obs["pointer"][:-1])
            for obs in setup["reads"]
        )
        assert len(setup["reads"]) == setup["connectors"]
        assert setup["completion_probability"] > 0.9999
        w = row["wave"]
        assert w["spatial"]["status"] == "pass"
        assert w["max_lattice_error"] < 0.004 and w["max_swapped_error"] > 0.35
        waves += len(w["spatial"]["even_trained_k"])
    assert waves == 34


def test_unfiltered_detector_with_every_conditional_word():
    r = json.loads((unsorted.HERE / "single_write.json").read_bytes())
    with np.load(unsorted.HERE / "single_write.npz", allow_pickle=False) as rolls:
        counts = rolls["counts"]
        assert counts.shape == (32, 5, 91, 4)
        assert np.all(counts.sum(axis=(0, 3)) == single_write.SHOTS)
        read = single_write.summarize(counts)
        np.testing.assert_allclose(read["pooled_norm"], r["pooled_norm"], atol=1e-9)
    eps = read["norm_error_bound"]
    np.testing.assert_allclose(np.asarray(read["pattern_norm"])[[0, -1]], 1, atol=eps)
    assert all(ns[-1] < 1 - eps for ns in read["pattern_norm"][1:-1])
    assert read["pooled_norm"][-1] + eps < read["mean_conditional_norm"][-1] - eps


def test_uniform_material_from_full_orientation_space_and_coherent_strand_agree():
    gaps, weights, scales, names = unsorted.spectrum()
    p = np.sqrt(2 / (unsorted.M + 1)) * np.sin(
        np.pi * np.outer(np.arange(1, unsorted.M + 1), [1, 2]) / (unsorted.M + 1)
    )
    for bit, idx in [(0, 0), (1, -1)]:
        u, e, _ = experiment.dual(unsorted.N, 8 + bit / 4)
        for t in (0.4, 2.3, 9.1):
            aa = (u * np.cos(e * t)) @ u.T
            v = p.T @ aa @ p
            a, b = (aa @ p).T
            expected = np.r_[
                np.linalg.det(v), v[0, 0] * b - v[0, 1] * a, v[0, 1] * b - v[1, 1] * a
            ]
            physical = np.cos(t * gaps[idx]) @ weights[idx, :, : len(expected)]
            np.testing.assert_allclose(physical, expected, atol=1e-11)
    assert len(scales) == 75 and len(names) == 5


def test_whole_profile_gate_keeps_failed_single_write_result():
    rows = json.loads((unsorted.HERE / "unsorted.json").read_bytes())["panels"]
    unsorted_row, *uniform = rows
    assert unsorted_row["spatial"]["status"] == "physical_fail"
    assert (
        unsorted_row["spatial"]["candidate_q"]
        > 600 * unsorted_row["spatial"]["quadratic_cutoff"]
    )
    assert "max_lattice_errors" not in unsorted_row
    for bit, row in enumerate(uniform):
        assert row["spatial"]["status"] == "pass"
        assert row["max_lattice_errors"][bit] < 0.004
        assert row["max_lattice_errors"][1 - bit] > 0.35


def test_conditional_profiles_keep_every_repair_history():
    receipt = json.loads((unsorted.HERE / "conditional_profiles.json").read_bytes())
    rows = receipt["rows"]
    level, route = single_write.enroll()
    assert receipt["level"] == level and receipt["route"] == [list(p) for p in route]
    assert [row["record"] for row in rows] == [format(i, "05b") for i in range(32)]
    assert all(len(row["count_sha256"]) == 17 for row in rows)
    assert all(row["launch_label_rejected"] for row in rows)
    passing = [r["record"] for r in rows if r["odd_profile_q"] <= r["quadratic_cutoff"]]
    assert passing == receipt["profile_pass_records"] == ["00000", "11111"]
    assert receipt["native_weight_passing_profile"] == len(passing) / len(rows)
    assert sum(r["pair_sum_closure"] < 0.0048 for r in rows[1:-1]) == 13


def test_free_seam_preserves_each_sharp_pointer_at_every_coupling():
    law = experiment.candidate.native_law
    exchange, contrast = law.local_term_parts()
    energies = np.diag([int(word != (1, 1)) for word in law.SRC])
    for bit in (0, 1):
        memory = growth.label(2, 1 - 2 * bit)
        embedding = np.kron(memory[:, None], np.eye(law.N_SRC))
        np.testing.assert_allclose(exchange @ embedding, bit * embedding @ energies, atol=1e-13)
        np.testing.assert_allclose(contrast @ embedding, 4 * embedding @ energies, atol=1e-13)
