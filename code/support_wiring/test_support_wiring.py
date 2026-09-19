from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from support_wiring import experiment as producer, readouts, verify
from support_wiring import check_reproduction


def test_exact_dyadic_encoding_and_overflow():
    values = np.array([0, 1, 2**63, 2**127 + 7, 2**192 - 1], dtype=object)
    assert np.array_equal(values, producer.unlimbs(producer.limbs(values)))
    with pytest.raises(AssertionError):
        producer.limbs(np.array([2**192], dtype=object))
    # The nearest-integer update used by the abandoned partial run is excluded.
    state = np.array([0, 1], dtype=object)
    numerator = state[0] + state[1]
    assert numerator == 1  # both outputs are 1 / 2, neither 0 nor 1


@pytest.mark.parametrize("level", [3, 4, 5])
def test_complete_matching_schedules_and_readback(level):
    g = dict(np.load(producer.HERE / f"geometry/geometry_l{level}.npz"))
    c = len(g["faces"])
    for wiring, expected in (("w12", 36*c-30), ("w3", 63*c//2), ("isolated", 30*c)):
        schedule = producer.phases(g, wiring)
        assert sum(len(e) for _, e in schedule) == expected
        for _, edges in schedule:
            assert len(np.unique(edges)) == 2 * len(edges)
    lap, projector, frame = producer.carrier(g["seams"])
    np.testing.assert_allclose(lap @ projector, (5-np.sqrt(5))*projector, atol=1e-13)
    np.testing.assert_allclose(frame@frame.T, 4*projector, atol=1e-13)
    assert np.linalg.matrix_rank(frame) == 3


def test_refinement_is_four_children_and_state_preserving():
    for level in (4, 5):
        g = dict(np.load(producer.HERE / f"geometry/geometry_l{level}.npz"))
        coarse = dict(np.load(producer.HERE / f"geometry/geometry_l{level-1}.npz"))
        parents, children, weights = g["parent"], g["children"], g["expectation_weights"]
        np.testing.assert_array_equal(parents[children], np.arange(len(children))[:, None].repeat(4, axis=1))
        np.testing.assert_allclose(weights[children].sum(axis=1), 1, atol=2e-15)
        np.testing.assert_allclose(g["areas"][children].sum(axis=1), coarse["areas"], atol=2e-15)
        np.testing.assert_allclose(weights, g["areas"] / coarse["areas"][parents], atol=1e-13)


@pytest.fixture
def small_trace(tmp_path, monkeypatch):
    """Small executor fixture; the full geodesic census is checked separately.

    Use two initial carriers and the same four-child copy convention. This
    exercises every serialized phase and law without copying the 60 MB archive
    for each deliberate corruption.
    """
    original = producer.HERE
    seams = np.load(original / "geometry/geometry_l3.npz")["seams"]
    (tmp_path / "geometry").mkdir()
    (tmp_path / "SPECIFICATION.md").write_bytes((original / "SPECIFICATION.md").read_bytes())
    (tmp_path / "geometry/geometry.json").write_text("{}\n")
    geometries = {}
    for level, count in ((3, 2), (4, 8), (5, 32)):
        glued = np.array([[a, p, a+1, p] for a in range(0, count, 2) for p in range(12)])
        g = {"faces": np.zeros((count, 3), dtype=int), "seams": seams, "w12": glued}
        if level > 3:
            g["parent"] = np.arange(count)//4
        np.savez_compressed(tmp_path / f"geometry/geometry_l{level}.npz", **g)
        geometries[level] = g
    monkeypatch.setattr(producer, "HERE", tmp_path)
    monkeypatch.setattr(verify, "HERE", tmp_path)
    trace = tmp_path / "trace"
    producer.build_trace(trace)
    return trace, geometries


def reseal(trace, manifest):
    """Adversary updates all hashes: semantic rejection must still hold."""
    import hashlib
    chain = hashlib.sha256(producer.canonical(manifest["binding"])).hexdigest()
    for item in manifest["chunks"]:
        item["previous"] = chain
        item["sha256"] = producer.sha(trace / item["file"])
        item.pop("chain", None)
        chain = hashlib.sha256(producer.canonical(item)).hexdigest()
        item["chain"] = chain
    manifest["final_chain"] = chain
    producer.save_json(trace / "trace.json", manifest)


def test_small_trace_independent_replay(small_trace):
    directory, g = small_trace
    replay = verify.verify_trace(directory, g)
    assert replay[-1]["exact_total_and_descent_checks"]
    # The two writes of an action read old versions, not one another.
    first = replay[0]["chunks"][1]["first_event"]
    assert np.array_equal(replay[1][first], replay[1][first+1])
    assert first not in replay[1][first+1]


@pytest.mark.parametrize("mutation, message", [
    ("mean", "noncanonical seam mean"),
    ("stale", "stale or unauthenticated read"),
    ("omit", "missing or changed seam attempts"),
    ("copy", "refinement componentwise copy"),
    ("join_writer", "refinement current parent writers"),
    ("position", "writer rank-three readback"),
    ("law", "trace contract binding"),
    ("missing_phase", "complete scheduled phase count"),
])
def test_rehashed_semantic_corruptions_fail(small_trace, mutation, message):
    import json
    trace, g = small_trace
    manifest = json.loads((trace / "trace.json").read_text())
    item = next(x for x in manifest["chunks"] if x["kind"] == ("copy" if mutation in ("copy", "join_writer") else "intra"))
    data = dict(np.load(trace / item["file"]))
    if mutation in ("mean", "copy"):
        data["values"][0, 0] ^= np.uint64(1)
    elif mutation in ("stale", "join_writer"):
        data["parents"].flat[0] = -1
    elif mutation == "omit":
        data["endpoints"] = data["endpoints"][1:]
    elif mutation == "position":
        data["positions"][0, 0] += .01
    elif mutation == "law":
        manifest["binding"]["law"] = "integer_nearest_agreement"
    elif mutation == "missing_phase":
        del manifest["chunks"][1]
    np.savez_compressed(trace / item["file"], **data)
    reseal(trace, manifest)
    with pytest.raises(ValueError, match=message):
        verify.verify_trace(trace, g)


def test_broken_hash_fails(small_trace):
    trace, g = small_trace
    path = trace / "phase_001.npz"
    with path.open("ab") as stream:
        stream.write(b"tampered")
    with pytest.raises(ValueError, match="trace chunk hash"):
        verify.verify_trace(trace, g)


def test_interval_algorithms_against_exhaustive_oracle():
    rng = np.random.default_rng(776)
    for size in range(2, 25):
        parents = [np.asarray([], dtype=int)]
        for v in range(1, size):
            parents.append(np.unique([v-1, int(rng.integers(v))]))
        for bottom in (0, size//2):
            members, row = readouts.interval(parents, bottom, size-1)
            verified, count, height = verify.count_interval(parents, bottom, size-1)
            reach = np.eye(size, dtype=bool)
            for v in range(size):
                for p in parents[v]:
                    reach[:, v] |= reach[:, p]
            expected = np.flatnonzero(reach[bottom] & reach[:, size-1])
            oracle = sum(bool(reach[a, b]) for a in expected for b in expected if a != b)
            assert np.array_equal(members, expected) and np.array_equal(verified, expected)
            assert count == row["strict_pairs"] == oracle
            verify.check_order_row(row, members, count, height)
            row["strict_pairs"] += 1
            with pytest.raises(ValueError, match="interval counts"):
                verify.check_order_row(row, members, count, height)


def test_dimension_references_and_undefined_cases():
    for fraction, dimension in ((1, 1), (.5, 2), (8/35, 3), (.1, 4)):
        assert readouts.dimension(fraction) == pytest.approx(dimension, abs=1e-12)
    assert readouts.dimension(0) is None
    assert readouts.dimension(-1) is None
    assert readouts.dimension(1.1) is None


@pytest.mark.parametrize("mutation", ["matrix", "share", "schedule"])
def test_kernel_and_confluence_mutations_fail(tmp_path, mutation):
    import json
    archive = producer.ARCHIVE
    rows = json.loads((archive/"controls.json").read_text())
    matrices = dict(np.load(archive/"kernels.npz"))
    if mutation == "matrix":
        key = rows[0]["kernels"][0]["matrix"]
        matrices[key][0, 0] += .1
        expected = "independent response Gram"
    elif mutation == "share":
        rows[0]["kernels"][0]["slow_share"] += .1
        expected = "slow-band share"
    else:
        rows[0]["confluence"][0]["squared_distance_by_sweep"][1] += 1
        expected = "finite confluence residuals"
    np.savez_compressed(tmp_path/"kernels.npz", **matrices)
    producer.save_json(tmp_path/"controls.json", rows)
    geometry = {level: dict(np.load(producer.HERE/f"geometry/geometry_l{level}.npz")) for level in (3, 4, 5)}
    with pytest.raises(ValueError, match=expected):
        verify.verify_controls(tmp_path, geometry)


@pytest.mark.parametrize("mutation", ["writer", "value", "fractional_value", "trailing_read", "offset", "interior"])
def test_q13_rehashed_mutations_fail(tmp_path, mutation):
    import json
    data = dict(np.load(producer.ARCHIVE/"q13_reads.npz"))
    report = json.loads((producer.ARCHIVE/"q13.json").read_text())
    if mutation == "writer":
        data["read_writers"][0] += 1
        expected = "q13 authenticated previous versions"
    elif mutation == "value":
        data["values"][1, 0] += 1
        expected = "q13 read law"
    elif mutation == "fractional_value":
        data["values"] = data["values"].astype(float)
        data["values"][4, 0] += .5
        expected = "q13 exact array format: values"
    elif mutation == "trailing_read":
        data["read_writers"] = np.r_[data["read_writers"], np.int32(0)]
        data["parent_offsets"] = np.r_[data["parent_offsets"], len(data["read_writers"])]
        report["reads"] += 1
        expected = "q13 exact array format"
    elif mutation == "offset":
        data["parent_offsets"][2198] += 1
        expected = "q13 complete read offsets"
    else:
        report["intervals"][-1]["interior"] = True
        expected = "q13 interior/clipping flag"
    np.savez_compressed(tmp_path/"q13_reads.npz", **data)
    report["trace_sha256"] = producer.sha(tmp_path/"q13_reads.npz")
    producer.save_json(tmp_path/"q13.json", report)
    geometry = dict(np.load(producer.HERE/"geometry/geometry_l3.npz"))
    with pytest.raises(ValueError, match=expected):
        verify.verify_q13(tmp_path, geometry)


@pytest.mark.parametrize("field", ["canonical_w3", "historical_integer_control"])
@pytest.mark.parametrize("mutation", ["empty", "omit_law", "omit_result"])
def test_missing_historical_comparison_fails(tmp_path, field, mutation):
    report = verify.read_json(producer.ARCHIVE/"support_wiring_receipt.json")
    if mutation == "empty":
        report["historical_L6"][field] = {}
    elif mutation == "omit_law":
        del report["historical_L6"][field]["law"]
    else:
        key = ("terminated" if field == "canonical_w3" else "quotient_hash_equals_expected_all")
        del report["historical_L6"][field][key]
    producer.save_json(tmp_path/"support_wiring_receipt.json", report)
    manifest = verify.read_json(producer.ARCHIVE/"trace/trace.json")
    # Exercise the receipt's comparison boundary using the already archived
    # replay outputs. Full replay is independently executed by CI.
    summary = {"mean_actions": report["execution"]["mean_actions"],
               "final_exponent": report["execution"]["final_denominator_exponent"],
               "final_chain": manifest["final_chain"]}
    with pytest.raises(ValueError, match="complete historical law-labelled comparison"):
        verify.verify_receipt(tmp_path, report["wiring"], (manifest, summary),
                              report["provenance_intervals"], report["record_metric_q13"],
                              report["canonical_controls"], report["record_metric_q13_controls"])


def test_q13_control_tapes_replay_and_reproduce(tmp_path):
    controls = verify.read_json(producer.ARCHIVE/"q13_controls.json")
    geometry = dict(np.load(producer.HERE/"geometry/geometry_l3.npz"))
    for dim, saved in enumerate(controls["families"], 1):
        fresh = readouts.record_metric_family(tmp_path, dim)
        verify.verify_metric_family(tmp_path, geometry, fresh, dim)
        filename = f"q13_control_d{dim}_reads.npz"
        check_reproduction.arrays(producer.ARCHIVE/filename, tmp_path/filename)
        fresh.pop("trace_sha256")
        saved.pop("trace_sha256")
        check_reproduction.same(saved, fresh)


@pytest.mark.parametrize("mutation, message", [
    ("missing", "complete q13 control families"),
    ("duplicate", "complete q13 control families"),
    ("target", "q13 control comparison boundary"),
    ("count", "interval counts"),
    ("historical", "q13 historical interval counts"),
    ("clipping", "q13 interior/clipping flag"),
])
def test_q13_control_mutations_fail(tmp_path, mutation, message):
    import shutil
    controls = verify.read_json(producer.ARCHIVE/"q13_controls.json")
    if mutation == "missing":
        controls["families"].pop()
    elif mutation == "duplicate":
        controls["families"][1] = controls["families"][0]
    elif mutation == "target":
        controls["reference_is_acceptance_target"] = True
    else:
        row = controls["families"][0]["intervals"][-1]
        if mutation == "count":
            row["strict_pairs"] += 1
        elif mutation == "historical":
            row["historical_counts_equal"] = False
        else:
            row["interior"] = True
    for dim in (1, 2):
        name = f"q13_control_d{dim}_reads.npz"
        shutil.copyfile(producer.ARCHIVE/name, tmp_path/name)
    producer.save_json(tmp_path/"q13_controls.json", controls)
    geometry = dict(np.load(producer.HERE/"geometry/geometry_l3.npz"))
    with pytest.raises(ValueError, match=message):
        verify.verify_q13_controls(tmp_path, geometry)


@pytest.mark.parametrize("raw, message", [('{'+'"x":1,"x":2}', "duplicate JSON key"),
                                          ('{'+'"x":NaN}', "nonfinite JSON constant")])
def test_strict_json_rejects_ambiguous_custody(tmp_path, raw, message):
    path = tmp_path/"ambiguous.json"
    path.write_text(raw)
    with pytest.raises(ValueError, match=message):
        verify.read_json(path)


def test_reproduction_keeps_integer_evidence_exact():
    check_reproduction.same({"count": 2**60, "readback": .1}, {"count": 2**60, "readback": .1+1e-12})
    with pytest.raises(ValueError, match="exact field"):
        check_reproduction.same(2**60, 2**60+1)
    with pytest.raises(ValueError, match="float"):
        check_reproduction.same(float("nan"), float("nan"))
    with pytest.raises(ValueError, match="keys"):
        check_reproduction.same({"x": 1}, {})


def test_reproduction_checks_all_array_payloads(tmp_path):
    a, b = tmp_path/"a.npz", tmp_path/"b.npz"
    np.savez(a, counts=np.array([2**60], dtype=np.uint64))
    np.savez(b, counts=np.array([2**60+1], dtype=np.uint64))
    with pytest.raises(ValueError, match="array values"):
        check_reproduction.arrays(a, b)
