"""Measure how many oscillatory response patterns independently prepared charge probes recover."""

import argparse
from collections import Counter
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"

import networkx as nx
import numpy as np
from scipy.linalg import eigh, expm
from scipy.optimize import minimize_scalar
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent / "electron_field"
REPO = HERE.parents[4]
TIMES = (0, 2, 8, 32, 128)
SHOTS = 65536
SEED = 321478
SOURCE_SLOTS = (0, 3, 6, 9)


def checked_archive(root):
    manifest = {}
    for line in (ARCHIVE / "records.sha256").read_text().splitlines():
        digest, name = line.split("  ", 1)
        if name in manifest or Path(name).is_absolute() or ".." in Path(name).parts:
            raise ValueError(name)
        manifest[name] = digest
    with ZipFile(ARCHIVE / "records.zip") as archive:
        if sorted(archive.namelist()) != sorted(manifest):
            raise ValueError("Archive and manifest disagree")
        for name, digest in manifest.items():
            if sha256(archive.read(name)).hexdigest() != digest:
                raise ValueError(name)
        archive.extractall(root)
    (root / "reverse-engineering-reality").symlink_to(REPO, target_is_directory=True)
    sys.path.insert(0, str(root / "experiments" / "interval_echo_electron_field"))
    return manifest


def station_matrix(sphere, instrument):
    from geometry import stations
    from release import instruments

    graph = nx.Graph((i, j) for i, js in enumerate(sphere.neighbours) for j in js)
    distances = dict(nx.all_pairs_shortest_path_length(graph))
    selected = [0]
    sources, receivers = {0: sphere.source}, {}
    skipped = 0
    for slot in range(1, 10):
        adding_source = slot in SOURCE_SLOTS
        choices = sorted(set(graph).difference(selected),
                         key=lambda f: (-min(distances[f][i] for i in selected), f))
        for f in choices:
            link = sphere.detectors[f]
            opposite = receivers if adding_source else sources
            pairs = [(link, e) if adding_source else (e, link) for e in opposite.values()]
            if any(a[0] // 12 == b[0] // 12 for a, b in pairs):
                skipped += 1
                continue
            arms = [stations(sphere.board, sphere.faces, a, b) for a, b in pairs]
            if any(arm is None for arm in arms):
                skipped += 1
                continue
            for (a, b), arm in zip(pairs, arms):
                parked = {x for _, r, m in arm.values() for x in (r, m) if x is not None}
                assert parked.isdisjoint(sphere.core)
                instruments(instrument, a, b, arm)
            (sources if adding_source else receivers)[f] = link
            selected.append(f)
            break
        else:
            raise RuntimeError("No reachable independent station for slot " + str(slot))
    return sources, receivers, skipped, distances


def verify_source_cuts(sphere, original, sources):
    """Check circulation on each complete physical head word after reversing each chosen edge."""
    from geometry import edge, segments, support

    bits = {e: 1 << i for i, e in enumerate(sorted(sphere.head))}
    masks = [sum(bits[e] for e in support(face)) for face in sphere.faces]
    forward = [sum(bits[edge(a, b)] for a, b in segments(face) if a < b)
               for face in sphere.faces]
    base = sum(bits[e] * value for e, value in sphere.head.items())
    for f, e in sources.items():
        assert [j for j, mask in enumerate(masks) if mask & bits[e]] == [f]
    for i, w in enumerate(original.words):
        physical_word = base
        for f in np.flatnonzero(original.flipped[i]):
            physical_word ^= masks[f]
        targets = original.h.indices[original.h.indptr[i]:original.h.indptr[i + 1]]
        outgoing = {(w ^ original.words[j]).bit_length() - 1 for j in targets}
        for f, e in sources.items():
            moved = physical_word ^ bits[e]
            observed = {j for j, (mask, direction) in enumerate(zip(masks, forward))
                        if moved & mask in (direction, mask ^ direction)}
            assert observed == outgoing.difference({f})


def rank(singular):
    return int(np.count_nonzero(singular > max(1e-12, singular[0] * 1e-10)))


def resolved_poles(h, source_z, detector_z):
    size = h.shape[0]
    dense = h.toarray()
    subset = None if size < 100 else [0, 24]
    values, vectors = eigh(dense, subset_by_index=subset, driver="evr")
    residual = float(np.max(np.linalg.norm(h @ vectors - vectors * values, axis=0)))
    assert residual < 1e-9 and abs(values[0]) < 1e-12 and values[1] > 1e-8
    groups = np.split(np.arange(1, len(values)), np.flatnonzero(np.diff(values[1:]) > 1e-8) + 1)
    guard = None
    if subset:
        guard = float(values[groups[-1][0]])
        groups = groups[:-1]
    poles = []
    strongest = None
    for indices in groups:
        v = vectors[:, indices]
        residue = (detector_z.T @ v) @ (v.T @ source_z) / (2 * size)
        singular = np.linalg.svd(residue, compute_uv=False)
        if strongest is None or singular[0] > strongest[0]:
            strongest = (singular[0], residue, float(values[indices[0]]))
        poles.append({"angular_frequency": float(np.mean(values[indices])),
                      "eigenvalue_multiplicity": len(indices),
                      "singular_values": singular.tolist(), "response_rank": rank(singular)})
    assert strongest is not None
    # An independently addressed pair at EXACTLY the same frequency must register twice.
    tagged = np.kron(np.eye(2), strongest[1])
    split = np.linalg.svd(tagged, compute_uv=False)
    blind = np.linalg.svd(tagged[:detector_z.shape[1], :] + tagged[detector_z.shape[1]:, :],
                          compute_uv=False)
    assert rank(split) == 2 * rank(blind) and rank(blind) > 0
    return {"poles": poles, "rank_histogram": dict(sorted(Counter(p["response_rank"] for p in poles).items())),
            "max_eigenpair_residual": residual, "first_excluded_guard_frequency": guard,
            "supplied_two_copy_control": {"frequency": strongest[2], "tag_resolved_rank": rank(split),
                                           "tag_summed_rank": rank(blind)}}


def joint_record(source_sign, ei, uj, receiver, probe_sign):
    """Return joint probabilities with source first, preserving every material flag."""
    d = receiver if probe_sign == 1 else 1 - receiver
    size = len(ei)
    if source_sign == 0:
        stay_density = np.ones(size) / size
        move_density = np.zeros(size)
    else:
        # Moved words see the source face blocked. Their selector is constant on cut components.
        permitted = ei if source_sign == 1 else 1 - ei
        unmoved_amplitude = 1 - uj if source_sign == 1 else uj
        stay_density = abs(unmoved_amplitude)**2 / size
        move_density = permitted / size
    probs = np.array([stay_density @ (1 - d), stay_density @ d,
                      move_density @ (1 - d), move_density @ d])
    np.testing.assert_allclose(probs.sum(), 1, atol=2e-11, rtol=0)
    return probs


def count_spectroscopy(h, x, d, poles):
    """Fit the lowest frequency and its unrestricted 6x4 residue using material counts."""
    rng = np.random.default_rng(757472)
    times = np.sort(rng.uniform(0, 4096, size=256))
    field_matrix = h.toarray()
    frequencies = np.array([p["angular_frequency"] for p in poles])
    coefficients = 1 - np.cos(np.outer(times, frequencies))
    fit_matrix = np.linalg.pinv(coefficients)
    shape = (len(times), d.shape[1], x.shape[1])
    ideal_data, count_data = np.zeros(shape), np.zeros(shape)
    all_counts = []
    for ti, t in enumerate(times):
        moved_free = expm((-1j * t) * field_matrix) @ x
        batch = []
        for di in range(d.shape[1]):
            for si in range(x.shape[1]):
                signed_counts = []
                for q in (1, -1):
                    probabilities = joint_record(q, x[:, si], moved_free[:, si], d[:, di], 1)
                    counts = rng.multinomial(SHOTS, probabilities / probabilities.sum())
                    signed_counts.append(counts.tolist())
                    assert counts.sum() == SHOTS
                    ideal_data[ti, di, si] += q * (probabilities[1] + probabilities[3])
                    count_data[ti, di, si] += q * (counts[1] + counts[3]) / SHOTS
                batch.append(signed_counts)
        all_counts.append(batch)
    weights = fit_matrix[0]
    ideal = np.tensordot(weights, ideal_data, axes=(0, 0))
    empirical = np.tensordot(weights, count_data, axes=(0, 0))
    exact_s = np.linalg.svd(ideal, compute_uv=False)
    empirical_s = np.linalg.svd(empirical, compute_uv=False)
    np.testing.assert_allclose(exact_s, poles[0]["singular_values"], rtol=0, atol=2e-11)

    # When higher rates are independently calibrated, eliminate their full unconstrained residue matrices.
    nuisance, _ = np.linalg.qr(coefficients[:, 1:], mode="reduced")
    measurements = count_data.reshape(len(times), -1)
    residual = measurements - nuisance @ (nuisance.T @ measurements)
    def loss(omega):
        shape = 1 - np.cos(times * omega)
        shape -= nuisance @ (nuisance.T @ shape)
        projection = np.outer(shape, shape @ residual) / (shape @ shape)
        return np.linalg.norm(residual - projection)**2
    trial = np.linspace(0.001, 0.008, 141)
    lowest = int(np.argmin([loss(v) for v in trial]))
    assert 0 < lowest < len(trial) - 1
    fit = minimize_scalar(loss, bounds=(trial[lowest - 1], trial[lowest + 1]),
                          method="bounded", options={"xatol": 1e-11})
    assert fit.success
    # Two independent binomial proportions per contrast, with every enrolled exposure kept.
    alpha = .01
    cell_bound = np.linalg.norm(weights) * np.sqrt(np.log(2 * empirical.size / alpha) / SHOTS)
    operator_bound = np.sqrt(empirical.size) * cell_bound
    lower_first = max(0., empirical_s[0] - operator_bound)
    upper_second = empirical_s[1] + operator_bound
    return {"time_window": [0, 4096], "seed": 757472, "arrival_times": times.tolist(),
            "enrolled_shots": len(times) * empirical.size * 2 * SHOTS,
            "counts_axes": ["arrival_time", "receiver_major_source_minor_pair", "source_sign_+1_then_-1",
                            "four_joint_outcomes_in_receipt_order"],
            "counts": all_counts, "calibrated_lowest_frequency": frequencies[0],
            "profile_frequency_using_counts": float(fit.x), "frequency_search_interval": [0.001, 0.008],
            "simultaneous_design_condition_number": float(np.linalg.cond(coefficients)),
            "fitted_lowest_residue": empirical.tolist(),
            "fitted_lowest_singular_values": empirical_s.tolist(), "ideal_singular_values": exact_s.tolist(),
            "familywise_confidence": 1 - alpha, "entrywise_error_bound": float(cell_bound),
            "operator_error_bound": float(operator_bound), "lower_bound_first_component": float(lower_first),
            "upper_bound_second_component": float(upper_second),
            "upper_bound_second_over_first": float(upper_second / lower_first) if lower_first > 0 else None}


def board(level, instrument, rng):
    from flight import LaplacianSector, Sector, Sphere
    from waves import federation

    sphere = Sphere(federation(level), "depth-first within shell")
    original = Sector(sphere)
    edges = original.verify_native_moves(sphere)
    sector = LaplacianSector(original)
    sources, detectors, skipped, dist = station_matrix(sphere, instrument)
    verify_source_cuts(sphere, original, sources)
    def heads(stations):
        return np.array([sector.flipped[:, f] ^ bool(sphere.head[e]) for f, e in stations.items()],
                        dtype=float).T
    x, d = heads(sources), heads(detectors)
    z, y = 2 * x - 1, 2 * d - 1
    h = sector.h
    size = len(sector.words)
    np.testing.assert_allclose(h @ np.ones(size), 0, atol=1e-12)
    for source_index, f in enumerate(sources):
        np.testing.assert_allclose(sector.hamiltonian(cut=(f,)) @ z[:, source_index], 0, atol=1e-12)
    spectrum = resolved_poles(h, z, y)
    blocked = {}
    batches = []
    max_closure = 0.0
    for t in TIMES:
        u = expm_multiply((-1j * t) * h, x) if t else x.copy()
        contrast = (y.T @ (z - (2 * u.real - 1))) / (2 * size)
        local = []
        for sj, f in enumerate(sources):
            for dj, g in enumerate(detectors):
                by_sign = {}
                for q in (1, 0, -1):
                    for probe in (1, -1):
                        p = joint_record(q, x[:, sj], u[:, sj], d[:, dj], probe)
                        counts = rng.multinomial(SHOTS, p / p.sum())
                        assert counts.sum() == SHOTS
                        local.append([sj, dj, q, probe, p.tolist(), counts.tolist()])
                        by_sign[q, probe] = p[1] + p[3]
                for q in (1, 0, -1):
                    np.testing.assert_allclose(by_sign[q, 1] + by_sign[q, -1], 1, atol=2e-11)
                signed = by_sign[1, 1] - by_sign[-1, 1]
                max_closure = max(max_closure, abs(signed - contrast[dj, sj]))
                if t == 0:
                    for q in (1, -1):
                        np.testing.assert_allclose(by_sign[q, 1], by_sign[0, 1], atol=2e-11)
                    blocked[str((f, g))] = abs(signed)
        batches.append({"field_time": t, "peak_absolute_signed_contrast": float(np.max(abs(contrast))),
                        "settings": local})
    assert max_closure < 2e-11 and max(blocked.values()) < 2e-11
    count_fit = {"count_spectroscopy": count_spectroscopy(h, x, d, spectrum["poles"])} if level == 1 else {}
    return {**count_fit, "cells": sphere.board.cells, "reachable_words": size, "directed_transitions": edges,
            "exterior": sphere.exterior_check(), "skipped_station_proposals": skipped,
            "source_faces_and_seams": list(sources.items()), "receiver_faces_and_seams": list(detectors.items()),
            "source_receiver_face_distances": [[dist[f][g] for g in detectors] for f in sources],
            "max_record_contrast_reconstruction_error": float(max_closure),
            "max_free_flight_disabled_contrast_all_times": float(max(blocked.values())),
            "spectrum": spectrum, "record_batches": batches}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    with TemporaryDirectory() as temp:
        root = Path(temp)
        manifest = checked_archive(root)
        from waves import Instrument, board_test, pulse_oracle
        instrument = Instrument()
        pulse_oracle(instrument)
        folder = root / "experiments" / "interval_echo_electron_field"
        saved = json.loads((folder / "waves.json").read_text())
        rng = np.random.default_rng(SEED)
        boards, guided = [], []
        for level in (1, 2):
            print(f"Resolving complete sphere at refinement {level}", flush=True)
            boards.append(board(level, instrument, rng))
            record = board_test(level, instrument)
            # Compare the actual crossing vectors and calibrated rates of all archived harmonics.
            original = next(x for x in saved["boards"] if x["cells"] == record["cells"])
            assert record["faces_on_path"] == original["faces_on_path"]
            assert [list(e) for e in record["underlying_material_seams"]] == original["underlying_material_seams"]
            np.testing.assert_allclose(record["frequency_second_over_first"], original["frequency_second_over_first"])
            for mode, expected in zip(record["modes"], original["modes"], strict=True):
                for key in mode:
                    np.testing.assert_allclose(mode[key], expected[key], rtol=2e-11, atol=2e-12)
            guided.append({"cells": record["cells"], "words": record["independent_gauss_words"],
                           "frequencies": [m["recovered_angular_frequency"] for m in record["modes"]],
                           "frequency_second_over_first": record["frequency_second_over_first"]})
        receipt = {"apparatus": "archived positive ring-Laplacian written-field comparison",
                   "units": "angular frequency per finite Hamiltonian time unit; wave number per guided face",
                   "archive_manifest_sha256": sha256((ARCHIVE / "records.sha256").read_bytes()).hexdigest(),
                   "verified_archived_files": len(manifest),
                   "design_sha256": sha256((HERE / "DESIGN.md").read_bytes()).hexdigest(),
                   "code_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
                   "shots_per_setting": SHOTS, "seed": SEED,
                   "settings_columns": ["source_index", "receiver_index", "source_sign", "probe_sign",
                                        "joint_probabilities", "joint_counts"],
                   "joint_outcomes": ["source_stays/probe_stays", "source_stays/probe_crosses",
                                      "source_moves/probe_stays", "source_moves/probe_crosses"],
                   "guided_dispersion_replay": guided, "boards": boards}
        with args.output.open("x") as out:
            json.dump(receipt, out, separators=(",", ":"), allow_nan=False)
            out.write("\n")


if __name__ == "__main__":
    main()
