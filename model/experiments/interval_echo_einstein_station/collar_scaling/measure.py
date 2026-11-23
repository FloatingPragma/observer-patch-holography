"""Spatial scaling and locality of a native two-word connected modular charge."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra, shortest_path

from model.executable_model import cell, federation
from model.executable_model.refinement import load_refinement
from model.executable_model.symmetry import rotations
from .. import shared_variance, volume
from ..proper_ball_ruler.measure import pulse_calibration
from .replay import verify_replay

HERE = Path(__file__).resolve().parent
SEED, SHOTS, REFERENCE_SHOTS = 762251, 2**48, 2**38


class ProtectedRanges:
    """Contract each protected cell template to three native gateway stations."""

    def __init__(self, board):
        self.gates = np.sort(board.external.ravel()).reshape(board.cells, 3)
        if np.any(self.gates % 12 == 4):
            raise ValueError("A parked source centre serves an external gateway")
        graph = cell().adjacency().toarray()
        graph[4, :] = graph[:, 4] = 0
        self.local = shortest_path(graph, directed=False, unweighted=True)
        a, b = np.array(list(combinations(range(3), 2))).T
        ids = np.arange(3 * board.cells).reshape(board.cells, 3)
        x, y = ids[:, a].ravel(), ids[:, b].ravel()
        w = self.local[self.gates.ravel()[x] % 12, self.gates.ravel()[y] % 12]
        u, v = np.searchsorted(self.gates.ravel(), board.external).T
        self.routes = coo_matrix(
            (np.r_[w, w, np.ones(len(u) * 2)], (np.r_[x, y, u, v], np.r_[y, x, v, u])),
            shape=(3 * board.cells,) * 2,
        ).tocsr()

    def survey(self, source_cell):
        seed = 3 * source_cell + np.arange(3)
        d = dijkstra(self.routes, indices=seed)
        depart = self.local[0, self.gates[source_cell] % 12]
        reach = np.min(d + depart[:, None], axis=0).reshape(-1, 3)
        result = np.min(reach + self.local[self.gates % 12, 0], axis=1)
        result[source_cell] = 0
        if not np.all(np.isfinite(result)):
            raise ValueError("Protected clock stations are disconnected")
        return result.astype(int)


def measure_support(directory, level):
    board = load_refinement(directory / f"FEDERATION_LEVEL{level}.npz", level=level)
    ranges = ProtectedRanges(board)
    a = ranges.survey(0)
    distal = int(np.argmax(a))
    b = ranges.survey(distal)
    radius = 6 * 2 ** ((level - 1) // 2)
    classes = []
    for scale in (1, 2):
        inner, outer = radius * scale / 2, radius * scale
        aa = (a > inner) & (a <= outer)
        bb = (b > inner) & (b <= outer)
        if np.any(aa & bb):
            raise ValueError("Experimental annuli overlap")
        classes.append([int(aa.sum()), int(bb.sum()), int((~(aa | bb)).sum())])
    return dict(
        level=level,
        cells=board.cells,
        board_sha256=board.sha256,
        home=0,
        distal=distal,
        calibration_range=int(a.max()),
        radius=radius,
        annulus_counts=classes,
    )


def slot_membership(counts, colocated):
    n = np.asarray(counts, dtype=float)
    if colocated:
        return np.diag(n / n.sum())
    return (np.outer(n, n) - np.diag(n)) / (n.sum() * (n.sum() - 1))


def acquire(law, counts, colocated, rng):
    placement = slot_membership(counts, colocated)
    trials = rng.multinomial(SHOTS, placement.ravel()).reshape(3, 3)
    histograms = {}
    for a, b in np.ndindex(3, 3):
        if a == 2 and b == 2:
            histogram = [int(trials[a, b])]
        else:
            distribution = (
                law.ravel() if max(a, b) < 2 else law.sum(axis=0 if a == 2 else 1)
            )
            histogram = rng.multinomial(trials[a, b], distribution).tolist()
        histograms[f"{a}{b}"] = histogram
    return histograms


def modular_score(reference):
    return np.log(np.outer(reference.sum(axis=1), reference.sum(axis=0)) / reference)


def summarize(
    shared, fresh, reversed_arm, counts, colocated, score, conditional, relative_radius
):
    placement = slot_membership(counts, colocated)
    result = {}
    for label, entries in (
        ("home", ((0, 0),)),
        ("distal", ((1, 1),)),
        ("cross", ((0, 1), (1, 0))),
    ):
        chance = sum(placement[x, y] for x, y in entries)

        def moment(arm):
            return (
                sum(np.asarray(arm[f"{x}{y}"]) @ score.ravel() for x, y in entries)
                / SHOTS
            )

        measured = float(moment(shared) - moment(fresh))
        reverse = float(moment(reversed_arm) - moment(fresh))
        # Conservative Bernstein union over 96 contrasts, conditional on the independently
        # frozen score. Each outcome is zero or bounded by max|score| on an admission.
        m = float(np.max(abs(score)))
        log_tail = np.log(2 * 384 / 0.01)
        half_width = (
            (
                2 * np.sqrt(2 * chance * m**2 * log_tail / SHOTS)
                + 8 * m * log_tail / (3 * SHOTS)
            )
            if chance
            else 0.0
        )
        exact = float(chance * conditional)
        result[label] = dict(
            admission=float(chance),
            population_shared_change=exact,
            measured_shared_change=measured,
            measured_reverse_change=reverse,
            simultaneous_radius=float(half_width),
            normalized_population=exact / relative_radius**4,
            normalized_measurement=measured / relative_radius**4,
        )
    return result


def source_reader():
    _, _, _, source = volume.radar_source_and_hand(federation(0))
    # A zero-energy written mode and a blank seam are orthogonal native outcomes.
    # The presence projector therefore implements the admission flag without a
    # separately carried metadata bit. It commutes with the native seam generator.
    memory_modes = np.column_stack(
        [volume.label(0, 1), volume.label(2, 1), volume.label(2, -1)]
    )
    np.testing.assert_array_equal(memory_modes[0], np.zeros(3))
    np.testing.assert_allclose(memory_modes.T @ memory_modes, np.eye(3), atol=2e-15)
    native = volume.candidate.native_law
    blank = np.zeros_like(native.seam_term(volume.candidate.KAPPA))
    blank[: len(native.SRC), : len(native.SRC)] = np.eye(len(native.SRC))
    generator = native.seam_term(volume.candidate.KAPPA)
    np.testing.assert_array_equal(blank @ generator, generator @ blank)
    laws, _ = shared_variance.source_words(source)
    ref_rng = np.random.default_rng([SEED, 0])
    counts = ref_rng.multinomial(REFERENCE_SHOTS, laws[0].ravel()).reshape(9, 9)
    score = modular_score(counts / counts.sum())
    delta = laws[2] - laws[1]
    return (
        laws,
        score,
        dict(
            counts=counts.tolist(),
            score_when_both_shared_change=float((score * delta).sum()),
            ideal_reference_contrasts=[
                float((modular_score(w) * delta).sum()) for w in laws
            ],
        ),
    )


def experiment(supports):
    laws, score, reference = source_reader()
    rotation = min(r for r in rotations() if r[1] == 4)
    assert all(
        4 in (rotation[a], rotation[b]) for a, b in ((0, 1), (1, 5), (1, 7), (1, 8))
    )
    panels = []
    for j, support in enumerate(supports):
        result = []
        for k, n in enumerate(support["annulus_counts"]):
            rho = (k + 1) * support["radius"] / support["calibration_range"]
            arm = {}
            for pos, colocated in enumerate((False, True)):
                readings = [
                    acquire(
                        law,
                        n,
                        colocated,
                        np.random.default_rng([SEED, 1, j, k, pos, i]),
                    )
                    for i, law in enumerate(laws)
                ]
                contrasts = summarize(
                    readings[2],
                    readings[1],
                    readings[0],
                    n,
                    colocated,
                    score,
                    reference["score_when_both_shared_change"],
                    rho,
                )
                arm["colocated" if colocated else "separated"] = dict(
                    reads=readings, contrasts=contrasts
                )
            result.append(dict(relative_radius=rho, arms=arm))
        panels.append(dict(geometry=support, settings=result))
    return dict(
        protocol="protected cell clock, native protected three-level writer",
        seed=SEED,
        shots=SHOTS,
        reference_shots=REFERENCE_SHOTS,
        decoded_swap=pulse_calibration()["decoded_single_step"],
        source_rotation=list(rotation),
        reference=reference,
        panels=panels,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--directory",
        type=Path,
        help="Recompute route geometry from the pinned carrier files",
    )
    parser.add_argument("--output", type=Path, help="Write a fresh acquisition")
    args = parser.parse_args()
    prior = HERE / "results.json"
    if args.directory:
        supports = [measure_support(args.directory, level) for level in (1, 3, 5, 7)]
    else:
        supports = [
            item["geometry"] for item in json.loads(prior.read_text())["panels"]
        ]
    report = experiment(supports)
    if args.output:
        with args.output.open("x") as stream:
            json.dump(report, stream, separators=(",", ":"), allow_nan=False)
            stream.write("\n")
    else:
        verify_replay(report, json.loads(prior.read_text()), name="reads")
    for panel in report["panels"]:
        geometry = panel["geometry"]
        print(
            geometry["cells"], geometry["annulus_counts"], geometry["calibration_range"]
        )
        for row in panel["settings"]:
            for mode, arm in row["arms"].items():
                hh = arm["contrasts"]["home"]
                print(
                    mode,
                    row["relative_radius"],
                    hh["normalized_population"],
                    hh["measured_shared_change"],
                    hh["simultaneous_radius"],
                    arm["contrasts"]["cross"]["normalized_population"],
                )


if __name__ == "__main__":
    main()
