"""Physical discrimination of blank and vacuum-written memory on a known material input."""

import argparse
from collections import deque
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm

from model.experiments.count_intervals import clopper_pearson
from .. import measure as release

HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "results.json.gz"
SHOTS, SEED = 2**22, 7311931


def detector(time):
    law = release.candidate.native_law
    echo = law.echo()
    basis = np.eye(law.N_SRC)
    host = basis[:, law.SRC.index((0, 0))]
    written = (echo @ np.kron(np.eye(11)[0], host)).reshape(11, 9)[:, 0]
    frames = np.column_stack(
        [
            np.kron(mem, basis[:, law.SRC.index((1, 0))])
            for mem in (np.eye(11)[0], written)
        ]
    )
    h0, h1 = law.local_term_parts()
    final = expm(-1j * time * (h0 + release.candidate.KAPPA * h1)) @ frames
    blocks = final.reshape(11, 9, 2)
    effects = np.asarray([b.conj().T @ b for b in blocks.transpose(1, 0, 2)])
    release.check(
        np.linalg.norm(effects.sum(axis=0) - np.eye(2)) < 2e-12,
        "incomplete public material read",
    )
    correct = np.kron(np.eye(11)[0], basis[:, law.SRC.index((0, 1))])
    if time == np.pi:
        release.check(
            np.linalg.norm(final[:, 0] - correct) < 3e-12
            and np.linalg.norm(final[:, 1] - frames[:, 1]) < 3e-12,
            "native contact failed to transfer the flag to public material",
        )
    release.check(
        np.max(abs(effects[:, 0, 1])) < 2e-13,
        "memory coherence affects terminal word probabilities",
    )
    return effects


def paths(board, seam):
    dest = seam[0]
    parent = {dest: None}
    todo = deque([dest])
    while todo:
        p = todo.popleft()
        for q in sorted(board.neighbours[p]):
            if q not in parent and {p, q} != set(seam):
                parent[q] = p
                todo.append(q)
    release.check(
        len(parent) == board.n,
        "an orbital address cannot approach without crossing the core",
    )
    return [parent[i] for i in range(board.n)]


def observe(a, b, refusal, head, rng):
    # After the sharp position read each column has a two-component memory
    # state. Every subsequent operation uses only that public address.
    p = np.array(
        [np.real(e[0, 0]) * abs(a) ** 2 + np.real(e[1, 1]) * abs(b) ** 2 for e in head]
    ).T
    joint = np.r_[p.ravel(), refusal]
    release.check(
        np.min(joint) >= 0 and abs(joint.sum() - 1) < 3e-11,
        "position/head read lost a trial",
    )
    raw = release.sample(joint, rng)
    counts = np.asarray(raw["counts"])
    table = counts[:-1].reshape(len(a), 9)
    accepted = int(table.sum())
    flips = int(table[:, 1].sum())  # local word (0,1)
    lo, hi = clopper_pearson(flips, accepted, 0.01 / (2 * 2 * 2 * 4 * 2 * 3))
    return dict(
        flip_count=flips,
        admitted=accepted,
        spectral_refusal=int(counts[-1]),
        flip_probability_given_admission=float(p[:, 1].sum() / (1 - refusal)),
        simultaneous_99_interval=[float(lo), float(hi)],
        sampling_grid_l1=raw["sampling_grid_l1"],
        measured_words_per_position=[
            [int(i), int(w), int(table[i, w])] for i, w in zip(*np.nonzero(table))
        ],
    )


def report():
    release.check(release.SHOTS == SHOTS, "source enrollment changed")
    heads = {
        name: detector(t)
        for name, t in [
            ("native_pi", np.pi),
            ("half_time", np.pi / 2),
            ("omitted", 0.0),
        ]
    }
    rng = np.random.default_rng(SEED)
    boards = []
    for level in (0, 1):
        board = release.federation(level)
        seam = release.EMITTER_SEAM[level]
        h = release.transition.generator(board, seam, release.candidate.KAPPA)
        psi, refusal, preparation = release.prepare(h, seam[0])
        blank_spec = eigh(board.laplacian().toarray() / 2)
        panels = []
        for contrast in (8, 0):
            material = eigh(release.transition.generator(board, seam, contrast))
            for action in release.OPERATIONS:
                blank, written = release.branch_states(psi, seam, action)
                for t in (0, 64):
                    a = release.evolve(blank_spec, blank, t)
                    b = release.evolve(material, written, t)
                    records = {
                        name: observe(a, b, refusal, head, rng)
                        for name, head in heads.items()
                    }
                    panels.append(
                        dict(
                            flight_contrast=contrast,
                            operation=action,
                            delay=t,
                            measurements=records,
                        )
                    )
        boards.append(
            dict(
                cells=board.cells,
                ports=board.n,
                seam=list(seam),
                next_hop_to_meter_after_public_address=paths(board, seam),
                admission=preparation["acceptance"],
                panels=panels,
            )
        )
    return dict(
        physical_calibration={
            name: np.real(np.diagonal(h, axis1=1, axis2=2)).tolist()
            for name, h in heads.items()
        },
        shots=SHOTS,
        seed=SEED,
        sampling_weight_denominator=release.GRID,
        read_order="source flag, orbital address i, then public two-port load (a,b), encoded by 3a+b; unlisted counts zero",
        timing="source free contrast labelled per arm; isolated routing and terminal head use registered kappa=8",
        boards=boards,
    )


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path)
    args = cli.parse_args()
    result = report()
    if args.output:
        content = (json.dumps(result, indent=2, allow_nan=False) + "\n").encode()
        with args.output.open("xb") as stream:
            stream.write(
                gzip.compress(content, mtime=0)
                if args.output.suffix == ".gz"
                else content
            )
    else:
        release.transition.compare(
            result,
            json.loads(gzip.decompress(RECEIPT.read_bytes())),
            "public core readback",
        )
        print(
            "Replayed the material erasure detector on both orbital boards and all four interventions."
        )


if __name__ == "__main__":
    main()
