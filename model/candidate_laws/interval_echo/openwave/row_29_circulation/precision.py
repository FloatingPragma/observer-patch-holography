"""Count a weak distant direction signal with a second acquisition and an owner-number read."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import eigsh, expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
from model.executable_model import carrier
from .experiment import HERE, CURRENT, circulation, laboratories, guide, select_source

SHOTS = 2**32
TIMES = (22, 23, 24, 25, 26)


def run(times=TIMES, seed=291125):
    owners = laboratories()[-1]
    board = guide(owners)
    face = (27, 32, 33)
    assert face in board.faces()
    source, _, _ = select_source(board)
    rng = np.random.default_rng(seed)
    panels = []
    for load, localized, cut in (
        (2, False, False),
        (2, True, False),
        (2, True, True),
        (1, False, False),
    ):
        words = list(combinations(range(board.n), load))
        index = {w: i for i, w in enumerate(words)}
        h = (
            PairSector(board).hamiltonian(candidate.KAPPA)
            if load == 2
            else board.laplacian() / 2
        )
        if load == 2:
            _, v = eigsh(h, k=2, which="SA", v0=np.ones(len(words)), tol=2e-13)
            ground = v[:, 0]
        else:
            ground = np.ones(len(words)) / np.sqrt(len(words))
        if localized:
            ground = ground * np.array([max(w) < 12 for w in words])
        success = min(float(ground @ ground), 1.0)
        src = circulation(words, index, source)
        j = circulation(words, index, face) / CURRENT
        number = np.array([sum(p // 12 == 2 for p in w) for w in words])
        launch = np.column_stack(
            [expm_multiply(-1j * s * src, ground) for s in (1, -1)]
        )
        if cut:
            islands = carrier.Carrier(
                "isolated owners",
                board.n,
                [e for e in board.seams if e[0] // 12 == e[1] // 12],
            )
            h = PairSector(islands).hamiltonian(candidate.KAPPA)
        waves = expm_multiply(
            -1j * h, launch, start=0, stop=max(times), num=max(times) + 1
        )
        reads = []
        for t in times:
            x = waves[t]
            jx = j @ x
            joint = []
            for n in (0, 1, 2):
                mask = number == n
                pop = np.sum(abs(x[mask]) ** 2, axis=0)
                mean = np.sum(x[mask].conj() * jx[mask], axis=0).real
                nz = np.sum(abs(jx[mask]) ** 2, axis=0)
                joint.extend(((nz + mean) / 2, pop - nz, (nz - mean) / 2))
            joint.append(np.full(2, 1 - success))
            prob = np.array(joint).T
            assert prob.min() > -1e-10 and np.max(abs(prob.sum(axis=1) - 1)) < 1e-10
            obs = []
            for p in prob:
                # Put the dominant category last so that roundoff on small
                # probabilities cannot consume the multinomial remainder.
                # Retain the native order among exactly equivalent outcomes.
                dominant = int(np.argmax(p))
                order = [i for i in range(len(p)) if i != dominant] + [dominant]
                sampled = rng.multinomial(SHOTS, np.clip(p[order].round(12), 0, 1))
                obs.append(sampled[np.argsort(order)].tolist())
            reads.append(
                dict(time=t, counts=obs, probabilities=prob.round(12).tolist())
            )
        panels.append(
            dict(
                load=load,
                local_herald=localized,
                bridge_cut=cut,
                success=success,
                reads=reads,
            )
        )
    return dict(
        shots=SHOTS,
        face=face,
        source=source,
        times=times,
        simultaneous_error=float(
            np.sqrt(4 * np.log(2 * len(times) * len(panels) / 0.01) / SHOTS)
        ),
        panels=panels,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = run()
    if args.output:
        with args.output.open("x") as f:
            json.dump(data, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        compare(
            json.loads(json.dumps(data)),
            json.loads((HERE / "precision.json").read_text()),
            "precision",
        )
        print(
            "Verified held two-step direction reads, heralds and co-located material counts."
        )


if __name__ == "__main__":
    main()
