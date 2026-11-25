"""Six-state transmission through the native one-load polarization channel."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.stats import beta

from .. import experiment as dipole
from .. import survey

HERE = Path(__file__).resolve().parent
SHOTS = 2**24
TIMES = dipole.TIMES
ARMS = ("native", "dephased", "half_exchange")
TAIL = 0.01 / (2 * 3 * len(TIMES) * (len(ARMS) * 6 * 2 + 1))
INPUTS = np.array([[1, 0], [0, 1], [1, 1], [1, -1], [1, 1j], [1, -1j]], complex)
INPUTS /= np.linalg.norm(INPUTS, axis=1)[:, None]


def bound(n, total):
    return [
        float(beta.ppf(TAIL, n, total - n + 1)) if n else 0.0,
        float(beta.isf(TAIL, n + 1, total - n)) if n < total else 1.0,
    ]


def summary(rows):
    counts = np.array(rows)
    if counts.shape != (6, 5) or not np.all(counts.sum(axis=1) == SHOTS):
        raise ValueError("incomplete challenge outcomes")
    c, w = np.mean(counts[:, :2], axis=0) / SHOTS
    ci, wi = [np.mean([bound(n, SHOTS) for n in counts[:, i]], axis=0) for i in (0, 1)]
    return dict(
        score=round(c - 2 * w, 10),
        score_interval=np.round(ci - 2 * wi[::-1], 10).tolist(),
        conditional_fidelity=round(c / (c + w), 7) if c + w else None,
        transverse_acceptance=round(c + w, 10),
    )


def draw(rng, probabilities):
    probs = np.asarray(probabilities)
    if probs.min() < -1e-11 or abs(sum(probs) - 1) > 1e-11:
        raise ValueError("incomplete receiver")
    p = np.maximum(probs, 0).round(12)
    return rng.multinomial(SHOTS, p / p.sum()).tolist()


def detector(cell, analyzer, frame, source, axis, dephased):
    rho = np.outer(source, source.conj())
    if dephased:
        rho = np.diag(np.diag(rho))
    local = cell @ rho @ cell.conj().T
    dip = frame.T @ local @ frame
    plane = analyzer @ dip[:2, :2] @ analyzer.conj().T
    # Each setting has two fixed effects; the receiver never uses the source sign.
    signs = [float(np.vdot(v, plane @ v).real) for v in axis]
    trans = float(plane.trace().real)
    np.testing.assert_allclose(sum(signs), trans, rtol=0, atol=1e-12)
    longitudinal = float(dip[2, 2].real)
    arrived = float(local.trace().real)
    return [*signs, longitudinal, arrived - trans - longitudinal, 1 - arrived]


def run():
    board = dipole.carrier.federation(1)
    owners, gates = dipole.path_and_gates(board)
    h = dipole.generator(board, "unit", dipole.candidate.KAPPA)
    energies, vec = eigh(h, driver="evd")
    source = vec[:12].T @ dipole.frame(gates[0])
    rng = np.random.default_rng(292903)
    panels = []
    with (
        np.load(dipole.HERE / "counts.npz") as pilot,
        np.load(dipole.HERE / "flight.npz") as flight,
    ):
        for step, owner in enumerate(owners[1:], 1):
            raw = pilot[f"level1_full_unit_8.0_d{step}_quadratures"]
            sharp = flight[f"n80_d{step}_counts"][:, 0]
            dest = vec[owner * 12 : owner * 12 + 12]
            frame = dipole.frame(gates[step])
            readings = []
            for i, t in enumerate(TIMES):
                q = 2 * raw[i, :, :2, :2] / dipole.SHOTS - 1
                left, _, right_h = np.linalg.svd(q[0] + 1j * q[1])
                analyzer = right_h.conj().T @ left.conj().T
                np.testing.assert_allclose(
                    analyzer.conj().T @ analyzer, np.eye(2), atol=1e-12
                )
                sharp_click = int(sharp[int(2 * t)])
                rec = dict(
                    exposure=int(t),
                    sharp_clicks=sharp_click,
                    sharp_interval=bound(sharp_click, survey.SHOTS),
                    arms={},
                )
                for arm in ARMS:
                    dt = t / 2 if arm == "half_exchange" else t
                    cell = dest @ (np.exp(-1j * dt * energies)[:, None] * source)
                    rows = []
                    for index, word in enumerate(INPUTS):
                        start = 2 * (index // 2)
                        axis = INPUTS[start : start + 2]
                        p = detector(
                            cell[:, :2], analyzer, frame, word, axis, arm == "dephased"
                        )
                        # Archive signs relative to the recorded input sign. A permutation
                        # of the multinomial categories is a relabeling of this same POVM.
                        sign = index % 2
                        rows.append(draw(rng, [p[sign], p[1 - sign], *p[2:]]))
                    rec["arms"][arm] = dict(counts=rows, **summary(rows))
                    if arm == "native":
                        axial = abs(frame.T @ cell[:, 2]) ** 2
                        arrival = float(np.vdot(cell[:, 2], cell[:, 2]).real)
                        rec["axial_counts"] = draw(
                            rng,
                            [
                                axial[:2].sum(),
                                axial[2],
                                arrival - axial.sum(),
                                1 - arrival,
                            ],
                        )
                        for word in INPUTS[:2]:
                            np.testing.assert_allclose(
                                detector(
                                    cell[:, :2],
                                    analyzer,
                                    frame,
                                    word,
                                    INPUTS[:2],
                                    False,
                                ),
                                detector(
                                    cell[:, :2], analyzer, frame, word, INPUTS[:2], True
                                ),
                            )
                readings.append(rec)
            panels.append(dict(owner=owner, steps=step, readings=readings))
    return dict(
        route=owners, gates=gates, shots=SHOTS, sharp_shots=survey.SHOTS, panels=panels
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
        if data != json.loads((HERE / "records.json").read_text()):
            raise ValueError("count bank disagrees")
    final = data["panels"][-1]["readings"][-1]["arms"]
    if final["native"]["score_interval"][0] <= 0 or any(
        r["arms"]["dephased"]["score_interval"][0] > 0
        for p in data["panels"]
        for r in p["readings"]
    ):
        raise ValueError(
            "third receiver or dephasing control failed the coherence gate"
        )
    for p in data["panels"]:
        print(
            p["steps"],
            [
                (r["exposure"], r["arms"]["native"]["score_interval"])
                for r in p["readings"]
            ],
        )


if __name__ == "__main__":
    main()
