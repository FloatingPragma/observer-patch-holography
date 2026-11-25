"""Sharp and reciprocal phase-sensitive launch and detection on one fixed source."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from ...material_clocks.measure import blank_clock
from .. import measure as optics

HERE = Path(__file__).resolve().parent
PHASES = {
    "plus_plus": (1, 1),
    "minus_minus": (-1, -1),
    "plus_minus": (1, -1),
    "minus_plus": (-1, 1),
}
T = 4
DWELL = 2 * np.arctan2(3, 4)


def splitter(charge, sign=1):
    g = np.eye(5, dtype=complex)
    # A positive blank-seam activation of duration 2 pi - DWELL realizes
    # the conjugate splitter up to the relative phase on the two-port block.
    g[:2, :2] = np.array([[4, sign * 3j], [sign * 3j, 4]]) / 5
    dimer = blank_clock(charge, 8)
    # With sharp mutually exclusive inputs and outputs only the phase
    # difference between the two rows within the active block matters.
    gate = expm(-1j * DWELL * dimer)
    np.testing.assert_allclose(
        gate / gate[0, 0], np.array([[1, 0.75j], [0.75j, 1]]), atol=8e-14
    )
    np.testing.assert_allclose(g.conj().T @ g, np.eye(5), atol=1e-14)
    return g


def amplitudes(board, route, rings, charge, odd, strength):
    a, b = map(optics.beam.indices, rings)
    h = (
        optics.beam.front.native.generator(board, route, charge, odd) / 2
        if strength == 8
        else optics.beam.front.optics.ranging.generator(
            board, route, charge, strength, odd
        )
    )
    assert (h - h.T).count_nonzero() == 0
    ih = -1j * h
    x = np.zeros((h.shape[0], 20))
    x[a] = np.eye(20)
    wave = expm_multiply(T * ih, x, traceA=T * ih.diagonal().sum())
    u = []
    for i in range(3):
        u.append(wave[b].reshape(5, 4, 5, 4).transpose(0, 2, 1, 3))
        if i < 2:
            wave = ih @ wave
    return np.array(u)


def probabilities(waves, source, left, right):
    # left and right are complex symmetric native blank-seam splitters.
    # Contract the source only after forming each tag-resolved amplitude.
    p = []
    for wave in waves:
        w = np.einsum("ay,pyxij,xb->pabij", right, wave, left)

        def adj(m):
            return m.conj().swapaxes(-1, -2)

        v, d, e = w
        first = adj(v) @ d
        second = adj(v) @ e + adj(d) @ d
        p.append([adj(v) @ v, first + adj(first), second + adj(second)])
    return np.array(
        [optics.beam.probabilities(op, source) for op in np.array(p).swapaxes(0, 1)]
    )


def read(p):
    basis = np.c_[np.ones(5), optics.beam.BEARINGS]
    j = optics.matrix_jets(basis.T @ p @ basis)
    det, b, c = optics.optical_matrices(j)
    return dict(area=det, expansion=b.trace(), focusing=-c.trace())


def run():
    f = optics.beam.front.optics.pairing.formation
    board = f.federation(1)
    route = f.region(board)
    rings = optics.beam.screens(board, route)
    state = optics.beam.aperture.survey.preparations()["equilibrium"]
    assert list(route) == [0, 1, 38, 39]
    blanks = {tuple(sorted(ab)) for ab in board.seams} - {
        tuple(sorted(ab)) for ab in (route[:2], route[2:])
    }
    assert all(tuple(sorted(r[:2])) in blanks for r in rings)
    data = {}
    for strength in [8, 0]:
        rows = {}
        for q in [-1, 1]:
            waves = np.array(
                [
                    amplitudes(board, route, rings, q, odd, strength)
                    for odd in [False, True]
                ]
            )
            rev = waves.transpose(0, 1, 3, 2, 5, 4)
            g = splitter(q)
            labels = {}
            rng = np.random.default_rng(731927 + strength + q)
            shots = 2**38
            for side, w in [("forward", waves), ("reverse", rev)]:
                sharp = probabilities(w, state, np.eye(5), np.eye(5))
                both = probabilities(w, state, g, g)
                conjugate = probabilities(w, state, g.conj(), g.conj())
                avg = (
                    sum(
                        probabilities(w, state, a, b)
                        for a in [g, g.conj()]
                        for b in [g, g.conj()]
                    )
                    / 4
                )
                blur = abs(g) ** 2
                np.testing.assert_allclose(avg, blur @ sharp @ blur, atol=1e-13, rtol=0)
                inverse = np.linalg.inv(blur)
                restored = inverse @ avg @ inverse
                np.testing.assert_allclose(restored, sharp, atol=1e-13, rtol=0)
                # This one fixed sender/receiver record distinguishes a phase
                # change before taking any ratio, centroid or time derivative.
                channels = {}
                for tag, signs in PHASES.items():
                    sl, sr = signs if side == "forward" else signs[::-1]
                    p = probabilities(
                        w, state, g if sl == 1 else g.conj(), g if sr == 1 else g.conj()
                    )
                    count = int(rng.binomial(shots, p[0, 0, 0]))
                    channels[tag] = {
                        **read(p),
                        "probability00": p[0, 0, 0],
                        "click_and_complement": [count, shots - count],
                    }
                labels[side] = dict(
                    channels=channels,
                    sharp=read(sharp),
                    coherent=read(both),
                    conjugate=read(conjugate),
                    phase_averaged=read(avg),
                    corrected_phase_average=read(restored),
                    maximum_corrected_probability_jet_error=float(
                        np.max(abs(restored - sharp))
                    ),
                    pixel00=dict(
                        sharp=sharp[0, 0, 0],
                        coherent=both[0, 0, 0],
                        conjugate=conjugate[0, 0, 0],
                        dephased=avg[0, 0, 0],
                    ),
                    complete_coherent_probabilities=both[0].tolist(),
                    miss_per_launch=(1 - both[0].sum(0)).tolist(),
                )
            odd = {}
            for tag in ["sharp", "coherent", "conjugate", "phase_averaged"]:
                odd[tag] = {
                    k: labels["forward"][tag][k] - labels["reverse"][tag][k]
                    for k in ["expansion", "focusing"]
                }
            for tag in PHASES:
                odd[tag] = {
                    k: labels["forward"]["channels"][tag][k]
                    - labels["reverse"]["channels"][tag][k]
                    for k in ["expansion", "focusing"]
                }
            rows[str(q)] = dict(settings=labels, odd=odd)
        data[str(strength)] = rows
    return dict(
        native_dwell=T,
        splitter_probability=9 / 25,
        splitter_edges=[r[:2] for r in rings],
        strength=data,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if p.parse_args().write:
        (HERE / "results.json").write_text(out)
    else:
        print(out, end="")
