"""Fit local triplet frames and measure their held response on two bound sources."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo import candidate, currents
from model.candidate_laws.interval_echo.openwave.quark_winding.bound_face import (
    measure as bound,
)
from model.executable_model import carrier

from .layouts import PERMS, SIGNS
from .measure import HERE, SHOTS, native_pack, packet, spectrum, validate

EPS = np.zeros((3, 3, 3))
for p, v in zip(PERMS, SIGNS, strict=True):
    EPS[p] = v / np.sqrt(6)


def optimize(tensor, rng):
    best = -1
    runs = []
    for _ in range(16):
        frames = [np.linalg.qr(rng.normal(size=(n, 3)))[0] for n in tensor.shape]
        last = -1
        for iteration in range(240):
            for i in range(3):
                # Cyclic permutation of the three tensor indices preserves epsilon.
                current = tensor.transpose(i, (i + 1) % 3, (i + 2) % 3)
                mat = np.einsum(
                    "ijk,jb,kc,abc->ia",
                    current,
                    frames[(i + 1) % 3],
                    frames[(i + 2) % 3],
                    EPS,
                    optimize=True,
                )
                u, _, v = np.linalg.svd(mat, full_matrices=False)
                frames[i] = u @ v
            overlap = float(
                np.einsum("ijk,ia,jb,kc,abc", tensor, *frames, EPS, optimize=True)
            )
            if abs(overlap - last) < 1e-11:
                break
            last = overlap
        runs.append([iteration + 1, overlap**2])
        if overlap**2 > best:
            winner = [x.copy() for x in frames]
            best = overlap**2
    return winner, runs


def face_then_colour(source, idx, face, frames, admission, rng):
    selected = source * face
    f = float(np.linalg.norm(selected) ** 2)
    code = np.einsum("ijk,ia,jb,kc->abc", selected[idx], *frames, optimize=True)
    trip = float(np.linalg.norm(code) ** 2)
    det = float(abs(np.vdot(EPS, code)) ** 2)
    pulled = np.einsum("ia,jb,kc,abc->ijk", *frames, EPS, optimize=True)
    dephased = float(np.sum(abs(pulled * selected[idx]) ** 2))
    return dict(
        face=f,
        face_and_triplets=trip,
        face_and_singlet=det,
        conditional_singlet=det / trip,
        dephased_joint_singlet=dephased,
        counts=packet(
            rng,
            [
                admission * det,
                admission * (trip - det),
                admission * (f - trip),
                admission * (1 - f),
                1 - admission,
            ],
        ),
    )


def finite_filter(eigenvalues, eigenvectors, energy, launch, idx, face, frames, rng):
    # The central Fourier output has the time-average amplitude of Q delayed inputs.
    q = 2**20
    step = 1 / 32
    x = (eigenvalues - energy) * step / (2 * np.pi)
    filt = np.exp(-1j * np.pi * (q - 1) * x) * np.sinc(q * x) / np.sinc(x)
    coeff = eigenvectors[launch] * filt
    admission = float(np.linalg.norm(coeff) ** 2)
    coeff /= np.sqrt(admission)
    outside = float(np.linalg.norm(coeff[abs(eigenvalues - energy) > 1e-8]) ** 2)
    records = {}
    for time in (0, 8, 32):
        state = eigenvectors @ (coeff * np.exp(-1j * (eigenvalues - energy) * time))
        records[str(time)] = face_then_colour(state, idx, face, frames, admission, rng)
    return dict(
        timing_arms=q,
        step=step,
        source_acceptance=admission,
        energy_outside_target=outside,
        trace_distance_from_ideal_at_all_free_times=np.sqrt(outside),
        records=records,
    )


def assess(source, admission, idx, face, word, fit_rng, rng):
    t = source[idx]
    frames, runs = optimize(t, fit_rng)
    colour = currents._service(carrier.cell(), False, 1.0)[4][0]

    def complete(three):
        tail = np.linalg.qr(three, mode="complete")[0][:, 3:]
        return np.column_stack([three, tail])

    target = complete(colour)
    compiled = []
    for frame in frames:
        local = np.zeros((12, 3))
        local[: len(frame)] = frame
        change = target @ complete(local).T
        np.testing.assert_allclose(change @ local, colour, atol=1e-12)
        compiled.append(native_pack(change, candidate.KAPPA))
    code = np.einsum("ijk,ia,jb,kc->abc", t, *frames, optimize=True)
    p = float(np.vdot(code, code).real)
    singlet = float(abs(np.vdot(EPS, code)) ** 2)
    words = {}
    for name, u in word.items():
        for together in (True, False):
            ops = [u] * 3 if together else [u, np.eye(3), np.eye(3)]
            turned = np.einsum("ia,jb,kc,abc->ijk", *ops, code, optimize=True)
            dark = np.linalg.norm(turned - code) ** 2 / 4
            bright = np.linalg.norm(turned + code) ** 2 / 4
            # The finite compiled word has very small leakage into other triplets;
            # count it when it appears, along with both source/decoder failures.
            loss = max(0.0, p - bright - dark)
            words[name + ("_all" if together else "_one")] = dict(
                dark_given_three=float(dark / p),
                counts=packet(
                    rng,
                    [
                        admission * bright,
                        admission * dark,
                        admission * loss,
                        admission * (1 - p),
                        1 - admission,
                    ],
                ),
            )
    upper = []
    for k in range(3):
        mat = np.moveaxis(t, k, 0).reshape(t.shape[k], -1)
        s = np.linalg.svd(mat, compute_uv=False)
        upper.append(float(sum(s[:3]) ** 2 / 3))
    return dict(
        source_acceptance=admission,
        face_probability=float(np.linalg.norm(source[face]) ** 2),
        one_per_domain=float(np.linalg.norm(t) ** 2),
        three_triplets=p,
        singlet=singlet,
        singlet_given_three=singlet / p,
        frames=[u.tolist() for u in frames],
        optimization_history=runs,
        compiled_receiver_frames=compiled,
        schmidt_squared_overlap_upper_bounds=upper,
        face_then_colour=face_then_colour(source, idx, face, frames, admission, rng),
        singlet_counts=packet(
            rng,
            [
                admission * singlet,
                admission * (p - singlet),
                admission * (1 - p),
                1 - admission,
            ],
        ),
        held_turns=words,
    )


def acquisition_stream(position):
    # Positions from the deposited seed-133922 acquisition. Separate generators
    # keep decoder selection independent of modifications to detector sampling.
    rng = np.random.default_rng(133922)
    state = rng.bit_generator.state
    state["state"]["state"] = position
    rng.bit_generator.state = state
    return rng


def run():
    data = json.loads((HERE / "results.json").read_text())["strengths"]["8"][
        "all_bound_receiver_domains"
    ]
    pick = np.argmax(data["one_per_group_maxima"])
    groups = data["maximizing_groups"][pick]
    c = carrier.cell()
    sec = BinarySector(c, 3)
    energies, vectors, lines = spectrum(sec.hamiltonian(candidate.KAPPA, False))
    en, v = min(lines, key=lambda x: abs(x[0] - data["energies"][pick]))
    idx = np.array([sec.index[tuple(sorted(t))] for t in product(*groups)]).reshape(
        tuple(map(len, groups))
    )
    q = v[idx.ravel()].T @ v[idx.ravel()]
    vals, R = np.linalg.eigh(q)
    max_source = v @ R[:, -1]
    # At a repeated top eigenvalue choose the projection of the lowest-index sharp
    # face with nonzero weight in the maximizing subspace. This is basis invariant.
    space = v @ R[:, vals >= vals[-1] - 1e-9]
    face = sec.agreed == 3
    tag = next(i for i in np.flatnonzero(face) if np.linalg.norm(space[i]) > 1e-8)
    max_source = space @ (space[tag] / np.linalg.norm(space[tag]))
    face_indices = np.flatnonzero(face)
    probs = (v[face_indices] ** 2).sum(axis=1)
    f = face_indices[np.flatnonzero(probs >= max(probs) - 1e-10)[0]]
    acc = float(v[f] @ v[f])
    sharp = v @ v[f] / np.sqrt(acc)
    _, _, w, _, _, _ = bound.colour_words()
    bands = currents._service(c, False, 1.0)[4][0]
    words = {s: bands.T @ w[s] @ bands for s in ("sham", "centre", "half_turn")}
    ms = assess(
        max_source,
        1.0,
        idx,
        face,
        words,
        acquisition_stream(34797434909415434004239637893737029124),
        acquisition_stream(202255707923641132210905905550675206942),
    )
    ss = assess(
        sharp,
        acc,
        idx,
        face,
        words,
        acquisition_stream(65286932915696487277475557556206992540),
        acquisition_stream(176947738671518749063623659511712338092),
    )
    rng = acquisition_stream(193853122577654522354345806977955088235)
    finite = finite_filter(
        energies, vectors, en, f, idx, face, [np.array(g) for g in ss["frames"]], rng
    )
    return currents._portable(
        dict(
            finite_filter=finite,
            acquisitions=SHOTS,
            energy=en,
            groups=groups,
            face_for_energy_filter=sec.states[f],
            top_occupancy_multiplicity=space.shape[1],
            max_one_per_group=ms,
            sharp_energy_filtered=ss,
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = run()
    if args.output:
        with args.output.open("x") as handle:
            json.dump(data, handle, indent=2, allow_nan=False)
            handle.write("\n")
    else:
        validate(data, json.loads((HERE / "frames.json").read_text()))
        print("Reproduced the fitted colour frames and held singlet tests.")


if __name__ == "__main__":
    main()
