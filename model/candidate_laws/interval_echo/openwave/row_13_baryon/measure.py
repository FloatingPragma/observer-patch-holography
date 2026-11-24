"""Release and assess three independently readable orbital triplets."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier, rotations
from model.candidate_laws.interval_echo.openwave.quark_winding.bound_face import (
    measure as bound,
)

from .layouts import SIGNS, enumerate_layouts, singlet_support
from .partitions import maxima, transfer_certificate

HERE = Path(__file__).resolve().parent
SHOTS = 2**19
SEED = 133921


def detector(sector, groups):
    basis = np.zeros(len(sector.states))
    support = singlet_support(groups)
    basis[[sector.index[f] for f in support]] = SIGNS / np.sqrt(6)
    allowed = np.zeros(len(basis), bool)
    for addresses in product(*groups):
        allowed[sector.index[tuple(sorted(addresses))]] = True
    return basis, allowed


def packet(rng, probs):
    p = np.asarray(probs, dtype=float)
    if min(p) < -1e-10 or abs(p.sum() - 1) > 1e-10:
        raise ValueError("Invalid outcome probabilities")
    p = np.maximum(0, p).round(12)
    return rng.multinomial(SHOTS, p / p.sum()).tolist()


def read(rng, waves, s, cycle, allowed, face, admitted=1):
    # Each column is weighted by the square root of its mixture probability.
    weight = abs(waves) ** 2
    inside = weight[allowed].sum()
    f = weight[face].sum()

    def square(proj, w):
        return float(sum(abs(proj @ w) ** 2))

    good = square(s, waves)
    joint = square(s[face], waves[face])
    shifted = square(cycle, waves)
    return dict(
        face=f,
        code=inside,
        singlet=good,
        singlet_after_face=joint,
        one_receiver_colour_cycle=shifted,
        face_counts=packet(rng, [admitted * f, admitted * (1 - f), 1 - admitted]),
        singlet_counts=packet(
            rng,
            [
                admitted * good,
                admitted * (inside - good),
                admitted * (1 - inside),
                1 - admitted,
            ],
        ),
        after_face_counts=packet(
            rng,
            [
                admitted * joint,
                admitted * (f - joint),
                admitted * (1 - f),
                1 - admitted,
            ],
        ),
        cycle_counts=packet(
            rng,
            [
                admitted * shifted,
                admitted * (inside - shifted),
                admitted * (1 - inside),
                1 - admitted,
            ],
        ),
    )


def spectrum(h):
    values, vecs = np.linalg.eigh(h)
    cuts = np.flatnonzero(np.diff(values) > 1e-8) + 1
    return (
        values,
        vecs,
        [
            (float(values[idx].mean()), vecs[:, idx])
            for idx in np.split(np.arange(len(h)), cuts)
        ],
    )


def panel(sector, groups, k, rng, eig, fragment_min):
    evals, vecs, lines = eig
    face = sector.agreed == 3
    s, allowed = detector(sector, groups)
    shift = ((groups[0][1], groups[0][2], groups[0][0]), groups[1], groups[2])
    cycle, _ = detector(sector, shift)
    mask = face & (abs(s) > 0)
    birth = np.flatnonzero(mask)[0]
    sharp = np.eye(len(s))[:, birth]
    cut = s * face
    cut /= np.linalg.norm(cut)
    mixture = np.eye(len(s))[:, abs(s) > 0] / np.sqrt(6)
    prepared = {
        "sharp_face": sharp[:, None],
        "singlet": s[:, None],
        "face_filtered_singlet": cut[:, None],
        "unsigned_permutations": abs(s)[:, None],
        "incoherent_six": mixture,
    }
    flights = {}
    for name, source in prepared.items():
        admitted = (
            np.linalg.norm(s * face) ** 2 if name == "face_filtered_singlet" else 1.0
        )
        waves = vecs.T @ source
        flights[name] = {
            str(wait): read(
                rng,
                vecs @ (np.exp(-1j * evals * wait)[:, None] * waves),
                s,
                cycle,
                allowed,
                face,
                admitted,
            )
            for wait in (0, 1, 8, 32)
        }
    bound_lines = [(en, b) for en, b in lines if en < fragment_min - 1e-8]
    lookup = {"all": lines, "bound": bound_lines}
    rows = {}
    for name, pool in lookup.items():
        if not pool:
            continue
        en, v = max(pool, key=lambda line: np.linalg.norm(line[1].T @ s) ** 2)
        phi = v @ (v.T @ s)
        overlap = float(phi @ phi)
        phi /= np.sqrt(overlap)
        projection = v @ (v.T @ sharp)
        accepted = float(projection @ projection)
        meter = (
            read(
                rng,
                projection[:, None] / np.sqrt(accepted),
                s,
                cycle,
                allowed,
                face,
                accepted,
            )
            if accepted > 1e-13
            else None
        )
        rows[name] = dict(
            energy=en,
            line_dimension=v.shape[1],
            maximal_singlet=overlap,
            face_on_maximizer=float(np.sum(phi[face] ** 2)),
            sharp_filter_acceptance=accepted,
            filtered_sharp=meter,
        )
    h = sector.hamiltonian(k, False)
    spare = sorted(set(range(12)) - set(sum((tuple(g) for g in groups), ())))
    leak = np.array([bool(set(t).intersection(spare)) for t in sector.states])
    return dict(
        groups=groups,
        birth_face=sector.states[birth],
        spare_sites=spare,
        exact_spare_hopping_norm_squared=float(np.linalg.norm((h @ s)[leak]) ** 2),
        spectral_optima=rows,
        free_flights=flights,
    )


def native_pack(change, k):
    # Isolated binary contacts have a common endpoint phase and an i-swap.
    # Choose port phases on their two sides to make a real Givens rotation.
    reduced = change.copy()
    moves = []
    for col in range(len(reduced)):
        for row in range(len(reduced) - 1, col, -1):
            a, b = reduced[row - 1 : row + 1, col]
            if abs(b) < 1e-14:
                continue
            angle = np.arctan2(b, a)
            c, s = np.cos(angle), np.sin(angle)
            reduced[row - 1 : row + 1] = (
                np.array([[c, s], [-s, c]]) @ reduced[row - 1 : row + 1]
            )
            moves.append((row - 1, row, angle))
    rebuilt = np.diag(np.sign(np.diag(reduced))).astype(complex)
    for a, b, angle in reversed(moves):
        wait = (2 * angle) % (4 * np.pi)
        block = np.exp(-1j * (k + 0.5) * wait) * (
            np.cos(wait / 2) * np.eye(2)
            + 1j * np.sin(wait / 2) * np.array([[0, 1], [1, 0]])
        )
        block = (
            np.exp(1j * (k + 0.5) * wait) * np.diag([1, -1j]) @ block @ np.diag([1, 1j])
        )
        rebuilt[[a, b]] = block @ rebuilt[[a, b]]
    err = float(np.linalg.norm(rebuilt - change))
    np.testing.assert_allclose(rebuilt, change, atol=1e-12)
    return dict(givens_contacts=len(moves), operator_residual=err)


def meter_calibration(groups):
    cell = carrier.cell()
    _, _, _, _, bands, _, _, _, _, _ = currents._service(cell, False, 1.0)
    colour = bands[0]
    phase, anti, words, _, _, identical = bound.colour_words()
    full = np.linalg.qr(colour, mode="complete")[0]
    # Replace the QR sign choices of the first three columns by the calibrated axes.
    frame = np.column_stack([colour, full[:, 3:]])
    errors = []
    compiled = []
    for group in groups:
        order = [*group, *sorted(set(range(12)) - set(group))]
        perm = np.eye(12)[:, order]
        change = frame @ perm.T
        errors.append(float(np.linalg.norm(change[:, group] - colour)))
        compiled.append(native_pack(change, candidate.KAPPA))
        np.testing.assert_allclose(change.T @ change, np.eye(12), atol=1e-13)
    return {
        "pack_basis_residuals": errors,
        "isolated_native_conversion": compiled,
        "three_fundamental_singlet": bound.colour_read(anti, phase, anti, words),
        "three_identical_colours": bound.colour_read(identical, phase, anti, words),
    }


def run():
    cell = carrier.cell()
    tri = BinarySector(cell, 3)
    rng = np.random.default_rng(SEED)
    histogram, groups = enumerate_layouts(cell)
    rot = next(g for g in sorted(rotations()) if tuple(g) != tuple(range(12)))
    other = tuple(tuple(rot[x] for x in row) for row in groups)
    result = {
        "seed": SEED,
        "enrolled_per_panel": SHOTS,
        "geometry_histogram": histogram,
        "calibration": meter_calibration(groups),
        "strengths": {},
    }
    for k in (candidate.KAPPA, 0):
        h = tri.hamiltonian(k, False)
        es = spectrum(h)
        fragment = sum(
            np.linalg.eigvalsh(BinarySector(cell, n).hamiltonian(k, False))[0]
            for n in (1, 2)
        )
        data = dict(
            fragment_ground=float(fragment),
            trimer_ground=float(es[0][0]),
            original=panel(tri, groups, k, rng, es, fragment),
            rotated=panel(tri, other, k, rng, es, fragment),
        )
        if k:
            data["all_bound_receiver_domains"] = maxima(
                cell, tri, [(en, v) for en, v in es[2] if en < fragment - 1e-8]
            )
        data["modular_observability"] = transfer_certificate(cell, tri, h)
        result["strengths"][str(int(k))] = data
    return currents._portable(result)


def validate(fresh, stored):
    if isinstance(stored, dict):
        assert set(fresh) == set(stored)
        for k, v in stored.items():
            validate(fresh[k], v)
    elif isinstance(stored, list):
        assert len(fresh) == len(stored)
        for x, y in zip(fresh, stored, strict=True):
            validate(x, y)
    elif isinstance(stored, (int, float)):
        np.testing.assert_allclose(fresh, stored, rtol=0, atol=1e-8)
    else:
        assert fresh == stored


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
        validate(data, json.loads((HERE / "results.json").read_text()))
        print("Reproduced three-triplet colour, binding and partition records.")


if __name__ == "__main__":
    main()
