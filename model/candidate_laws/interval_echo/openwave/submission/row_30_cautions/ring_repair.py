"""Reconstruct a ring flip-sector correlation from two local Pauli measurements.

Each arm starts from independently randomized one-load dimers (I / 2 per dimer).
Write O for the first physical Pauli string, B for the terminal one and F for
the product of all dimer X's.  The arms measure O or f O F locally and then B.
On I/D their product-sign means sum to tr[B(t) O (1+fF)]/D.  The desired
sector exists in this ensemble identity. Individual measurements do disturb it.
"""

import argparse
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from scipy.linalg import eigh

from ...row_30_collective.experiment import K, verify, wedge_norm
from ...row_30_collective.formation import native
from ...row_30_collective.rest.commission import prepare_ring
from ...row_30_collective.rest.physics import propagation
from ...row_30_collective.rest.read import PAIRS, decode, profile

EXPOSED = (312142, 927016)
SHOTS = 2**21  # per arm per detector setting
HAND = 2**20
TIMES = 0.08 * np.arange(4096)
PAULI = {
    "I": np.eye(2),
    "X": np.array([[0, 1], [1, 0]]),
    "Y": np.array([[0, -1j], [1j, 0]]),
    "Z": np.diag([1, -1]),
}


def axes(n, left, right, signed=False):
    # i X_[left,right) product(-Z_(j-1) Z_j) in increasing physical site order.
    # Multiplication by F exchanges the two end axes and complements the X's.
    assert 0 <= left < right < n
    factors = ["I"] * n
    factors[(left - 1) % n], factors[right - 1] = "Z", "Y"
    factors[left : right - 1] = ["X"] * (right - left - 1)
    if signed:
        exchange = {"I": "X", "X": "I", "Y": "Z", "Z": "Y"}
        factors = [exchange[a] for a in factors]
    return factors, (-1) ** (right - left)


def local_projectors(factors, overall):
    # Enumerate local output words and their actual tensor-product projectors.
    # Site zero is the least significant bit of a physical orientation word.
    entries = [(overall, np.ones((1, 1)))]
    for a in factors:
        options = (
            [(1, PAULI["I"])]
            if a == "I"
            else [(s, (PAULI["I"] + s * PAULI[a]) / 2) for s in (-1, 1)]
        )
        entries = [(s * t, np.kron(p, q)) for t, p in options for s, q in entries]
    return entries


def physical_check(n, coupling):
    size = 2**n
    word = np.arange(size)
    z = 1 - 2 * ((word[:, None] >> np.arange(n)) & 1)
    h = np.diag(coupling * np.sum(z * np.roll(z, 1, axis=1), axis=1))
    for j in range(n):
        h[word, word ^ (1 << j)] = -0.5
    energies, eigen = eigh(h)
    identity = np.eye(size)
    flip = identity[::-1]
    err, losses, plain_biases = [], [], []
    times = (0.37, 2.15, 5.81)
    twists = {}
    for f in (1, -1):
        u, omega = propagation(n, coupling, f)
        twists[f] = [(u * np.cos(t * omega)) @ u.T for t in times]
    old_sector = (identity + flip) / size
    for left, right, dest, end in (
        (0, 1, 2, 3),
        (0, 2, 1, n - 1),
        (1, n - 1, 1, n - 1),
    ):
        plain = local_projectors(*axes(n, left, right))
        signed = local_projectors(*axes(n, left, right, signed=True))
        final = local_projectors(*axes(n, dest, end))
        # Weighted postmeasurement states sum s*P_s (I/D) P_s for each full word.
        src = [sum(s * (p @ p) for s, p in arm) / size for arm in (plain, signed)]
        terminal = sum(s * p for s, p in final)
        np.testing.assert_allclose(src[1], src[0] @ flip, atol=1e-14)
        np.testing.assert_allclose(sum(p for _, p in plain), identity, atol=1e-14)
        disturbed = sum(p @ old_sector @ p for _, p in plain)
        losses.append(float(1 - np.trace((identity + flip) @ disturbed).real / 2))
        for i, t in enumerate(times):
            evolution = (eigen * np.exp(-1j * t * energies)) @ eigen.T
            observed = [
                np.trace(terminal @ evolution @ s @ evolution.conj().T).real
                for s in src
            ]
            # Direct physical sector target and two actual unconditioned arms.
            for f in (-1, 1):
                direct = np.trace(
                    terminal
                    @ evolution
                    @ src[0]
                    @ (identity + f * flip)
                    @ evolution.conj().T
                ).real
                target, other = [
                    np.linalg.det(twists[sign][i][np.ix_([left, right], [dest, end])])
                    for sign in (f, -f)
                ]
                plain_biases.append(abs(observed[0] - direct))
                err.extend(
                    abs(x)
                    for x in (
                        direct - target,
                        observed[0] - (target + other) / 2,
                        f * observed[1] - (target - other) / 2,
                        observed[0] + f * observed[1] - direct,
                    )
                )
    assert max(err) < 1e-11
    np.testing.assert_allclose(losses, 0.5, atol=1e-14)
    assert max(plain_biases) > 0.001
    return {
        "max_trace_error": float(max(err)),
        "single_sector_loss": max(losses),
        "max_plain_arm_alone_error": float(max(plain_biases)),
    }


def responses(n, coupling, flip):
    # The antisymmetrized bilinears sample physical strings with probability
    # proportional to absolute coefficient, folding its sign into the count.
    u, spectrum = propagation(n, coupling, flip)
    waves = profile(n)
    cosines = np.cos(TIMES[:, None] * spectrum)
    evolved = [(cosines * (u.T @ wave)) @ u.T for wave in waves.T]
    for a, b in PAIRS:
        fa, fb = waves[:, a], waves[:, b]
        aa, ab, bb = evolved[a] @ fa, evolved[a] @ fb, evolved[b] @ fb
        scale = wedge_norm(fa, fb)
        yield f"{a}{b}_pair", (aa * bb - ab**2) / scale**2
        if b <= 2:
            for j, (ref, vector) in enumerate(
                (
                    (fa, aa[:, None] * evolved[b] - ab[:, None] * evolved[a]),
                    (fb, ab[:, None] * evolved[b] - bb[:, None] * evolved[a]),
                )
            ):
                yield f"{a}{b}_scan{j}", vector / (scale * (abs(ref).sum() - abs(ref)))


def acquire(n, seed):
    streams = np.random.SeedSequence([seed, 13608437]).spawn(5)
    hand, plain, signed, plain_scan, signed_scan = (
        np.random.default_rng(child) for child in streams
    )

    def counts(rng, mean, budget):
        assert np.max(abs(mean)) < 1 + 1e-10
        return rng.binomial(budget, np.clip((1 + mean) / 2, 0, 1)).astype("<u4")

    bank = dict(
        dimers=n,
        arm_shots=SHOTS,
        hand_shots=HAND,
        cos=counts(hand, np.cos(TIMES), HAND),
        sin=counts(hand, np.sin(TIMES), HAND),
    )
    f = -((-1) ** n)
    for tag, coupling in (("native", K), ("changed", K / 2)):
        target = responses(n, coupling, f)
        other = responses(n, coupling, -f)
        for (name, ct), (other_name, co) in zip(target, other, strict=True):
            assert name == other_name
            key = f"{tag}_{name}"
            plain_rng, signed_rng = (
                (plain_scan, signed_scan)
                if tag == "changed" and "_scan" in name
                else (plain, signed)
            )
            bank["plain_" + key] = counts(plain_rng, (ct + co) / 2, SHOTS)
            # The factor f is part of the measured signed string.
            bank["signed_" + key] = counts(signed_rng, (ct - co) / 2, SHOTS)
    return bank


def read_bank(path):
    # The stored reader estimates frequencies and normalized spatial planes.
    # With equal budgets the pooled sign mean is exactly half the sector
    # estimate. A positive constant gain preserves both reader outputs.
    with np.load(path, allow_pickle=False) as bank:
        n, reps = int(bank["dimers"]), int(bank["arm_shots"])
        pooled = dict(
            dimers=n,
            shots=2 * reps,
            hand_shots=int(bank["hand_shots"]),
            cos=bank["cos"],
            sin=bank["sin"],
        )
        digests = {}
        detector_shots = 0
        sigma_max = 0.0
        for key in sorted(bank.files):
            counts = bank[key]
            digests[key] = hashlib.sha256(
                np.asarray(counts, dtype="<u8").tobytes()
            ).hexdigest()
            if not key.startswith("plain_"):
                continue
            setting = key.removeprefix("plain_")
            signed = bank["signed_" + setting]
            assert np.all(counts <= reps) and np.all(signed <= reps)
            pooled[setting] = counts.astype(np.int64) + signed
            means = [2 * arm.astype(float) / reps - 1 for arm in (counts, signed)]
            reconstructed = means[0] + means[1]
            np.testing.assert_array_equal(
                reconstructed, 2 * (pooled[setting] / reps - 1)
            )
            # Unbiased plug-in variance of the sum of two independent sign means.
            variance = (2 - means[0] ** 2 - means[1] ** 2) / (reps - 1)
            sigma_max = max(sigma_max, float(np.sqrt(variance.max())))
            detector_shots += 2 * reps * counts.size

        native_read = decode(pooled, "native")
        changed_read = decode(pooled, "changed")
        halves = [
            decode(pooled, "native", start, start + 2048)["rest"] for start in (0, 2048)
        ]
        return {
            "readings": {"native": native_read, "changed": changed_read},
            "independent_half_rests": halves,
            "count_sha256_uint64_le": digests,
            "detector_trials": detector_shots,
            "hand_trials": 2 * int(bank["hand_shots"]) * bank["cos"].size,
            "max_sampled_sector_mean_standard_error": sigma_max,
        }


def run(directory):
    direct = {
        f"{n},{coupling:g}": physical_check(n, coupling)
        for n in (5, 6)
        for coupling in (K, K / 2)
    }
    rings = []
    for seed in EXPOSED:
        route, commissioning = prepare_ring(1, seed)
        n = len(route)
        path = directory / f"two_arm_{seed}.npz"
        with path.open("xb") as stream:
            np.savez_compressed(stream, **acquire(n, seed))
        reading = read_bank(path)
        rings.append(
            dict(
                dimers=n,
                cells=commissioning["ring_cells"],
                seed=seed,
                circuit_sha256=commissioning["circuit_sha256"],
                sector=-((-1) ** n),
                **reading,
            )
        )
    return dict(
        scope="Exploratory exposed closed rings, coherent unread connectors, mixed dimers, separate local measurements",
        fit="Stored frequency and spatial reader on equally pooled raw counts",
        per_arm_shots=SHOTS,
        sector_to_pooled_mean_gain=2,
        worst_case_sector_mean_standard_error=float(np.sqrt(2 / SHOTS)),
        native_input_residual=native(),
        direct_projector_checks=direct,
        rings=rings,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bank-dir", type=Path, help="Keep the complete generated integer banks"
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Write an exploratory receipt instead of checking it",
    )
    parser.add_argument(
        "--read-bank", type=Path, help="Read a saved counts-only bank and exit"
    )
    args = parser.parse_args()
    if args.read_bank:
        print(json.dumps(read_bank(args.read_bank), indent=2, allow_nan=False))
        raise SystemExit(0)
    with TemporaryDirectory() as scratch:
        directory = args.bank_dir or Path(scratch)
        directory.mkdir(parents=True, exist_ok=True)
        result = run(directory)
    if args.output:
        with args.output.open("x") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        verify(result, json.loads(Path(__file__).with_suffix(".json").read_text()))
        print("Reproduced both exposed ring measurements.")
