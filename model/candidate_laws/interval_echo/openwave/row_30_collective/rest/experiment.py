"""Detector experiment on the longest grown rings of two already exposed boards."""

import argparse
import json
from pathlib import Path

import numpy as np

from ..experiment import K, wedge_norm
from ..formation import native
from .commission import prepare_ring
from .physics import full_spin_check, propagation
from .read import PAIRS, decode, profile

HERE = Path(__file__).resolve().parent


def acquire(n, seed):
    times = 0.08 * np.arange(4096)
    shots = 2**22
    rng = np.random.default_rng([seed, 882001])
    hand = 2**20

    def counts(mean, reps):
        assert np.max(abs(mean)) < 1 + 1e-10
        return rng.binomial(reps, np.clip((1 + mean) / 2, 0, 1)).astype("uint32")

    bank = dict(
        dimers=n,
        shots=shots,
        hand_shots=hand,
        cos=counts(np.cos(times), hand),
        sin=counts(np.sin(times), hand),
    )
    f = profile(n)
    sector = -((-1) ** n)
    for tag, coupling in (("native", K), ("changed", K / 2)):
        u, s = propagation(n, coupling, sector)
        co = np.cos(times[:, None] * s)
        evolved = [(co * (u.T @ wave)) @ u.T for wave in f.T]
        for a, b in PAIRS:
            fa, fb = f[:, a], f[:, b]
            aa, ab, bb = evolved[a] @ fa, evolved[a] @ fb, evolved[b] @ fb
            scale = wedge_norm(fa, fb)
            key = f"{tag}_{a}{b}_"
            bank[key + "pair"] = counts((aa * bb - ab**2) / scale**2, shots)
            if tag == "native" and b <= 2:
                # Each counter estimates a product of two addressed Pauli
                # strings; setting probabilities implement its wedge weights.
                for i, (ref, y) in enumerate(
                    (
                        (fa, aa[:, None] * evolved[b] - ab[:, None] * evolved[a]),
                        (fb, ab[:, None] * evolved[b] - bb[:, None] * evolved[a]),
                    )
                ):
                    normalization = scale * (abs(ref).sum() - abs(ref))
                    bank[key + "scan" + str(i)] = counts(y / normalization, shots)
    return bank


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record", action="store_true")
    record = cli.parse_args().record
    results = dict(native_input_residual=native(), spin_checks={}, rings=[])
    for n in (5, 6, 7):
        for flip in (1, -1):
            error = full_spin_check(n, flip, K)
            results["spin_checks"][f"{n},{flip}"] = float(error)
    for seed in (312142, 927016):
        ring, commissioning = prepare_ring(1, seed)
        bank = acquire(len(ring), seed)
        path = HERE / f"counts_{seed}.npz"
        if record:
            with path.open("xb") as out:
                np.savez_compressed(out, **bank)
        with np.load(path, allow_pickle=False) as stored:
            assert set(stored) == set(bank)
            for key in bank:
                np.testing.assert_array_equal(stored[key], bank[key])
            readings = {tag: decode(stored, tag) for tag in ("native", "changed")}
            halves = [
                decode(stored, "native", start, stop)["rest"]
                for start, stop in ((0, 2048), (2048, 4096))
            ]
        results["rings"].append(
            dict(
                commission=commissioning,
                flip_sector=-((-1) ** len(ring)),
                readings=readings,
                independent_half_rests=halves,
            )
        )
    path = HERE / "results.json"
    if record:
        with path.open("x") as out:
            json.dump(results, out, indent=2)
            out.write("\n")
    else:
        reference = json.loads(path.read_text())
        np.testing.assert_allclose(
            results["native_input_residual"],
            reference["native_input_residual"],
            atol=1e-12,
            rtol=0,
        )
        assert results["spin_checks"].keys() == reference["spin_checks"].keys()
        for key, actual in results["spin_checks"].items():
            np.testing.assert_allclose(
                actual, reference["spin_checks"][key], atol=1e-11, rtol=0
            )
        for measured, saved in zip(results["rings"], reference["rings"], strict=True):
            assert (
                measured["commission"] == saved["commission"]
                and measured["flip_sector"] == saved["flip_sector"]
            )
            for tag in ["native", "changed"]:
                a, b = measured["readings"][tag], saved["readings"][tag]
                for key in (
                    "omega",
                    "held_sum_error",
                    "max_difference_line_error",
                    "clock_step",
                    "rest",
                ):
                    np.testing.assert_allclose(a[key], b[key], rtol=0, atol=1e-8)
                for key, values in a["pair_lines"].items():
                    np.testing.assert_allclose(
                        values, b["pair_lines"][key], rtol=0, atol=1e-8
                    )
                for m, s in zip(
                    a.get("spatial", []), b.get("spatial", []), strict=True
                ):
                    np.testing.assert_allclose(
                        list(m.values()), list(s.values()), rtol=0, atol=1e-8
                    )
            np.testing.assert_allclose(
                measured["independent_half_rests"],
                saved["independent_half_rests"],
                rtol=0,
                atol=1e-8,
            )
    for r in results["rings"]:
        print(
            r["commission"]["seed"],
            r["readings"]["native"]["rest"],
            r["readings"]["changed"]["rest"],
            r["readings"]["native"]["spatial"][0],
            r["independent_half_rests"],
        )


if __name__ == "__main__":
    main()
