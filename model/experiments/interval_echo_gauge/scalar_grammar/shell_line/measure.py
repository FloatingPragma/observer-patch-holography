"""Acquire locality-resolved estimates of the quartic and complex Yukawa line."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np

from model.experiments.interval_echo_consensus.signed_effects.sequential import compare

from . import apparatus as a

HERE = Path(__file__).resolve().parent
RADIUS, STRENGTH, YUKAWA = 1.0, 1.5, 0.5
SHELL = ((0, 0), (0, 1), (1, 0), (3 / 5, 4j / 5), (0, 1 / 2), (0, 3 / 2), (1, 1j))
VAC = (0, 1)
PARALLEL = (1 / np.sqrt(2), 1j / np.sqrt(2))
LINES = (
    ("vacuum_upper", (1, 0), VAC, 1),
    ("vacuum_lower", (0, 1), VAC, 1),
    ("vacuum_imaginary", (1j, 0), VAC, 1),
    ("singlet_phase", (1, 0), VAC, 1j),
    ("singlet_zero", (1, 0), VAC, 0),
    ("exchanged", (0, 1), (1, 0), 1),
    ("parallel", PARALLEL, PARALLEL, 1),
    ("held", (4 / 5, -3j / 5), (3 / 5, 4j / 5), (3 + 4j) / 5),
)


def pair_complex(z):
    return [float(z.real), float(z.imag)]


def shell_bank(rng, cut_rng, h, correlated):
    local = a.read(h)
    # In the copy arm a sharp local observation controls preparation of
    # the remote frame in that very address before its three-rail transfer.
    remote = a.read(h, site="receiver")
    transfer = np.eye(3)[:, a.transport()["connected_columns"]]
    two = (
        local[:, None] * abs(transfer.T) ** 2 if correlated else a.joint(local, remote)
    )
    blocked = a.acquire(cut_rng, a.joint(local, a.read(h, site="cut"), a.coin(True)))
    raw = a.acquire(rng, a.joint(two, a.coin(True)))
    signs = (1, 1, -1)
    num = a.sample_mean(raw, lambda i, j, k: (1, 2, 0)[k] * signs[i] * signs[j])
    den = a.sample_mean(raw, lambda i, j, k: int(i == 2 and j == 2))
    expected = STRENGTH * (np.vdot(h, h).real - RADIUS**2) ** 2
    # The score range is [-2,2], the probability denominator range [0,1].
    interval = a.ratio_interval(num, den, 4 * a.UNIT_ERROR)
    return dict(
        h=[pair_complex(z) for z in h],
        repeated_address=correlated,
        raw=raw,
        blocked=blocked,
        numerator=num,
        reference_coincidence=den,
        measured=num / den if den else None,
        interval=interval,
        independent_prediction=float(expected),
    )


def coefficient(i, target):
    return 1 if i == target else -1 if i == 2 else 0


def line_bank(rng, cut_rng, name, psi, h, chi):
    singular = (chi, 0)
    exposures = []
    sharp = a.acquire(
        rng, a.joint(a.read(psi), a.read(h, site="receiver"), a.read(singular))
    )
    cut_sharp = a.acquire(
        cut_rng, a.joint(a.read(psi), a.read(h, site="cut"), a.read(singular))
    )
    den = a.sample_mean(sharp, lambda i, j, k: int(i == j == k == 2))
    estimates = []
    for fair in (False, True):
        terms = []
        for first, second in ((0, 1), (1, 0)):
            total = 0j
            for axes in product("XY", repeat=3):
                near = a.read(psi, first, axes[0])
                singlet = a.read(singular, 0, axes[2])
                remote = a.read(h, second, axes[1], "receiver")
                raw = a.acquire(rng, a.joint(near, remote, singlet, a.coin(fair)))
                blocked = a.acquire(
                    cut_rng,
                    a.joint(
                        near, a.read(h, second, axes[1], "cut"), singlet, a.coin(fair)
                    ),
                )
                value = a.sample_mean(
                    raw,
                    lambda i, j, k, c: coefficient(i, first)
                    * coefficient(j, second)
                    * coefficient(k, 0)
                    * (1, -1, 0)[c],
                )
                total += 1j ** axes.count("Y") * value / 8
                exposures.append(
                    [int(fair), first, second, "".join(axes), raw, blocked]
                )
            terms.append(total)
        for orientation, signal in (
            ("epsilon", terms[0] - terms[1]),
            ("symmetric", terms[0] + terms[1]),
        ):
            intervals = [
                a.ratio_interval(v, den, 2 * a.UNIT_ERROR)
                for v in (signal.real, signal.imag)
            ]
            estimates.append(
                dict(
                    fair_gain=fair,
                    orientation=orientation,
                    ratio=pair_complex(signal / den) if den else None,
                    intervals=intervals,
                )
            )
    truth = YUKAWA * (psi[0] * h[1] - psi[1] * h[0]) * chi
    return dict(
        name=name,
        psi=list(map(pair_complex, psi)),
        h=list(map(pair_complex, h)),
        chi=pair_complex(chi),
        reference_raw=sharp,
        blocked_reference=cut_sharp,
        reference_coincidence=den,
        shots_by_setting=exposures,
        estimates=estimates,
        prediction=pair_complex(truth),
    )


def run():
    geometry = a.transport()
    rng = np.random.default_rng(a.SEED)
    cut_rng = np.random.default_rng(a.SEED + 1)
    gain = [a.acquire(rng, a.coin(fair)) for fair in (True, False)]
    shell = [
        shell_bank(rng, cut_rng, h, control) for h in SHELL for control in (False, True)
    ]
    vacuum_phase = {
        axis: a.acquire(rng, a.read(VAC, 1, axis, "receiver")) for axis in "XY"
    }
    lines = [line_bank(rng, cut_rng, *spec) for spec in LINES]
    for panel in shell:
        if not panel["repeated_address"]:
            assert (
                panel["interval"][0]
                <= panel["independent_prediction"]
                <= panel["interval"][1]
            )
    for panel in lines:
        for estimate in panel["estimates"]:
            if estimate["orientation"] == "epsilon":
                target = [0, 0] if estimate["fair_gain"] else panel["prediction"]
                assert all(
                    low <= t <= high
                    for t, (low, high) in zip(
                        target, estimate["intervals"], strict=True
                    )
                )
    return dict(
        schema="Two-station scalar shell and oriented line 1",
        geometry=geometry,
        shots=a.SHOTS,
        seed=a.SEED,
        cut_seed=a.SEED + 1,
        unit_error=a.UNIT_ERROR,
        grid=a.GRID,
        choices=dict(radius=RADIUS, quartic=STRENGTH, yukawa=YUKAWA),
        gain_calibration=gain,
        quartic=shell,
        vacuum_phase=vacuum_phase,
        lines=lines,
        raw_encoding="Each row is the tuple of terminal address indices, followed by its count. Absent tuples have zero count.",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "results.json"
    if args.write:
        with path.open("x") as f:
            f.write("{\n")
            for i, (key, val) in enumerate(result.items()):
                if key in ("quartic", "lines"):
                    body = (
                        "[\n"
                        + ",\n".join("    " + json.dumps(r) for r in val)
                        + "\n  ]"
                    )
                else:
                    body = json.dumps(val)
                f.write(
                    "  "
                    + json.dumps(key)
                    + ": "
                    + body
                    + (",\n" if i + 1 < len(result) else "\n}\n")
                )
    else:
        compare(json.loads(path.read_text()), json.loads(json.dumps(result)))
    print(
        "Verified 7 scalar preparations, 8 complex Yukawa panels and both control families"
    )


if __name__ == "__main__":
    main()
