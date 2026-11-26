"""Measure source size, static contacts, bath susceptibility and signed energy fluctuations."""

import argparse
import gzip
import json
from math import comb, log, sqrt
from pathlib import Path
import numpy as np
from model.executable_model import EMITTER_SEAM, federation
from ...row_31_orbitals.written_transitions.measure import compare
from . import native

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "records.json.gz"
SHOTS, SUS_SHOTS, GRID, SEED = 2**24, 2**29, 2**42, 823173
TAIL = log(2 * 2000 / 0.01)


class Sampler:
    def __init__(self, reference=None):
        self.rng = np.random.Generator(np.random.PCG64(SEED))
        self.reference = reference
        self.grid = {}

    def read(self, key, prob, shots=SHOTS):
        p = np.asarray(prob, dtype=float)
        np.testing.assert_allclose(p.sum(), 1, atol=2e-11, rtol=0)
        if p.min() < -2e-13 or key in self.grid:
            raise ArithmeticError("Invalid probability or repeated setting")
        p = np.clip(p, 0, 1)
        v = np.rint(p * GRID).astype(np.int64)
        v[np.argmax(p)] += GRID - v.sum()
        if self.reference is not None:
            v = np.asarray(self.reference[key])
            np.testing.assert_allclose(v / GRID, p, rtol=0, atol=2e-9)
        self.grid[key] = v.tolist()
        draws = self.rng.multinomial(shots, v / GRID)
        return dict(probability=p.tolist(), counts=draws.tolist(), starts=shots)


def allowance(values, p, n):
    values = np.asarray(values, dtype=float)
    p = np.asarray(p)
    var = float(p @ ((values - p @ values) ** 2))
    return (
        sqrt(2 * var * TAIL / n)
        + 2 * np.ptp(values) * TAIL / (3 * n)
        + np.ptp(values) * len(p) / GRID
    )


def source_and_force(board, home, sampler):
    h, free = native.vacancy_operator(board, 8, {})
    state = np.zeros(len(free))
    state[home] = 1
    d = native.distances(board, home)
    panels = {}
    for age in (0, 1, 3):
        psi = native.propagation(h, state, age)
        p = abs(psi) ** 2
        p = p / p.sum()
        tag = f"{board.cells}_age{age}"
        profile = sampler.read(tag + "_source", p)
        observed = np.array(profile["counts"]) / SHOTS
        r2 = float(observed @ d**2)
        dr2 = allowance(d**2, p, SHOTS)
        row = dict(
            profile=profile,
            rms=float(np.sqrt(r2)),
            native_rms=float(np.sqrt(p @ d**2)),
            sampler_rms_interval=[sqrt(max(0, r2 - dr2)), sqrt(r2 + dr2)],
            employed={},
        )
        for shell in (1, 2, 3, 4):
            queries, aperture, rejected = native.classification(board, home, shell)
            hits, bad = float(p @ aperture), float(p @ rejected)
            forecast = float(observed @ aperture)
            da = allowance(aperture, p, SHOTS)
            for k in (8, 0):
                bins, _ = native.clock_bins(k)
                direct = bins[1] / 4
                meter = sampler.read(
                    f"{tag}_shell{shell}_k{k}", [hits, 1 - hits - bad, bad]
                )
                fraction = meter["counts"][0] / SHOTS
                dh = allowance([1, 0, 0], meter["probability"], SHOTS)
                ratio = (
                    [
                        max(0, fraction - dh) / (forecast + da),
                        (fraction + dh) / (forecast - da),
                    ]
                    if forecast > da
                    else None
                )
                row["employed"][f"{shell}_{k}"] = dict(
                    shell=shell,
                    native_contrast=k,
                    probe_ports=list(map(int, queries)),
                    energy_bins=[direct, 0, "coincident guard refusal"],
                    clock=meter,
                    mean_observed=direct * fraction,
                    mean_forecast=direct * forecast,
                    spectral_radius_here=shell,
                    independent_profile_hit=forecast,
                    ratio_to_calibrated_contact=ratio,
                    signed_residual=direct * (fraction - forecast),
                    residual_sampler_halfwidth=abs(direct) * (da + dh),
                    in_radius_window=bool(
                        r2 > 0 and 0.5 * sqrt(r2) <= shell <= 2 * sqrt(r2)
                    ),
                )
        panels[str(age)] = row
    return panels


def susceptibility(board, location, sampler):
    distance = native.distances(board, location)
    initial = native.uniform_control(board, location)
    rows = {}
    for k in (8, 0):
        for sign in (-1, 1):
            h, free = native.vacancy_operator(board, k, {location: 1 + sign})
            shells = distance[free]
            for t in (0, 0.125, 0.5, 1):
                wave = native.propagation(h, initial, t)
                sharp = np.bincount(shells, weights=abs(wave) ** 2)
                uniform = np.bincount(shells) / len(free)
                for state, p in (("mixed", uniform), ("coherent", sharp)):
                    key = f"{board.cells}_s{location}_k{k}_{sign}_{t}_{state}"
                    rows[f"{k}_{sign}_{t}_{state}"] = sampler.read(key, p, SUS_SHOTS)
    contrasts = {}
    for k in (8, 0):
        for t in (0, 0.125, 0.5, 1):
            for state in ("mixed", "coherent"):
                a, b = [rows[f"{k}_{sign}_{t}_{state}"] for sign in (1, -1)]
                scores = np.zeros(len(a["counts"]))
                scores[1] = 1
                y = 0.5 * (a["counts"][1] - b["counts"][1]) / SUS_SHOTS
                dy = (
                    sum(allowance(scores, z["probability"], SUS_SHOTS) for z in (a, b))
                    / 2
                )
                contrasts[f"{k}_{t}_{state}"] = dict(
                    first_shell_half_difference=y,
                    native_half_difference=(a["probability"][1] - b["probability"][1])
                    / 2,
                    sampler_halfwidth=dy,
                )
    return dict(port=location, cohorts=rows, readings=contrasts)


def occupations(cap_a, cap_b, overlap, M, fixed):
    one, two = cap_a - overlap, cap_b - overlap
    local = one + two + overlap
    h = M // 2
    p = np.zeros((cap_a + 1, cap_b + 1))
    for i in range(one + 1):
        for j in range(two + 1):
            for s in range(overlap + 1):
                m = i + j + s
                placements = comb(one, i) * comb(two, j) * comb(overlap, s)
                if fixed:
                    prob = (
                        comb(M - local, h - m) / comb(M, h)
                        if 0 <= h - m <= M - local
                        else 0
                    )
                else:
                    prob = (h / M) ** m * (1 - h / M) ** (local - m)
                p[i + s, j + s] += placements * prob
    return p


def fluctuations(board, home, sampler):
    d = native.distances(board, home)
    rows = {}
    for shell in (2, 3, 4):
        probe = int(np.flatnonzero(d == shell)[0])
        a, b = set(board.neighbours[home]), set(board.neighbours[probe])
        common = len(a & b)
        M = board.n - 2
        for fixed in (False, True):
            p = occupations(len(a), len(b), common, M, fixed)
            ia, ib = np.indices(p.shape)
            fa, fb = ia - len(a) / 2, ib - len(b) / 2
            key = f"{board.cells}_noise_{shell}_{int(fixed)}"
            raw = sampler.read(key, p.ravel())
            obs = np.array(raw["counts"]).reshape(p.shape) / SHOTS
            covariance = np.sum(obs * fa * fb) - np.sum(obs * fa) * np.sum(obs * fb)
            v = 0.25
            c = -v / (M - 1) if fixed else 0.0
            expected = (v - c) * common + c * len(a) * len(b)
            np.testing.assert_allclose(np.sum(p * fa * fb), expected, atol=1e-12)
            error = allowance((fa * fb).ravel(), p.ravel(), SHOTS) + allowance(
                fa.ravel(), p.ravel(), SHOTS
            ) * allowance(fb.ravel(), p.ravel(), SHOTS)
            gains = {}
            for k in (8, 0):
                _, local = native.clock_bins(k)
                gain = (local[0] - local[1]) / 4
                gains[str(k)] = dict(
                    signed_energy_variance=2 * gain**2 * float(covariance),
                    exact=2 * gain**2 * expected,
                    sampler_halfwidth=2 * gain**2 * error,
                )
            rows[f"{shell}_{'fixed' if fixed else 'product'}"] = dict(
                probe=probe,
                overlap=common,
                shape=list(p.shape),
                outcomes=raw,
                readings=gains,
            )
    return rows


def run(reference=None):
    source = Sampler(reference)
    timers = []
    for t, p in native.instruments():
        timers.append([t, source.read(f"timer_{t}", [p, 1 - p])])
    direct = {}
    for k in (8, 0):
        first, boundary = native.clock_bins(k)
        direct[str(k)] = dict(
            four_sign_clock_bins=first,
            boundary_read_times_four=boundary,
            adjacent_calibration=source.read(f"adjacent_{k}", [1, 0]),
            separated_calibration=source.read(f"separated_{k}", [0, 1]),
        )
    boards = {}
    for level in (0, 2):
        b = federation(level)
        home = EMITTER_SEAM[level][0]
        second = min(b.neighbours[home])
        shape = source_and_force(b, home, source)
        mobile = {str(s): susceptibility(b, s, source) for s in (home, second)}
        for age, row in shape.items():
            counter = row["profile"]["counts"]
            row["susceptibility_address_counts"] = [
                counter[home],
                counter[second],
                SHOTS - counter[home] - counter[second],
            ]
        boards[str(b.cells)] = dict(
            home=home,
            source=shape,
            susceptibility=mobile,
            fluctuation=fluctuations(b, home, source),
        )
    # Fresh original starts concatenate the two enrolled source addresses and
    # every terminal shell; all other captured addresses are an explicit refusal.
    for level in (0, 2):
        board = federation(level)
        panel = boards[str(board.cells)]
        panel["terminal_reader"] = native.location_instrument(board)
        mobile = panel["susceptibility"]
        locations = list(mobile)
        for age, row in panel["source"].items():
            source_p = row["profile"]["probability"]
            joined = {}
            for setting in mobile[locations[0]]["cohorts"]:
                groups = []
                for location in locations:
                    local = mobile[location]["cohorts"][setting]["probability"]
                    groups.append(source_p[int(location)] * np.array(local))
                unused = 1 - sum(source_p[int(s)] for s in locations)
                p = np.concatenate([*groups, [unused]])
                joined[setting] = source.read(
                    f"joint_{board.cells}_{age}_{setting}", p, SUS_SHOTS
                )
            row["joint_susceptibility"] = dict(
                locations=list(map(int, locations)),
                shell_counts_per_location=[
                    len(mobile[s]["cohorts"]["8_1_1_mixed"]["counts"])
                    for s in locations
                ],
                final_bin="all other captured material addresses",
                cohorts=joined,
            )
    if reference is not None and source.grid.keys() != reference.keys():
        raise ArithmeticError("Sampler inputs differ")
    return dict(
        provenance="Born probabilities from calibrated native-law reductions. PCG64-seeded multinomial counts on complete dyadic cohorts; quoted bounds are conditional sampler resolution.",
        seed=SEED,
        clock=timers,
        contact_control=direct,
        boards=boards,
        sampling_weights_over_2_to_42=source.grid,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    opts = parser.parse_args()
    old = None if opts.output else json.loads(gzip.decompress(OUTPUT.read_bytes()))
    ans = run(old["sampling_weights_over_2_to_42"] if old else None)
    if opts.output:
        with opts.output.open("xb") as f:
            f.write(
                gzip.compress(json.dumps(ans, separators=(",", ":")).encode(), mtime=0)
            )
    else:
        compare(ans, old, "free bath")
        print(
            "Replayed both native source profiles, interaction and polarization cohorts, and bath fluctuations."
        )


if __name__ == "__main__":
    main()
