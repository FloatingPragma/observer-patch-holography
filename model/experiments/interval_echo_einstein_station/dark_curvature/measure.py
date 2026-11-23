"""Measure an angular-coordinate tetrahedron across the native correlated-collar arms."""
import argparse
from hashlib import sha256
from itertools import combinations
import json
import os
from pathlib import Path

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
from scipy.sparse.linalg import expm_multiply

from .. import cap, shared_variance as collar

HERE = Path(__file__).resolve().parent
TIMES = np.array([0, 1 / 128, 1 / 64, 1 / 4, 1, 2])
SHOTS = 2**48
SEED = 2026112321


def oriented(vertices):
    return np.linalg.det(vertices[..., 1:, :] - vertices[..., :1, :]) / 6


def geometry(base, source):
    coords = cap.chord_coordinates()
    excluded = {1, *source["hand_seam"], *collar.TAG}
    degree = collar.federation(0).degree
    eligible = [rail for rail in combinations(range(12), 2) if not excluded.intersection(rail)
                and abs(base[rail] + 0.5) < 1e-12 and 5 in degree[list(rail)]]
    fixed = (4, 5)
    eligible.remove(fixed)
    setups = [(fixed, *others) for others in combinations(eligible, 3)]
    volumes = np.array([oriented(coords[np.array(setup)].mean(axis=1)) for setup in setups])
    idx = np.flatnonzero(abs(volumes) >= abs(volumes).max() - 1e-10)[0]
    sign = float(np.sign(volumes[idx]))
    assert sign != 0
    setup = setups[idx]
    psi = np.zeros((base.shape[0], 4))
    for j, rail in enumerate(setup):
        psi[list(rail), j] = 2**-0.5
    # Quantize only the final volume label. Marginal port probabilities stay unrounded.
    terminals = np.vstack([coords, np.zeros(3)])
    words = np.indices((13, 13, 13, 13)).reshape(4, -1).T
    raw = sign * oriented(terminals[words])
    bins, labels = np.unique(np.round(raw, 12), return_inverse=True)
    return psi, coords, setup, sign, bins, labels


def conditional(h, psi, coords, sign, bins, labels):
    # Four resets with the same parked word, independent conditional on that word.
    p0 = psi[:12]**2
    v = p0.T @ coords
    velocity = h @ psi
    d2_probability = 2 * (velocity**2 - psi * (h @ velocity))
    d2 = d2_probability[:12].T @ coords
    edges, accel_edges = v[1:] - v[0], d2[1:] - d2[0]
    acceleration = 0.0
    for i in range(3):
        fixed = edges.copy()
        fixed[i] = accel_edges[i]
        acceleration += sign * np.linalg.det(fixed) / 6
    rows = []
    for t in TIMES:
        wave = psi if t == 0 else expm_multiply((-1j * t) * h, psi)
        p = abs(wave[:12])**2
        p = np.vstack([p, np.maximum(0, 1 - p.sum(axis=0))])
        p /= p.sum(axis=0)
        joint = np.einsum("i,j,k,l->ijkl", *p.T).ravel()
        histogram = np.bincount(labels, weights=joint, minlength=len(bins))
        expected = sign * oriented(p[:12].T @ coords)
        np.testing.assert_allclose(histogram @ bins, expected, atol=2e-12, rtol=0)
        rows.append(histogram)
    return acceleration, np.array(rows)


def run():
    base, _, _, source = collar.shared_station()
    weights, _ = collar.source_words(source)
    psi, coords, setup, sign, bins, labels = geometry(base, source)
    v0 = sign * oriented((psi[:12]**2).T @ coords)
    common = (psi, coords, sign, bins, labels)
    blank_accel, blank = conditional(base, *common)
    cond = [conditional(h, *common) for h in collar.word_generators(base, source["source_word_energies"])]
    accel, histograms = zip(*cond, strict=True)
    w = weights.reshape(3, 81)
    laws = np.concatenate([blank[None], np.einsum("sw,wtb->stb", w, histograms)])
    exact_focus = -np.r_[blank_accel, w @ accel] / v0
    np.testing.assert_allclose(exact_focus[1:], exact_focus[2], atol=1e-11, rtol=0)
    signed_population = laws @ bins / v0
    rng = np.random.default_rng(SEED)
    assert source["hand_gap"] == 1 and TIMES[-1] < np.pi
    hand = rng.binomial(SHOTS, (1 + np.stack([np.cos(TIMES), np.sin(TIMES)], axis=1)) / 2)
    hand_xy = 2 * hand / SHOTS - 1
    hand_times = np.arctan2(hand_xy[:, 1], hand_xy[:, 0])
    hand_times[0] = 0   # The release appointment is defined as the common origin.
    def stencil(ts):
        return np.linalg.solve(np.array([np.ones(3), ts[:3]**2, ts[:3]**4]), [0, 2, 0])
    coeff = stencil(hand_times)
    counts = np.array([[rng.multinomial(SHOTS, p / p.sum()) for p in arm] for arm in laws])
    record_means = counts @ bins / (SHOTS * v0)
    var_means = np.maximum(0, counts @ (bins / v0)**2 / SHOTS - record_means**2) / (SHOTS - 1)
    focus = -record_means[:, :3] @ coeff
    focus_se = np.sqrt(var_means[:, :3] @ coeff**2)
    finite = -signed_population[:, :3] @ stencil(TIMES)
    # Simultaneous two-sided Hoeffding for all 4 arms x 6 appointments.
    halfwidth = float(np.ptp(bins) * np.sqrt(np.log(2 * 4 * len(TIMES) / 0.01) / (2 * SHOTS)) / v0)
    contrast_bound = 2 * halfwidth
    # For the conditional clock, apply Hoeffding directly to each weighted, two-arm sum.
    # A union bound gives joint 99% coverage for the 3 non-fresh contrasts.
    focus_bound = float(np.ptp(bins) / v0 * np.sqrt(np.dot(coeff, coeff) * np.log(6 / 0.01) / SHOTS))
    conditional_bias = -signed_population[:, :3] @ coeff - exact_focus
    return {"protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
            "source_table_sha256": sha256((HERE.parent / "shared_variance.json").read_bytes()).hexdigest(),
            "arms": ["vacant", *collar.ARMS], "times": TIMES.tolist(), "shots": SHOTS, "seed": SEED,
            "launch_rails": setup, "initial_volume": float(v0),
            "signed_volume_outcomes": bins.tolist(), "histogram_counts": counts.tolist(),
            "unit_clock_yes_counts_cos_sin": hand.tolist(), "clock_times": hand_times.tolist(),
            "exact_signed_volume_over_release_by_arm_time": signed_population.tolist(),
            "record_signed_volume_over_release_by_arm_time": record_means.tolist(),
            "simultaneous_99pct_radius_for_record_volume_difference": contrast_bound,
            "exact_launch_focusing": exact_focus.tolist(),
            "finite_stencil_focusing": finite.tolist(), "record_focusing": focus.tolist(),
            "one_se_given_hand": focus_se.tolist(), "conditional_stencil_bias": conditional_bias.tolist(),
            "simultaneous_99pct_radius_for_focus_contrast_given_clock": focus_bound,
            "corrected_record_focus_contrasts_vs_fresh": (focus - focus[2] - conditional_bias + conditional_bias[2]).tolist(),
            "signed_score_covariance_shared_minus_fresh": float((np.log(3) / 12)**2),
            "local_proper_radius_readout": None}


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, default=HERE / "results.json")
    args = cli.parse_args()
    pkt = run()
    args.output.write_text(json.dumps(pkt, indent=2) + "\n")
    for key in ("launch_rails", "exact_launch_focusing", "record_focusing",
                "corrected_record_focus_contrasts_vs_fresh", "simultaneous_99pct_radius_for_focus_contrast_given_clock"):
        print(key, pkt[key])


if __name__ == "__main__":
    main()
