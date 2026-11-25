"""Local tomography after a timing draw is deliberately forgotten on every launch."""

from argparse import ArgumentParser
import json

import numpy as np

from .. import detectors
from . import measure as energy

SHOTS = 2**30
WIDTHS = (
    0.0,
    1 / 16,
    1 / 4,
    1 / 2,
    1.0,
    2.0,
    4.0,
    8.0,
    64.0,
    4096.0,
    131072.0,
    1048576.0,
)


def averaged(e, mats, flight, width):
    n = len(e)
    bell = mats.transpose(0, 1, 3, 2).reshape(n, 8, 4) / np.sqrt(2)
    if width is None:
        kernel = np.eye(n)
    else:
        g = e[:, None] - e[None, :]
        kernel = np.exp(-1j * g * (flight + width / 2)) * np.sinc(
            g * width / (2 * np.pi)
        )
    rho = np.einsum("ij,ica,jcb->ab", kernel, bell, bell.conj(), optimize=True)
    np.testing.assert_allclose(rho, rho.conj().T, atol=1e-12)
    assert np.linalg.eigvalsh(rho)[0] > -1e-12
    return rho


def tomography(rho, rng, radius, shot_radius):
    projections = detectors.spin_readers()
    tensor = rho.reshape(2, 2, 2, 2)
    banks = []
    for src in projections:
        out = 2 * np.einsum("ij,iajb->ab", src, tensor)
        row = []
        for a in (0, 2, 4):
            signal = [np.trace(m @ out).real for m in projections[a : a + 2]]
            p = np.maximum([*signal, 1 - sum(signal)], 0)
            row.append(rng.multinomial(SHOTS, np.round(p / p.sum(), 13)).tolist())
        banks.append(row)
    pt = tensor.transpose(0, 3, 2, 1).reshape(4, 4)
    truth = float(np.linalg.eigvalsh(pt)[0])
    fit = float(detectors.measured_partial(banks, SHOTS)[0])
    kept = (sum(banks[0][0][:2]) + sum(banks[1][0][:2])) / (2 * SHOTS)
    # In a two-qubit system the partial transpose has at most one negative
    # eigenvalue. A uniform eigenvalue bound controls its magnitude and the
    # normalized accepted-state value after dividing by the acceptance rate.
    return dict(
        pauli_counts=banks,
        mixed_arrival=float(np.trace(rho).real),
        observed_mixed_arrival=kept,
        theory_pt_min=truth,
        measured_pt_min=fit,
        simultaneous_pt_interval=[fit - radius, fit + radius],
        accepted_negativity_lower_99=max(0, -fit - radius) / (kept + shot_radius)
        if kept
        else None,
        accepted_negativity_upper_99=min(
            0.5, max(0, radius - fit) / max(kept - shot_radius, 1e-30)
        )
        if kept
        else None,
    )


def run():
    evidence = json.loads((energy.HERE / "results.json").read_text())["cases"]["1"]
    with np.load(energy.HERE / "lines.npz") as f:
        e = f["energies_1"]
        mats = f["amplitudes_1"]
    a, b = evidence["main_pair"]
    gap = e[b] - e[a]
    # Anchor the long setting to a disjoint prior bank's measured beat,
    # before doing the randomized-wait tomography.
    beat_period = 2 * np.pi / evidence["count_spectroscopy"]["first_half_frequency"]
    widths = [*WIDTHS, beat_period]
    trials = len(widths) + 3
    # 0.5% allocated to all partial-transpose bounds, 0.5% to all arrival
    # bounds. Doubling the appointment count applies that allocation in the
    # existing six-source, three-axis instrument's formula.
    radius = float(detectors.simultaneous_bound(SHOTS, 2 * trials))
    shot_radius = float(np.sqrt(np.log(2 * trials / 0.005) / (2 * SHOTS)))
    rng = np.random.default_rng(672705)
    rows = []
    for w in widths:
        rho = averaged(e, mats, energy.FLIGHT, w)
        mid = averaged(e, mats, energy.FLIGHT + w / 2, 0)
        mid_sample = rho_mid(mid)
        rows.append(
            dict(
                width=w,
                slow_pair_visibility=float(abs(np.sinc(gap * w / (2 * np.pi)))),
                fixed_midpoint=mid_sample,
                **tomography(rho, rng, radius, shot_radius),
            )
        )
    now = averaged(e, mats, 0, 0)
    assert np.linalg.norm(now) < 1e-12
    parked = tomography(np.zeros((4, 4)), rng, radius, shot_radius)
    # Rank-one input support on every exact line proves this limiting Choi
    # state separable before numerical diagonalization or detector sampling.
    forever = tomography(averaged(e, mats, 0, None), rng, radius, shot_radius)
    sync = tomography(
        averaged(e, mats, energy.FLIGHT + 32, 0),
        np.random.default_rng(672706),
        radius,
        shot_radius,
    )
    return dict(
        synchronized_center_of_width_64=sync,
        shots_per_read=SHOTS,
        source_parity=1,
        base_flight=energy.FLIGHT,
        prior_fitted_beat_period=float(beat_period),
        simultaneous_pt_radius=radius,
        simultaneous_arrival_radius=shot_radius,
        rows=rows,
        no_flight=parked,
        infinite_wait=forever,
    )


def rho_mid(rho):
    pt = rho.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    return dict(
        mixed_arrival=float(np.trace(rho).real), pt_min=float(np.linalg.eigvalsh(pt)[0])
    )


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    s = json.dumps(run(), indent=1) + "\n"
    if parser.parse_args().write:
        (energy.HERE / "waiting.json").write_text(s)
    else:
        print(s, end="")


if __name__ == "__main__":
    main()
