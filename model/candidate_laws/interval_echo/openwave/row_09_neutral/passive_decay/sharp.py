"""Complete passive-decay census on nine single sharp writes without feedback."""

import argparse
import json
import numpy as np
from scipy.linalg import eigh
from model.candidate_laws.interval_echo.openwave.row_30_collective.sharp_growth.growth import (
    instrument,
)
from model.experiments.interval_echo_consensus.signed_effects.sequential import (
    GRID_BITS,
    grid_binomial,
    grid_counts,
)
from .material import commission, hardware, edge
from .measure import HERE, SHOTS, hamiltonians, physical, compare_records

ARMS = ["joined_parent", "joined_ground", "cut_parent", "cut_ground"]


def run():
    board, d, template = hardware()
    addresses = sorted(template)
    size = 2 ** len(addresses)
    rng = np.random.default_rng(26092611)
    uniform = np.full(size, 1 / size)
    # Integer draws on the 2**-40 grid give the same counts on every platform.
    starts = np.array(
        [[grid_counts(rng, SHOTS, uniform) for _ in range(2)] for _ in range(4)]
    )
    # Every row stores preparation counts in joint, distal banks and their hits.
    banks = np.zeros((size, 4, 4), dtype=np.int64)
    far = np.kron(np.eye(8), np.array([[0.5, -0.5], [-0.5, 0.5]]))
    _, _, split, _ = hamiltonians(4)
    flip = split - np.diag(np.diag(split))
    n_ac = np.diag([int((w >> 3) & 1) * int((w >> 1) & 1) for w in range(16)])
    double = np.zeros(size, dtype=bool)
    masses = []
    for mask in range(size):
        heads = {pair: 32 + ((mask >> i) & 1) for i, pair in enumerate(addresses)}
        v = physical(board, d, heads)
        h = flip + np.diag(v - v[0])
        coupling = heads[edge(d[0][0], d[2][1])] - heads[edge(d[0][0], d[2][0])]
        far_coupling = heads[edge(d[2][0], d[3][1])] - heads[edge(d[2][0], d[3][0])]
        double[mask] = bool(coupling * far_coupling)
        bare = h - coupling * n_ac
        s = bare[::4, ::4]
        r = bare[:4, :4]
        np.testing.assert_allclose(bare, np.kron(s, np.eye(4)) + np.kron(np.eye(4), r))
        se, sv = eigh(s)
        re, rv = eigh(r)
        basis = np.kron(sv, rv)
        m = se[-1] - se[0]
        allowed = (se - se[0] > 0.5 * m) & (se[-1] - se > 1e-8)
        light = (re - re[0] > 1e-8) & (re - re[0] < 0.1 * m)
        masses.append([m, re[1] - re[0]])
        for mode, flight in enumerate((h, bare)):
            energies, frame = eigh(flight)
            transfer = (frame * np.exp(-16j * energies)) @ frame.T
            for src, level in enumerate((3, 0)):
                arm = 2 * mode + src
                psi = transfer @ basis[:, 4 * level]
                pop = abs(basis.T @ psi).reshape(4, 4) ** 2
                joint = pop[np.ix_(allowed, light)].sum()
                distal = np.vdot(psi, far @ psi).real
                n_joint, n_far = starts[arm, :, mask]
                banks[mask, arm] = [
                    n_joint,
                    grid_binomial(rng, n_joint, joint),
                    n_far,
                    grid_binomial(rng, n_far, distal),
                ]
    m = np.array(masses)
    summaries = {
        name: dict(
            zip(
                ["starts", "coincidences", "distal_starts", "distal_minus"],
                map(int, banks[:, i].sum(axis=0)),
                strict=True,
            )
        )
        for i, name in enumerate(ARMS)
    }
    for i, name in enumerate(ARMS):
        summaries[name]["two_link_coincidences"] = int(banks[double, i, 1].sum())
    return dict(
        seed=26092611,
        shots=SHOTS,
        sampling_grid_bits=GRID_BITS,
        address_order=addresses,
        arms=ARMS,
        even_odd_pointer=instrument(),
        commission=commission(board, d, template, sharp=True),
        material_histories=size,
        both_links_nonzero=int(sum(double)),
        parent_formation_range=[round(m[:, 0].min(), 9), round(m[:, 0].max(), 9)],
        receiver_threshold_range=[round(m[:, 1].min(), 9), round(m[:, 1].max(), 9)],
        totals=summaries,
    ), banks


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()
    summary, banks = run()
    if args.save:
        (HERE / "sharp.json").write_text(json.dumps(summary, indent=1) + "\n")
        np.savez_compressed(HERE / "sharp_counts.npz", banks=banks)
    else:
        compare_records(json.loads((HERE / "sharp.json").read_text()), summary)
        np.testing.assert_array_equal(
            np.load(HERE / "sharp_counts.npz")["banks"], banks
        )
        print(
            "Verified all 512 repair histories, complete launch budgets and distal records."
        )
