"""Reconstruct the home rest inventory with a clean bond clock and a protected echo."""

import math

import numpy as np

from . import peer_register_lab, sector_sampling


def spectral_lines(four_h):
    # Exact exchange-parity eigenvectors give the energy measured by the
    # 64-tick ordinary Fourier clock on a fresh copy of the sharp bond input.
    lines = np.zeros((3, 2), dtype=int)
    eye = np.eye(len(four_h), dtype=int)
    assert not np.any(four_h @ eye[0])
    for load in (1, 2):
        for j, sign in enumerate((1, -1)):
            v = eye[load] + sign * eye[3 * load]
            Hv = four_h @ v
            assert Hv[load] % 4 == 0 and np.array_equal(Hv, Hv[load] * v)
            lines[load, j] = Hv[load] // 4
    assert np.array_equal(lines, [[0, 0], [8, 9], [32, 33]])
    return lines


def acquire(tables, subset, native, home, lines, ticks):
    rng = np.random.default_rng(806173)
    slots = 256
    banks = []
    statistics = []
    for pulse, register in ((726109, 915127), (328219, 734143)):

        def observe(members, ordinary, bank_seed=pulse):
            batches = len(members)
            selected = members[
                np.arange(batches)[:, None],
                rng.integers(members.shape[1], size=(batches, slots)),
            ].ravel()
            assert len(selected) == len(ordinary)
            rows = []
            for start in (selected, ordinary):
                end, waited = sector_sampling.report_from(
                    tables, subset, native, rng, 8, start
                )
                load = home[end]
                reading = lines[load, rng.integers(2, size=len(load))]
                echo = ticks[load]
                rest = 6 * reading - 3 * ((echo - 2) // 2)
                obs = np.column_stack((reading, echo, rest))
                rows.append(obs.reshape(batches, slots, 3).mean(1))
                banks.append(
                    {
                        "pulse_seed": bank_seed,
                        "reference": len(rows) == 2,
                        "maximum_sweeps": waited,
                        "joint_echo_energy_counts": np.bincount(
                            (echo - 2) // 2 * 34 + reading, minlength=68
                        )
                        .reshape(2, 34)
                        .tolist(),
                    }
                )
            statistics.append(rows)

        peer_register_lab.run(observe=observe, pulse_seed=pulse, register_seed=register)
    # Each prepared register is one independent observation. Ordinary shots
    # are independent individually. Recover their variance from joint counts.
    prepared = np.concatenate([row[0] for row in statistics])
    hist = np.sum(
        [row["joint_echo_energy_counts"] for row in banks if row["reference"]], axis=0
    )
    occ, E = np.indices(hist.shape)
    records = np.array([E, 2 + 2 * occ, 6 * E - 3 * occ])
    nr = hist.sum()
    rm = (records * hist).sum((1, 2)) / nr
    rv = (((records - rm[:, None, None]) ** 2) * hist).sum((1, 2)) / (nr - 1)
    bound = np.array([33, 2, 195])
    # Six two-sided empirical Bernstein bounds with a union error of 0.01.
    ell = math.log(4 * 6 / 0.01)

    def width(var, n):
        return np.sqrt(2 * var * ell / n) + 7 * bound * ell / (3 * (n - 1))

    error = width(prepared.var(0, ddof=1), len(prepared)) + width(rv, nr)
    diff = prepared.mean(0) - rm
    return {
        "detectors": ["bond_energy", "roundtrip_ticks", "inferred_rest_energy"],
        "contrasts": diff.tolist(),
        "simultaneous_99_percent": np.column_stack(
            (diff - error, diff + error)
        ).tolist(),
        "trials_per_arm": int(nr),
        "prepared_registers": len(prepared),
        "slots_per_register": slots,
        "prepared_register_variance": prepared.var(0, ddof=1).tolist(),
        "raw": banks,
    }
