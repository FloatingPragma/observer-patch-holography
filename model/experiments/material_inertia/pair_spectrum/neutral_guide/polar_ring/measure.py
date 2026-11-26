"""Return hands and remote dimer records from a closed polar ring."""

import argparse
from functools import cache
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_22_confinement.mobile_string import (
    material,
)
from model.candidate_laws.interval_echo.openwave.row_30_collective.formation import (
    native,
)
from .stock import charts

HOLDS = (1 / 64, 1 / 16, 1 / 4, 1, 4, 16, 64, 256, 1024)
FLIGHTS = (0, 16, 32, 64, 128)
SHOTS = 2**22
SEED = 191908


def microscopic():
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    eye = np.eye(law.DIM)
    blank = eye[:, :9]
    columns = blank[:, [law.SRC.index((1, 0)), law.SRC.index((0, 1))]]
    swap_block = np.array([[0.5, -0.5], [-0.5, 0.5]])
    degree = blank[:, [law.SRC.index((0, 0)), law.SRC.index((1, 0))]]
    np.testing.assert_allclose(exchange @ columns, columns @ swap_block, atol=1e-14)
    np.testing.assert_allclose(contrast @ degree, degree @ np.diag([0, 1]), atol=1e-14)
    return dict(
        commission_max_error=native(),
        contact_calibration=material.calibration_energies().tolist(),
    )


def meter():
    """The dimer's native signed echo and the private zero-reference hand relative to a separate unit load."""
    law = candidate.native_law
    basis = np.eye(law.DIM)
    U = law.echo()
    _, contrast = law.local_term_parts()
    labels = [("blank", 0), *law.LABELS]
    out = {}
    for private, gateway in law.SRC:
        relative = private - gateway
        start = basis[:, law.SRC.index((private, gateway))]
        written = (U @ start).reshape(len(labels), 9)
        np.testing.assert_allclose(
            sum(
                np.linalg.norm(written[i]) ** 2
                for i in range(1, len(labels))
                if labels[i][1] == relative
            ),
            1,
            atol=1e-14,
        )
        standard = basis[:, law.SRC.index((private, 0))]
        np.testing.assert_allclose(
            contrast @ standard, private**2 * standard, atol=1e-14
        )
        taps = np.exp(-2j * np.pi * np.arange(8) * (private**2 - 1) / 8)
        p = abs(np.fft.ifft(taps)) ** 2
        bin_ = int(np.argmax(p))
        assert p[bin_] > 1 - 1e-13 and bin_ in (7, 0, 3)
        reading = {7: 0, 0: 1, 3: 2}[bin_]
        other = reading - relative
        a = (reading - 1) + (other - 1) + 1
        b = (reading - 1) ** 2 + (other - 1) ** 2 - 1
        out[f"{private}{gateway}"] = dict(
            contrast=relative,
            degree_bin=bin_,
            private=reading,
            gateway=other,
            signed_change=a,
            squared_change=b,
            certain_count=SHOTS,
        )
    return out


def translation_blocks(n):
    mask = (1 << n) - 1
    orbits = []
    visited = set()
    for word in range(1 << n):
        if word in visited:
            continue
        orbit = []
        while word not in visited:
            visited.add(word)
            orbit.append(word)
            word = ((word << 1) & mask) | (word >> (n - 1))
        orbits.append(orbit)
    for k in range(n):
        columns, rows, data = [], [], []
        col = 0
        for orbit in orbits:
            if len(orbit) * k % n:
                continue
            rows.extend(orbit)
            columns.extend([col] * len(orbit))
            data.extend(
                np.exp(2j * np.pi * k * np.arange(len(orbit)) / n) / np.sqrt(len(orbit))
            )
            col += 1
        yield sparse.csc_matrix((data, (rows, columns)), shape=(1 << n, col))


@cache
def spectrum(n, kappa):
    orientations = material.bits(n)
    walls = np.sum(orientations != np.roll(orientations, 1, axis=1), axis=1)
    diagonal = 2 * kappa * walls + kappa * orientations.sum(axis=1)
    h = material.generate(diagonal.astype(float), n, material.H, True)
    blocks = []
    ground = None
    for k, trans in enumerate(translation_blocks(n)):
        reduced = (trans.conj().T @ h @ trans).toarray()
        energies, exact_modes = eigh(reduced, driver="evd")
        origin = np.asarray(trans.getrow(1).toarray()).ravel().conj()
        local = exact_modes.conj().T @ origin
        if k == 0:
            ground = np.asarray(trans @ exact_modes[:, 0])
        # This line is uniquely specified by its ordered eigenvalue, before any use of overlaps.
        line = 1 if k == 0 else 0
        blocks.append((trans, energies, exact_modes, local, line))
    return ground, blocks


def counts_clock(energies, weights, rng):
    banks = []
    for t in HOLDS:
        z = np.dot(weights, np.exp(-1j * energies * t))
        p = np.clip([(1 + z.real) / 2, (1 - z.imag) / 2], 0, 1)
        banks.append(rng.binomial(SHOTS, p).tolist())
    return banks


def read_clock(counts):
    est = None
    for t, (real, sine) in zip(HOLDS, counts, strict=True):
        phase = np.angle(complex(2 * real / SHOTS - 1, 2 * sine / SHOTS - 1))
        if est is not None:
            phase += 2 * np.pi * round((est * t - phase) / (2 * np.pi))
        est = float(phase / t)
    return est


def get_banks(n, kappa, rng, arrays, prefix):
    ground, blocks = spectrum(n, kappa)
    e0 = blocks[0][1][0]
    size = 1 << n
    starters = np.column_stack(
        (
            np.eye(size, 2)[:, 1],
            np.eye(size, 1)[:, 0],
            ground[np.arange(size) ^ 1],
            ground,
        )
    )
    patterns = [np.zeros_like(starters, dtype=complex) for _ in FLIGHTS]
    clock_rows = []
    band = []
    min_weights = []
    for k, (trans, e, v, local, line) in enumerate(blocks):
        weights = n * abs(local) ** 2
        np.testing.assert_allclose(weights.sum(), 1, atol=4e-13)
        counts = counts_clock(e - e0, weights, rng)
        clock_rows.append(
            dict(
                k=k,
                plus_counts=counts,
                energy=read_clock(counts),
                one_island_line_weight=float(weights[line]),
                ord_energy=float(e[line] - e0),
            )
        )
        band.append(float(e[line] - e0))
        min_weights.append(float(weights[line]))
        c = v.conj().T @ (trans.conj().T @ starters)
        for i, t in enumerate(FLIGHTS):
            patterns[i] += trans @ (v @ (np.exp(-1j * t * (e - e0))[:, None] * c))
    summary = []
    for amp, t in zip(patterns, FLIGHTS, strict=True):
        probs = abs(amp.T) ** 2
        np.testing.assert_allclose(probs.sum(axis=1), 1, atol=2e-12)
        freq = np.array([rng.multinomial(SHOTS, p / p.sum()) for p in probs])
        arrays[f"{prefix}_t{t}"] = freq
        # The two neighboring dimers at indices 0 and 1 share an owner.
        away = sum(freq[:, 1 << i] for i in range(2, n))
        summary.append(
            dict(
                time=t,
                unique_upper_other_owner=away.tolist(),
                any_unique_upper=freq[:, 1 << np.arange(n)].sum(axis=1).tolist(),
            )
        )
    return dict(
        dimer_count=n,
        kappa=kappa,
        rest_gap=band[0],
        band_interval=[min(band), max(band)],
        isolation_lower_above=6 * kappa - n / 2 - (max(band) + e0),
        minimum_one_island_weight=min(min_weights),
        ground_energy=float(e0),
        clocks=clock_rows,
        flight_summary=summary,
    )


def run(output):
    boards = charts()
    levels = {20: 0, 80: 1, 320: 2}
    raw = {}
    panels = []
    for case in boards:
        if not case["stock"]:
            continue
        n = len(case["stock"]["ring"])
        for kappa in (float(candidate.KAPPA), float(candidate.KAPPA) / 2):
            level = levels[case["cells"]]
            rng = np.random.default_rng([SEED, level, n, int(kappa)])
            prefix = f"b{case['cells']}_n{n}_g{int(kappa)}"
            panels.append(
                dict(
                    cells=case["cells"],
                    prefix=prefix,
                    **get_banks(n, kappa, rng, raw, prefix),
                )
            )
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / "word_counts.npz", **raw)
    result = dict(
        seed=SEED,
        shots_per_bank=SHOTS,
        clock_holds=HOLDS,
        flight_times=FLIGHTS,
        flight_sources=[
            "local single upper",
            "all lower",
            "local flip of the ground",
            "ground",
        ],
        microscopic=microscopic(),
        meter=meter(),
        boards=boards,
        panels=panels,
    )
    (output / "records.json").write_text(
        json.dumps(result, indent=1, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args().output)
