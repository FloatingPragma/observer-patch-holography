"""Read the static connected energy of two held material charges on a unit host."""

import argparse
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import federation
from model.experiments.interval_echo_encounter.conjugation import free_seam

HERE = Path(__file__).resolve().parent
PORTS = (1, 82, 72, 182)
BINS = 4096
PERIOD = 128
SPACING = PERIOD / BINS
SHOTS = 2**20
GRID = 2**36
SIGNS = (-1, 0, 1)
INTERIOR = ((1, 1), (2, 0), (0, 2))
SETTINGS = ("direct", "native", "conversion", "cut_conversion")
# Union over every marginal bin and both signs of the sample polarizations.
SAMPLE_ERR = np.sqrt(np.log(2 * (36 * BINS + 72) / 0.01) / (2 * SHOTS))
# Closing the rounded array by its residual changes total variation by at most K/(2*GRID).
ROUND_ERR = 3 * BINS / (2 * GRID)
# Stable spectral arithmetic over all 4096 clock hands. The largest hold is
# <202; energy rounding at 2^-30 contributes <1e-7 in state norm.
CLOCK_ERR = 3e-7


def check_board():
    board = federation(0)
    if any(board.degree[port] != 6 for port in PORTS):
        raise ValueError("The four equal-degree ports must have zero degree correction")
    for i, a in enumerate(PORTS):
        for j in range(i + 1, len(PORTS)):
            if (PORTS[j] in board.neighbours[a]) != (j == i + 1):
                raise ValueError("Expected an induced three-edge path")
    if PORTS[0] // 12 == PORTS[1] // 12:
        raise ValueError("The adjacent pair must cross an owner boundary")


def blank_term(conversion):
    h = free_seam(candidate.KAPPA).copy()
    for pair in INTERIOR[1:]:
        i, j = law.SRC.index(pair), law.SRC.index((1, 1))
        h[i, j] -= conversion / np.sqrt(2)
        h[j, i] -= conversion / np.sqrt(2)
    return h


def source_block(setting, left, right):
    conversion = 8 if setting in ("conversion", "cut_conversion") else 0
    seam = blank_term(conversion)
    if setting == "direct":
        i = law.SRC.index((left, right))
        return seam[i : i + 1, i : i + 1]
    # Source guards leave only the diagonal of either outside bond. The inner
    # bond keeps its full, unguarded two-load block in the total-two sector.
    indices = [law.SRC.index(pair) for pair in INTERIOR]
    h = seam[np.ix_(indices, indices)].copy()
    if setting == "cut_conversion":
        h.fill(0)
    for j, (a, b) in enumerate(INTERIOR):
        x, y = law.SRC.index((left, a)), law.SRC.index((b, right))
        h[j, j] += seam[x, x] + seam[y, y]
    return h


def energy_pointer(h):
    energy, vectors = np.linalg.eigh(h)
    if energy[0] < -16 or energy[-1] >= 112:
        raise ValueError("A spectral line aliases outside the declared energy interval")
    energy = np.rint(energy * 2**30) / 2**30
    projectors = np.rint((vectors * vectors[0, :]) * 2**36) / 2**36
    hands = np.arange(BINS)
    waves = np.exp(-2j * np.pi * energy[:, None] * hands / PERIOD)
    returned = projectors @ waves
    # The inverse transform recombines the coherent hold branches for each
    # final sharp material word, including interference between close lines.
    probabilities = np.abs(np.fft.ifft(returned, axis=1)) ** 2
    if abs(probabilities.sum() - 1) > CLOCK_ERR:
        raise ArithmeticError("Incomplete clock and material effects")
    return probabilities / probabilities.sum()


def line_certificate(h):
    if len(h) == 1 or np.all(h[0, 1:] == 0):
        weight = 1.0
        gap = None
    else:
        b = h[1:, 1:]
        radii = np.abs(b).sum(axis=1) - np.abs(np.diag(b))
        gap = float(min(np.diag(b) - radii) - h[0, 0])
        if gap <= 0:
            raise ArithmeticError("Ground separation could not be certified")
        weight = 1 / (1 + np.linalg.norm(h[0, 1:]) ** 2 / gap**2)
    # The nearest bin of the target line has mass >=4/pi^2. Every
    # bin >=2 spacings away has mass <=1/16 from that line.
    low = weight * 4 / np.pi**2
    high = 1 - weight + weight / 16
    certified = low - high > 2 * (SAMPLE_ERR + ROUND_ERR + CLOCK_ERR)
    return dict(
        prepared_line_weight_lower=round(float(weight), 10),
        block_gap_lower=gap,
        modal_energy_error=2 * SPACING if certified else None,
    )


def sample(probabilities, rng):
    units = np.rint(probabilities.ravel() * GRID).astype(np.int64)
    units[np.argmax(units)] += GRID - units.sum()
    if np.any(units < 0):
        raise ArithmeticError("Negative quantized probability")
    return rng.multinomial(SHOTS, units / GRID).reshape(probabilities.shape)


def decode(hist):
    marginal = hist.sum(axis=0)
    mode = int(np.argmax(marginal))
    energy = mode * SPACING
    if energy >= 112:
        energy -= PERIOD
    if hist.sum() != SHOTS:
        raise ValueError("Lost original launches")
    record = dict(energy=energy, modal_count=int(marginal[mode]))
    if len(hist) == 3:
        words = hist.sum(axis=1)
        record.update(
            central_counts=words.tolist(),
            polarization=float((words[1] - words[2]) / SHOTS),
        )
    return record


def run():
    check_board()
    rng = np.random.default_rng(230640)
    archives = {}
    rows = {}
    for setting in SETTINGS:
        spectra = {}
        energy = np.zeros((3, 3))
        for i, left in enumerate(SIGNS):
            for j, right in enumerate(SIGNS):
                h = source_block(setting, 1 + left, 1 + right)
                key = f"{setting}_{left}_{right}"
                hist = sample(energy_pointer(h), rng)
                archives[key] = hist
                result = {**decode(hist), **line_certificate(h)}
                spectra[f"{left},{right}"] = result
                energy[i, j] = result["energy"]
        u = energy + energy[1, 1] - energy[:, [1]] - energy[[1], :]
        corners = u[::2, ::2]
        all_certified = all(
            r["modal_energy_error"] is not None for r in spectra.values()
        )
        rows[setting] = dict(
            sources=spectra,
            signed_interaction=float(np.mean(corners * np.outer([-1, 1], [-1, 1]))),
            signed_error=2 * SPACING if all_certified else None,
            sign_even_interaction=float(corners.mean()),
            sign_even_error=8 * SPACING if all_certified else None,
        )
    return dict(
        ports=list(PORTS),
        launches_per_bank=SHOTS,
        clock_bins=BINS,
        clock_spacing=SPACING,
        seed=230640,
        tables=rows,
        simultaneous_bin_error=round(float(SAMPLE_ERR + ROUND_ERR + CLOCK_ERR), 10),
        polarization_error=round(float(2 * (SAMPLE_ERR + ROUND_ERR + CLOCK_ERR)), 10),
    ), archives


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path, help="create a new output directory")
    out = parser.parse_args().record
    result, archives = run()
    if out is not None:
        out.mkdir(parents=True, exist_ok=False)
        (out / "results.json").write_text(
            json.dumps(result, indent=1, allow_nan=False) + "\n"
        )
        np.savez_compressed(out / "counts.npz", **archives)
    else:
        if result != json.loads((HERE / "results.json").read_text()):
            raise ValueError("Decoded records changed")
        with np.load(HERE / "counts.npz", allow_pickle=False) as saved:
            if set(saved.files) != set(archives) or any(
                not np.array_equal(hist, saved[key]) for key, hist in archives.items()
            ):
                raise ValueError("Raw records changed")
    print("Verified 36 complete static source banks.")


if __name__ == "__main__":
    main()
