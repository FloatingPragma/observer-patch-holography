"""Ionization-referenced clock and radial spectra of two freely exchanging amounts."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import EMITTER_SEAM, cell, federation
from model.experiments.interval_echo_encounter.conjugation import free_seam
from model.experiments.material_inertia.neutral_defect import generators, ordered_index

HERE = Path(__file__).parent
K = candidate.KAPPA
BINS = 128
PERIOD = 64
SHOTS = 2**20
SEED = 311941


def spectral_bounds(board):
    """Geometry-only count on the adjacent and separated configuration blocks."""
    edges = np.asarray(board.seams)
    a, b = edges.T
    common = np.array(
        [len(set(board.neighbours[i]) & set(board.neighbours[j])) for i, j in edges]
    )
    degree_sum = board.degree[a] + board.degree[b]
    records = []
    for left, right in ((1, 1), (1, 2), (2, 2)):
        distinct = int(left != right)
        gain = 2 * K * left * right
        cap = float(np.max(degree_sum / 2 + common - 1 + distinct)) - gain
        rank = (1 + distinct) * len(edges)
        outside = board.n * (board.n - 1) * (1 + distinct) // 2 - rank
        cut_moves = (1 + distinct) * np.sum(degree_sum - 2 - 2 * common)
        sep = float(cut_moves / (2 * outside))
        predictions = {}
        for n in (2, 3):
            predictions[str(n)] = [-gain / n**2, cap / n**2 + (1 - 1 / n**2) * sep]
        records.append(
            dict(
                amounts=[left, right],
                negative_levels=rank,
                bound_band=[-gain, cap],
                separation_threshold=[0, sep],
                hydrogenic_windows=predictions,
            )
        )
    return dict(cells=board.cells, ports=board.n, sectors=records)


def reference_meter(board, rng, strength=K):
    """Verify a single load's symmetric blank seam phase and its zero-energy storage."""
    law = candidate.native_law
    term, cost = law.local_term_parts()
    full = term + strength * cost
    vacuum = law.SRC.index((0, 0))
    empty_record = (law.echo()[:, vacuum].reshape(-1, law.N_SRC))[:, vacuum]
    readings = {}
    for load in (1, 2):
        state = np.eye(law.N_SRC)[:, law.SRC.index((load, 0))]
        ordered = [law.SRC.index((load, 0)), law.SRC.index((0, load))]
        energy = float((free_seam(strength)[np.ix_(ordered, ordered)] @ np.ones(2))[0])
        np.testing.assert_allclose(full @ np.kron(empty_record, state), 0, atol=2e-13)
        per_port = board.degree * energy + (6 - board.degree) * strength * load**2
        np.testing.assert_allclose(per_port, per_port[0], atol=1e-12)
        arm = np.exp(-2j * np.pi * np.arange(256) * per_port[0] / 512)
        prob = abs(np.fft.ifft(arm)) ** 2
        counts = rng.multinomial(SHOTS, prob / sum(prob))
        readings[str(load)] = dict(
            clock_bins=counts.tolist(), peak_energy=2 * int(np.argmax(counts))
        )
    return readings


def count_joint(waves, distances, rng):
    """Recombine the 128 clock histories, then read the two positions on each copy."""
    energy_bin = np.fft.fftshift(np.fft.fftfreq(BINS) * PERIOD)
    tagged = np.abs(np.fft.fftshift(np.fft.ifft(waves, axis=0), axes=0)) ** 2
    joint = np.stack(
        [tagged[:, distances == d].sum(axis=1) for d in range(1, distances.max() + 1)],
        axis=1,
    )
    np.testing.assert_allclose(joint.sum(), 1, atol=1e-10)
    p = joint.ravel().round(13)
    counts = rng.multinomial(SHOTS, p / p.sum()).reshape(joint.shape)
    window = (energy_bin >= -9) & (energy_bin <= -2)
    bound = (energy_bin >= -32) & (energy_bin <= -23)
    return dict(
        counts_by_clock_then_distance=counts.tolist(),
        held_window_clicks=int(counts[window].sum()),
        bound_window_clicks=int(counts[bound].sum()),
        adjacent_in_bound_window=int(counts[bound, 0].sum()),
    )


def all_panels(board, rng, isolated=False):
    x, y, adjacent, kinetic, _, _ = generators(board)
    distances = shortest_path(board.adjacency(), unweighted=True, directed=False)[
        x, y
    ].astype(int)
    sep = reference_meter(board, rng)
    surveys = {}
    times = np.arange(BINS) * 2 * np.pi / PERIOD
    if isolated:
        radial = np.column_stack([distances == i for i in (1, 2, 3)]).astype(float)
        radial /= np.linalg.norm(radial, axis=0)
    for name, attraction in (("registered", 4 * K), ("exchange_only", 0)):
        h = kinetic - sparse.diags(attraction * adjacent, dtype=float)
        if isolated:
            reduced = radial.T @ (h @ radial)
            np.testing.assert_allclose(h @ radial, radial @ reduced, atol=1e-13)
            energies, vectors = eigh(reduced)
            panels = {}
            for launch in (0, 1):
                waves = (
                    np.exp(-1j * times[:, None] * energies) * vectors[launch]
                ) @ vectors.T
                panels[str(launch + 1)] = count_joint(waves, np.arange(1, 4), rng)
            surveys[name] = dict(
                relative_generator=reduced.tolist(),
                energies=energies.tolist(),
                radial_probabilities_of_lines=(vectors**2).T.tolist(),
                launches=panels,
            )
        else:
            source = np.zeros(len(x))
            source[ordered_index(*EMITTER_SEAM[0], board.n)] = 1
            waves = expm_multiply(-1j * h, source, start=0, stop=times[-1], num=BINS)
            surveys[name] = count_joint(waves, distances, rng)
    if isolated:
        # An assigned three-line reference exercises the unaltered clock and radial receiver.
        lines = np.array([-32, -8, -32 / 9])
        v = np.linalg.qr(
            np.array([[1.0, 1.0, 1.0], [1.0, 0.0, -1.0], [1.0, -2.0, 1.0]]).T
        )[0]
        waves = (np.exp(-1j * times[:, None] * lines) * v[0]) @ v.T
        surveys["assigned_detector_control"] = count_joint(waves, np.arange(1, 4), rng)
    return dict(storage_referenced_uniform_single_loads=sep, settings=surveys)


def leakage_bound(band_hi):
    """Worst-case Dirichlet sidelobe contribution to the predeclared union of bins."""
    bins = np.arange(-9, -1.99, PERIOD / BINS)
    distance = np.minimum(bins - band_hi, -bins)
    # For these clocks the periodic images of the spectrum are farther than these bounds.
    return float(np.sum((BINS * np.sin(np.pi * distance / PERIOD)) ** -2))


def unlike_contact():
    seam = free_seam(K)
    law = candidate.native_law
    ids = [law.SRC.index(w) for w in ((0, 2), (2, 0))]
    # In two independent copies both excursions across this contact enter
    # coincident-address states, removed by compression to distinct addresses.
    separate = sum(seam[law.SRC.index(w), law.SRC.index(w)] for w in ((0, 1), (1, 2)))
    difference = seam[np.ix_(ids, ids)] - separate * np.eye(2)
    return dict(
        correction=difference.tolist(),
        eigenvalues=np.linalg.eigvalsh(difference).tolist(),
    )


def run():
    rng = np.random.default_rng(SEED)
    bounds = [spectral_bounds(federation(level)) for level in range(3)]
    return dict(
        shots_per_setting=SHOTS,
        seed=SEED,
        clock_bins=BINS,
        period=PERIOD,
        relative_bin_energies=np.fft.fftshift(np.fft.fftfreq(BINS) * PERIOD).tolist(),
        boards=bounds,
        filled_host_vacancy_excess=unlike_contact(),
        twenty_cell_clock_leakage_bound=leakage_bound(
            bounds[0]["sectors"][1]["bound_band"][1]
        ),
        isolated_cell=all_panels(cell(), rng, isolated=True),
        twenty_cells=all_panels(federation(0), rng),
        exchange_only_single_load_reference=reference_meter(federation(0), rng, 0),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(
            result, json.loads((HERE / "results.json").read_text()), "relative_orbitals"
        )
        print(
            "Reproduced the relative spectra, ionization standard and three-board bounds."
        )


if __name__ == "__main__":
    main()
