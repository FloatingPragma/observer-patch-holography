"""Test every on-site port relabeling of the freely released Interval Echo pair."""

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.executable_model import CELL_SEAMS, cell, federation


HERE = Path(__file__).resolve().parent
LABELS = tuple(itertools.permutations(range(3), 2))
TIMES = (0, 0.25, 1, 4)
SHOTS = 65_536


def verify_reader(background, core):
    law = candidate.native_law
    echo = law.echo()
    mask = np.zeros(law.DIM, bool)
    for memory_index, (_, contrast) in enumerate(law.LABELS, start=1):
        if contrast == 0:
            mask[memory_index * law.N_SRC + law.SRC.index((core, core))] = True
    for a, b in itertools.product((background, core), repeat=2):
        ket = np.eye(law.DIM)[:, law.SRC.index((a, b))]
        encoded = echo @ ket
        yes = (a == core and b == core)
        np.testing.assert_allclose(echo.T @ (encoded * mask), ket * yes, atol=1e-13)
        np.testing.assert_allclose(echo.T @ (encoded * ~mask), ket * (not yes), atol=1e-13)


def contrast_diagonal(sector, background, core):
    """Full contrast minus the uniform background energy, before multiplying by κ."""
    r = 6 - sector.car.degree
    s = core - background
    end_references = r[sector.aa] + r[sector.bb]
    values = s**2 * (12 - 2 * sector.adjacent) + 2 * background * s * end_references
    # Directly compare every distinct geometric case against the unsimplified native sum.
    edges = np.array(sector.car.seams)
    classes = np.column_stack((sector.adjacent, end_references))
    for index in np.unique(classes, axis=0, return_index=True)[1]:
        loads = np.full(sector.n, background, dtype=int)
        loads[[sector.aa[index], sector.bb[index]]] = core
        seam_total = np.sum((loads[edges[:, 0]] - loads[edges[:, 1]])**2)
        reference_total = r @ (loads**2 - background**2)
        if seam_total + reference_total != values[index]:
            raise ArithmeticError("Reduced contrast disagrees with native contacts and references")
    return values


def record_laws(sector, potential, reader, source):
    # Subtract the initial diagonal potential. This contributes only a common phase.
    h = sector.exchange() + sparse.diags(potential - potential[source], dtype=float)
    state = np.zeros(sector.P, complex)
    state[source] = 1
    quiet = np.ones(sector.P, bool)
    quiet[reader] = False
    rows = []
    for index, t in enumerate(TIMES):
        if index:
            state = expm_multiply((-1j * (t - TIMES[index - 1])) * h, state)
        p = abs(state)**2
        if abs(p.sum() - 1) > 1e-10:
            raise ArithmeticError("Free sector lost normalization")
        rows.append(np.r_[p[reader], p[quiet].sum()])
    return np.array(rows)


def panel(car, kappa, seed):
    sector = PairSector(car)
    source = sector.index(0, 1)
    reader = np.array([sector.index(*e) for e in CELL_SEAMS])
    baseline = kappa * contrast_diagonal(sector, 0, 1)
    # For two independent batches Hoeffding gives P(|d|>x) <= 2 exp(-N x²).
    comparisons = (len(reader) + 1) * (len(LABELS) - 1) * (len(TIMES) - 1)
    threshold = np.sqrt(np.log(2 * comparisons / .01) / SHOTS)
    rows = []
    for background, core in LABELS:
        verify_reader(background, core)
        diagonal = kappa * contrast_diagonal(sector, background, core)
        law = record_laws(sector, diagonal, reader, source)
        rng = np.random.default_rng([60617, seed, background, core])
        counts = np.array([rng.multinomial(SHOTS, p / p.sum()) for p in law])
        if (background, core) == (0, 1):
            baseline_law, baseline_counts = law, counts
        distance = np.abs(law - baseline_law).sum(axis=1) / 2
        count_difference = np.abs(counts - baseline_counts).max(axis=1) / SHOTS
        spread = float(np.ptp(diagonal - baseline))
        if spread == 0 and max(distance) > 1e-7:
            raise ArithmeticError("Equivalent generators give distinct records")
        rows.append(dict(background=background, core=core,
                         diagonal_difference_range=spread,
                         record_total_variation=distance.round(9).tolist(),
                         measured_max_address_difference=count_difference.round(9).tolist(),
                         resolved_at_any_positive_exposure=bool(max(count_difference[1:]) > threshold),
                         probabilities=law.round(9).tolist(), counts=counts.tolist()))
    return dict(cells=car.cells, kappa=kappa, times=list(TIMES), shots_per_setting=SHOTS,
                simultaneous_count_threshold=round(float(threshold), 9),
                addresses=[list(e) for e in CELL_SEAMS] + ["quiet"], arms=rows)


def run():
    return dict(candidate=candidate.NAME, panels=[
        panel(federation(0), candidate.KAPPA, 0),
        panel(federation(0), 0, 1),
        panel(cell(), candidate.KAPPA, 2),
    ])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Create a separate receipt")
    output = parser.parse_args().output
    data = run()
    if output:
        with output.open("x") as f:
            json.dump(data, f, indent=2, allow_nan=False)
            f.write("\n")
    elif data != json.loads((HERE / "results.json").read_text()):
        raise SystemExit("Port conjugation receipt differs from archive")
    for arm in data["panels"][0]["arms"]:
        print(arm["background"], arm["core"], arm["diagonal_difference_range"],
              arm["record_total_variation"], arm["resolved_at_any_positive_exposure"])
