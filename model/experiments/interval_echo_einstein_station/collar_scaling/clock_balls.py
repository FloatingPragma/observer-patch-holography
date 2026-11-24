"""Protected clocks and the independently read source on complete station balls."""

import argparse
from fractions import Fraction
from itertools import combinations
import json
from math import fsum, log, pi

import numpy as np
from scipy.sparse import diags
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import spsolve

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from model.executable_model.refinement import Refinement
from model.experiments.interval_echo_space import instrument as space
from .. import shared_variance, volume
from . import measure


def clock_diameter():
    """Check the native blank generator, including both probe load sectors."""
    law = candidate.native_law
    swap = np.eye(law.N_SRC, dtype=int)[[law.SRC.index((b, a)) for a, b in law.SRC]]
    twice_exchange = np.eye(law.N_SRC, dtype=int) - swap
    phase = np.diag([8 * (a - b) ** 2 for a, b in law.SRC])
    np.testing.assert_array_equal(law.seam_term(candidate.KAPPA), twice_exchange / 2 + phase)
    np.testing.assert_array_equal(swap @ phase, phase @ swap)
    np.testing.assert_array_equal(phase.diagonal() % 2, 0)
    allowed = np.arange(12) != 4
    graph = cell().adjacency()[allowed][:, allowed]
    return int(shortest_path(graph, directed=False, unweighted=True).max())


def source_contrast():
    _, _, _, source = volume.radar_source_and_hand(federation(0))
    laws, _ = shared_variance.source_words(source)
    denominator = 2592
    counts = np.rint(laws * denominator).astype(int)
    np.testing.assert_array_equal(counts / denominator, laws)
    reverse, fresh, shared = counts
    np.testing.assert_array_equal(shared - fresh, fresh - reverse)
    np.testing.assert_array_equal(np.outer(fresh.sum(axis=1), fresh.sum(axis=0)),
                                  denominator * fresh)
    if np.any(reverse <= 0):
        raise ValueError("The reversed reference lacks full support")
    gaps = [(int(a), int(b)) for a, b in zip(fresh.flat, reverse.flat, strict=True)]
    divergence = fsum((a - b) * log(a / b) for a, b in gaps) / denominator
    lower = sum(Fraction((a - b)**2, denominator * max(a, b)) for a, b in gaps)
    if lower <= 0:
        raise ValueError("The native writer gives zero correlation contrast")
    return divergence, lower


def boundary_clock(board, inside_cells, seed):
    """Enter through a uniform cut seam and clock every swap up to the first exit."""
    ports = np.flatnonzero(np.arange(board.n) % 12 != 4)
    graph = board.adjacency()[ports][:, ports].tocsr()
    inside = inside_cells[ports // 12]
    interior = np.flatnonzero(inside)
    degrees = np.diff(graph.indptr)
    internal = graph[inside][:, inside]
    cut = degrees[inside] - np.diff(internal.indptr)
    grounded = diags(degrees[inside], dtype=float) - internal
    np.testing.assert_array_equal(grounded @ np.ones(len(interior)), cut)
    exit_times = spsolve(grounded.tocsc(), degrees[inside])
    boundary_size = int(cut.sum())
    exact = Fraction(int(degrees[inside].sum()), boundary_size)
    np.testing.assert_allclose(cut @ exit_times / boundary_size, float(exact), rtol=2e-13)
    # Entry edges with a common receiver appear once per edge. Each contact
    # after arrival is one clock step, including the contact that crosses out.
    rng = np.random.default_rng(seed)
    positions = rng.choice(np.repeat(interior, cut), size=8192)
    active = np.arange(len(positions))
    durations = np.zeros(len(positions), dtype=int)
    while len(active):
        durations[active] += 1
        ports_here = positions[active]
        outgoing = graph.indptr[ports_here] + rng.integers(degrees[ports_here])
        positions[active] = graph.indices[outgoing]
        active = active[inside[positions[active]]]
    return dict(cut_seams=boundary_size, mean_steps_exact=str(exact), trials=len(durations),
                sum_steps=int(durations.sum()), sum_squared_steps=int(durations @ durations),
                longest_completed=int(durations.max()))


def station_ball(level, contrast):
    board = federation(level)
    routes = measure.ProtectedRanges(Refinement(level, board.cells, "", np.asarray(board.external)))
    distances = routes.survey(0)
    height = int(distances.max())
    # Relative scale decreases while the number of microscopic services increases.
    radius = 4 + 2 * level
    inside = distances <= radius
    groups = np.bincount(np.flatnonzero(inside) // 4**level, minlength=20)
    population = int(inside.sum())
    denominator = board.cells * (board.cells - 1)
    globe = Fraction(population * (population - 1), denominator)
    same_ancestor = Fraction(int(groups @ (groups - 1)), denominator)
    ell = radius / height
    return dict(level=level, cells=board.cells, origin=0, radius=radius,
                farthest_station=height, stations_in_closed_ball=population,
                ancestor_counts=groups.tolist(), global_pair_probability=str(globe),
                same_ancestor_pair_probability=str(same_ancestor),
                boundary_exit_clock=boundary_clock(board, inside, 837041 + level),
                global_residual_per_coupling=contrast * float(globe) / (2 * pi * ell**4),
                ancestor_residual_per_coupling=contrast * float(same_ancestor) / (2 * pi * ell**4))


def local_gram_volumes():
    """Reconstruct local record-space volume from fresh short-dwell port clicks."""
    rng = np.random.default_rng(137081)
    volumes = []
    readings = []
    with np.load(space.HERE / "records.npz", allow_pickle=False) as records:
        for key in sorted(k for k in records if k.endswith("_home")):
            # Each column is one anonymous port launch. The thirteenth output is quiet.
            launches = []
            for column in records[key].T:
                cuts = np.rint(2**32 * np.cumsum(column / column.sum())).astype(np.int64)
                draws = rng.integers(2**32, size=2**20, dtype=np.uint32)
                outcomes = np.searchsorted(cuts, draws, side="right")
                launches.append(np.bincount(outcomes, minlength=len(column)))
            counts = np.array(launches).T
            frame, edges = space.slow_frame(counts[:12])
            if frame.shape != (12, 3) or edges != 30:
                raise ValueError("The local clicks do not resolve the repair triplet")
            gram = 4 * frame @ frame.T
            site = []
            for triple in combinations(range(12), 3):
                face = gram[np.ix_(triple, triple)]
                if min(face[np.triu_indices(3, 1)]) > 0:
                    site.append(float(np.sqrt(np.linalg.det(face)) / 6))
            if len(site) != 20:
                raise ValueError("The positive Gram triples do not resolve twenty faces")
            volumes.extend(site)
            readings.append(dict(bank=key, counts_by_launch=counts.T.tolist()))
    return dict(shots_per_launch=2**20, probability_ticks=2**32,
                tetrahedra=len(volumes), minimum_volume=min(volumes),
                maximum_volume=max(volumes), home_readings=readings)


def run():
    diameter = clock_diameter()
    charge, lower = source_contrast()
    density = Fraction(1, 28800 * (diameter + 1)**2)
    return dict(protected_cell_diameter=diameter, ideal_source_contrast_nats=charge,
                rational_contrast_lower=str(lower), station_fraction_lower_per_radius_squared=str(density),
                global_charge_lower_per_radius_fourth=str(lower * density**2 / 2),
                ancestor_charge_lower_per_radius_fourth=str(lower * density**2 / 40),
                local_gram=local_gram_volumes(),
                scope="Exact source-law populations on full routed balls. Uniform bounds assume "
                      "uniform midpoint incidence, the protected cell template and one station per cell.",
                bundled_boards=[station_ball(level, charge) for level in range(3)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = run()
    dest = measure.HERE / "clock_balls.json"
    if parser.parse_args().write:
        dest.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        from ..collar_pairing.measure import compare
        compare(data, json.loads(dest.read_text()))
    print([(r["cells"], r["stations_in_closed_ball"], r["global_residual_per_coupling"])
           for r in data["bundled_boards"]])
