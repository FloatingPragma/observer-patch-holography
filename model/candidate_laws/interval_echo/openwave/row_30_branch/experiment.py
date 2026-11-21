"""Measure free pair transfer and audit a proposed low-energy band on the 20-cell screen."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import eigsh

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_04.__main__ import evolve, verify_local_reader
from model.executable_model import carrier, federation


HERE = Path(__file__).resolve().parent
SHOTS = 262144
TIMES = (0, 16, 32, 64, 128, 256)
SEED = 302619


def lifts_and_geometry(sector):
    car = sector.car
    inner = PairSector(carrier.cell())
    _, orbitals = np.linalg.eigh(inner.hamiltonian(candidate.KAPPA).toarray())
    ground = orbitals[:, 0] * np.sign(orbitals[0, 0])
    uniform = inner.adjacent.astype(float) / np.sqrt(inner.adjacent.sum())
    lifts = {}
    for name, profile in (("uniform_seams", uniform), ("cell_ground", ground)):
        lift = np.zeros((sector.P, car.cells))
        for cell in range(car.cells):
            indices = [sector.index(12*cell + a, 12*cell + b) for a, b in zip(inner.aa, inner.bb, strict=True)]
            lift[indices, cell] = profile
        np.testing.assert_allclose(lift.T @ lift, np.eye(car.cells), atol=1e-12)
        lifts[name] = lift

    graph = np.zeros((car.cells, car.cells))
    for a, b in car.external:
        ca, cb = a // 12, b // 12
        graph[ca, cb] = graph[cb, ca] = 1.
    assert int(graph.sum()) == 2*len(car.external)
    lap = np.diag(graph.sum(axis=1)) - graph
    return lifts, graph, np.linalg.eigh(lap)


def dispersion(lift, waves, energy, states):
    spatial_energy, spatial_vectors = waves
    classes = np.unique(spatial_energy.round(8))
    # Center before calculating variances because the 20 modes are nearly degenerate.
    relative = energy - min(energy)
    overlap = states.T @ lift
    rows = []
    for value in classes:
        space = spatial_vectors[:, abs(spatial_energy - value) < 1e-7]
        weight = np.sum(abs(overlap @ space)**2, axis=1) / space.shape[1]
        retained = weight.sum()
        mean = weight @ relative / retained
        width = np.sqrt(weight @ (relative - mean)**2 / retained)
        rows.append(dict(lambda_graph=float(max(0., value)), multiplicity=space.shape[1],
                         low_band_weights=weight.tolist(), retained=retained,
                         centroid_above_bottom=mean, linewidth=width))
    rest = rows[0]["centroid_above_bottom"]
    anchor = rows[1]["centroid_above_bottom"] - rest
    for index, row in enumerate(rows):
        beta = row["centroid_above_bottom"] - rest
        row["shift_from_rest"] = beta
        row["width_over_first_shift"] = row["linewidth"] / abs(anchor)
        if index < 2 or anchor <= 0 or beta <= 0:
            continue
        ratio = row["lambda_graph"] / rows[1]["lambda_graph"]
        lower, upper = np.sqrt(ratio)*anchor, ratio*anchor
        row["allowed_shift_any_rest_mass"] = [lower, upper]
        row["relative_error_lower_bound"] = max(lower - beta, beta - upper, 0.) / beta
    return dict(classes=rows, sufficient_band_weight=min(r["retained"] for r in rows) >= .8,
                positive_shifts=all(r["shift_from_rest"] > 0 for r in rows[1:]),
                narrow_at_every_momentum=all(r["width_over_first_shift"] < .1 for r in rows))


def band_audit(sector, h, lifts, waves):
    # Oversample the requested 21 modes to resolve the entire cluster at its upper boundary.
    k = sector.car.cells
    start = np.random.default_rng(SEED).normal(size=sector.P)
    energies, states = eigsh(h, k=k+3, ncv=80, which="SA", v0=start, tol=1e-12)
    residual = np.max(np.linalg.norm(h @ states - energies*states, axis=0))
    orthogonality = np.max(abs(states.T @ states - np.eye(len(energies))))
    if residual > 1e-9 or orthogonality > 1e-10:
        raise ArithmeticError("free bound-state eigensolver failed its numerical checks")
    bandwidth = energies[k-1] - energies[0]
    gap = energies[k] - energies[k-1]
    low_states = states[:, :k]
    panels = {name: dispersion(lift, waves, energies[:k], low_states) for name, lift in lifts.items()}
    return dict(energies=energies[:k+1].tolist(), isolation_gap=gap, bandwidth=bandwidth,
                gap_above_full_width=bool(gap > bandwidth),
                min_adjacent_weight=float(np.min(np.sum(low_states[sector.adjacent]**2, axis=0))),
                spectral_residual=residual, orthogonality_error=orthogonality, lifts=panels)


def transfer(sector, graph, h):
    verify_local_reader()
    cells = [0, *map(int, np.flatnonzero(graph[0]))]
    internal_seams = [[sector.index(a, b) for a, b in sector.car.seams if a//12 == b//12 == cell]
                      for cell in cells]
    assert all(len(seams) == 30 for seams in internal_seams)
    source = np.zeros(sector.P)
    source[sector.index(0, 1)] = 1
    rng = np.random.default_rng(SEED + 1)
    panels = {}
    for name, strength, generator in (("native", candidate.KAPPA, h), ("exchange", 0, sector.hamiltonian(0))):
        rows = []
        for exposure, state in zip(TIMES, evolve(generator, strength, source, TIMES), strict=True):
            norm = np.linalg.norm(state)
            if abs(norm-1) > 1e-10:
                raise ArithmeticError("free evolution changed the norm")
            readers = []
            for cell, indices in zip(cells, internal_seams, strict=True):
                probability = abs(state[indices])**2
                total = probability.sum()
                quiet = max(0., 1 - total)
                counts = rng.multinomial(SHOTS, np.r_[probability, quiet] / (total + quiet))
                readers.append(dict(cell=cell, p=total, pair_clicks=int(sum(counts[:-1])),
                                    seam_or_quiet=counts.tolist()))
            rows.append(dict(exposure=exposure, patches=readers))
        panels[name] = rows
    # With evolution parked the distant effects have disjoint support from the source, exactly.
    parked = [float(sum(source[indices]**2)) for indices in internal_seams]
    return dict(home_and_neighbors=cells, parked=parked, panels=panels,
                detectable_transport=any(r["pair_clicks"] >= 25 for t in panels["native"][1:] for r in t["patches"][1:]))


def portable(value):
    if isinstance(value, dict):
        return {key: portable(v) for key, v in value.items()}
    if isinstance(value, list):
        return [portable(v) for v in value]
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        return round(float(value), 11) + 0.
    return value


def run():
    sector = PairSector(federation(0))
    h = sector.hamiltonian(candidate.KAPPA)
    lifts, graph, waves = lifts_and_geometry(sector)
    return portable(dict(candidate=candidate.NAME, sphere=carrier.seam_digest(sector.car),
                         cells=sector.car.cells, active_seams=len(sector.car.seams),
                         dimensions=sector.P, kappa=candidate.KAPPA, shots=SHOTS, seed=SEED,
                         band=band_audit(sector, h, lifts, waves), transfer=transfer(sector, graph, h)))


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, help="write a separate receipt to a new path")
    output = cli.parse_args().output
    result = run()
    if output is not None:
        with output.open("x") as handle:
            json.dump(result, handle, indent=2, allow_nan=False)
            handle.write("\n")
    elif result != json.loads((HERE/"results.json").read_text()):
        raise SystemExit("free bound-pair experiment differs from its archive")
    for name, panel in result["band"]["lifts"].items():
        print(name, "sufficient weight", panel["sufficient_band_weight"], "narrow lines", panel["narrow_at_every_momentum"])
    print("detectable free pair transfer", result["transfer"]["detectable_transport"])
