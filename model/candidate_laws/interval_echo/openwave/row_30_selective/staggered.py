"""Replay the alternative staggered wavelength convention on the same free branch."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo import candidate
from .analyzer import eigengroups, spatial_basis
from .experiment import N, TIMES, fly, portable, pulse, reverse_read, settings
from .source import HERE, LAYERS, STEP, Source


def alternate_geometry(source):
    adjacency = source.h < 0
    np.fill_diagonal(adjacency, False)
    depth = shortest_path(adjacency, directed=False, indices=0).astype(int)
    signs = (-1.)**depth
    degrees = adjacency.sum(axis=1)
    laplacian = np.diag(degrees) - adjacency
    potential = 6 - source.car.degree[source.sites]
    direct = signs[:, None] * source.h * signs[None, :]
    if not np.allclose(direct, candidate.KAPPA * (np.diag(potential) - laplacian/2), atol=1e-12):
        raise ArithmeticError("Staggered free operator differs from the signed path operator")
    if not np.array_equal(source.h.sum(axis=1), candidate.KAPPA * (potential - degrees)):
        raise ArithmeticError("The free generator acts differently on the uniform preparation")
    return signs


def prepare(source):
    signs = alternate_geometry(source)
    _, basis = spatial_basis(source)
    pulses = {"staggered_rest": source.compile(signs * basis[:, 0], 83157),
              "mode_26": source.compile(source.modes[:, 26], 83158)}
    return dict(ports=source.ports, layers=LAYERS, step=STEP, pulses=pulses)


def measure(source):
    signs = alternate_geometry(source)
    ks, basis = spatial_basis(source)
    staggered = signs[:, None] * basis
    spatial_groups = eigengroups(ks)
    rows = []
    for group in eigengroups(source.energies):
        projected = staggered.T @ source.modes[:, group]
        weights = [np.linalg.svd(projected[g, :], compute_uv=False)[0]**2 for g in spatial_groups]
        best = int(np.argmax(weights))
        rows.append(dict(indices=group, energy=source.energies[group[0]], best_k2=ks[spatial_groups[best][0]],
                         best_weight=weights[best], rest_weight=np.linalg.norm(projected[0, :])**2))
    alt = settings(source, HERE/"staggered_controls.json")
    base = settings(source, HERE/"controls.json")
    rest = pulse(source, alt["staggered_rest"])
    rng = np.random.default_rng(675130)
    energy_counts = []
    for i, spec in ((31, base["mode_31"]), (26, alt["mode_26"])):
        p = reverse_read(source, spec, rest)
        clicks = int(rng.binomial(N, p))
        energy_counts.append(dict(index=i, clicks=clicks, quiet=N-clicks, expected=p,
                                  analyzer_fidelity=abs(source.modes[:, i] @ pulse(source, spec))**2))
    p_initial = abs(fly(source, rest, TIMES) @ rest)**2
    h_uniform = source.h.sum(axis=1)
    shifted = signs * (source.h @ signs)
    return portable(dict(shots_per_setting=N, rest_source_fidelity=abs(staggered[:, 0] @ rest)**2,
                         energy_counts=energy_counts, free_return_minimum=min(p_initial),
                         uniform_rest_energy_variance=np.var(h_uniform),
                         staggered_rest_energy_variance=np.var(shifted), eigenspaces=rows))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_mutually_exclusive_group()
    commands.add_argument("--compile", type=Path, help="save two fresh native preparations")
    commands.add_argument("--output", type=Path, help="write a fresh experimental receipt")
    args = parser.parse_args()
    source = Source()
    result = prepare(source) if args.compile else measure(source)
    destination = args.compile or args.output
    if destination:
        with destination.open("x") as output:
            json.dump(result, output, separators=(",", ":"), allow_nan=False)
            output.write("\n")
    elif result != json.loads((HERE/"staggered.json").read_text()):
        raise SystemExit("Staggered geometry receipt differs from the archive")
    if not args.compile:
        print("Staggered rest", result["energy_counts"], "minimum return", result["free_return_minimum"])
