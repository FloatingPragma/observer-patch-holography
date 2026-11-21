"""Test each gateway geometry in the free bound-pair band and count intercell releases."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import eigsh

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_04.__main__ import evolve, verify_local_reader
from model.executable_model import carrier, federation
from .experiment import HERE, SEED, SHOTS, dispersion, lifts_and_geometry, portable


def gateway_classes(car):
    distance = shortest_path(carrier.cell().adjacency(), directed=False, unweighted=True)
    groups = {"lower_122": [], "upper_222": []}
    for cell in range(car.cells):
        ports = [p%12 for e in car.external for p in e if p//12 == cell]
        distances = sorted(distance[a,b] for i, a in enumerate(ports) for b in ports[i+1:])
        if distances == [1,2,2]:
            groups["lower_122"].append(cell)
        elif distances == [2,2,2]:
            groups["upper_222"].append(cell)
        else:
            raise ValueError("unexpected geometry on the enrolled board")
    return groups


def spatial_operators(graph, cells):
    block = graph[np.ix_(cells, cells)]
    lap = np.diag(graph.sum(axis=1)) - graph
    rest = sorted(set(range(len(graph))) - set(cells))
    cross = lap[np.ix_(cells, rest)]
    kron = lap[np.ix_(cells, cells)] - cross @ np.linalg.solve(lap[np.ix_(rest, rest)], cross.T)
    np.testing.assert_allclose(kron.sum(axis=1), 0., atol=1e-12)
    return {"induced": np.diag(block.sum(axis=1)) - block, "harmonic": kron}


def calibrate_half_gap(panels):
    lower, upper = panels["lower_122"], panels["upper_222"]
    mass = (upper["reference_energy"] - lower["reference_energy"])/2
    first = upper["classes"][1]
    beta = first["shift_from_rest"]
    c2 = beta*(2*mass + beta)/first["lambda_graph"]
    if c2 <= 0 or mass <= 0:
        return dict(rest_half_gap=mass, speed_squared=c2, calibrated=False)
    for name, data in panels.items():
        sense = -1 if name == "lower_122" else 1
        allowed = bool(data["has_held_wavenumber"] and data["sufficient_band_weight"]
                       and data["narrow_at_every_momentum"])
        resolved_curvature = False
        for row in data["classes"][1:]:
            k2 = row["lambda_graph"]
            measured = sense*row["shift_from_rest"]
            forecast = c2*k2/(np.sqrt(mass**2+c2*k2)+mass)
            tangent = c2*k2/(2*mass)
            error = abs(measured - forecast)/forecast
            separation = abs(tangent-forecast)/forecast
            row.update(signed_shift_from_midgap=measured, kg_forecast=forecast,
                       shift_error_over_forecast=error, tangent_forecast=tangent,
                       tangent_difference_over_forecast=separation)
            allowed = allowed and error < .05 and measured > 0.
            resolved_curvature |= separation > .05
        data["passes_with_common_speed_and_half_gap"] = allowed
        data["curvature_resolvable_at_enrolled_precision"] = resolved_curvature
    return dict(rest_half_gap=mass, speed_squared=c2, calibrated=True)


def spectra(h, lift, graph, groups):
    start = np.random.default_rng(SEED+4).normal(size=h.shape[0])
    energies, states = eigsh(h, k=24, ncv=84, which="SA", v0=start, tol=1e-12)
    if np.max(np.linalg.norm(h @ states-states*energies, axis=0)) > 1e-9:
        raise ArithmeticError("inaccurate class eigensolve")
    bands = {"lower_122": slice(0,8), "upper_222": slice(8,20)}
    panels = {"induced": {}, "harmonic": {}}
    for family, cells in groups.items():
        band = bands[family]
        for name, lap in spatial_operators(graph, cells).items():
            data = dispersion(lift[:, cells], np.linalg.eigh(lap), energies[band], states[:,band])
            data["reference_energy"] = energies[band][0] + data["classes"][0]["centroid_above_bottom"]
            data["has_held_wavenumber"] = len(data["classes"]) > 2
            panels[name][family] = data
    return {name: dict(calibration=calibrate_half_gap(panel), by_class=panel) for name,panel in panels.items()}


def records(sector, h, groups):
    verify_local_reader()
    addresses = [[sector.index(a,b) for a,b in sector.car.seams if a//12 == b//12 == cell]
                 for cell in range(sector.car.cells)]
    times = (0,32,128,256)
    rng = np.random.default_rng(SEED+5)
    result = {}
    for name, cells in groups.items():
        home = min(cells)
        source = np.zeros(sector.P)
        source[addresses[home][0]] = 1.
        by_exposure = []
        for exposure, state in zip(times, evolve(h, candidate.KAPPA, source, times), strict=True):
            if abs(np.linalg.norm(state)-1.) > 1e-10:
                raise ArithmeticError("free flight changed the source norm")
            detectors = []
            for indices in addresses:
                rates = abs(state[indices])**2
                probability = rates.sum()
                complete = np.r_[rates, max(0.,1-probability)]
                counts = rng.multinomial(SHOTS, complete/complete.sum())
                detectors.append(dict(probability=probability, clicks=int(sum(counts[:-1])),
                                      seam_or_quiet=counts.tolist()))
            by_exposure.append(dict(exposure=exposure, detectors_in_cell_order=detectors))
        result[name] = dict(source_cell=home, rows=by_exposure)
    return result


def run():
    sector = PairSector(federation(0))
    h = sector.hamiltonian(candidate.KAPPA)
    lifts, graph, _ = lifts_and_geometry(sector)
    groups = gateway_classes(sector.car)
    return portable(dict(candidate=candidate.NAME, seam_hash=carrier.seam_digest(sector.car),
                         groups=groups, spectral_comparison=spectra(h, lifts["uniform_seams"], graph, groups),
                         shots_per_cell_setting=SHOTS, patch_records=records(sector,h,groups)))


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, help="write a separate receipt to a new path")
    output = cli.parse_args().output
    record = run()
    if output:
        with output.open("x") as dst:
            json.dump(record,dst,indent=2,allow_nan=False)
            dst.write("\n")
    elif record != json.loads((HERE/"classes.json").read_text()):
        raise SystemExit("gateway class dispersion differs from the archive")
    for metric, panel in record["spectral_comparison"].items():
        for family,data in panel["by_class"].items():
            print(metric, family, "held momentum", data["has_held_wavenumber"],
                  "narrow lines", data["narrow_at_every_momentum"])
