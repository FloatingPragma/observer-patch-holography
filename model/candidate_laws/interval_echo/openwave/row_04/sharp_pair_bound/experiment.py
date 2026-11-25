"""Integer source/boundary certificates and a complete blank-memory control."""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import isqrt
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.executable_model.carrier import HERE as CARRIER_HOME
from model.experiments.interval_echo_consensus.adjacency_bound import validate_structure
from ..__main__ import verify_local_reader
from ..topology_census import native_blocks


HERE = Path(__file__).parent
ARCHIVE = HERE / "results.json"
SOURCE = (156, 157)
TOLERANCE = 2e-10


def boundary_certificate(board):
    """Enumerate every P-to-Q departure on the complete board using integers."""
    validate_structure(board)
    adjacent = set(board.seams)
    rows, columns, separated = [], [], {}
    departures = []
    for column, (left, right) in enumerate(board.seams):
        targets = []
        for moving, fixed in ((left, right), (right, left)):
            for arrival in board.neighbours[moving]:
                if arrival == fixed:
                    continue
                pair = tuple(sorted((arrival, fixed)))
                if pair not in adjacent:
                    targets.append(pair)
        if len(targets) != len(set(targets)):
            raise ArithmeticError("A boundary destination has more than one half-rate hop")
        departures.append(len(targets))
        for pair in targets:
            rows.append(separated.setdefault(pair, len(separated)))
            columns.append(column)
    incidence = sparse.csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, columns)),
                                  shape=(len(separated), len(adjacent)))
    gram = incidence.T @ incidence
    maximum_row = int(np.asarray(gram.sum(axis=1)).max())
    maximum_departure = max(departures)
    if maximum_row > 20 or maximum_departure > 10:
        raise ArithmeticError("The complete boundary violates the proved Gram or half-degree bound")
    return {"adjacent_configurations": len(adjacent),
            "separated_boundary_configurations": len(separated),
            "maximum_gram_row_sum": maximum_row,
            "maximum_P_to_Q_departures": maximum_departure,
            "cross_norm_squared_upper": str(Fraction(maximum_row, 4)),
            "P_boundary_half_degree_upper": str(Fraction(maximum_departure, 2))}


def source_certificate(board, strength=Fraction(8)):
    if SOURCE not in board.seams:
        raise ArithmeticError("The adopted sharp source is not an actual native seam")
    moves = Counter()
    for moving, fixed in (SOURCE, SOURCE[::-1]):
        for arrival in board.neighbours[moving]:
            if arrival != fixed:
                moves[tuple(sorted((arrival, fixed)))] += 1
    degree = sum(int(board.degree[p]) for p in SOURCE)
    if (degree != 11 or len(moves) != 9 or set(moves.values()) != {1}
            or SOURCE in moves):
        raise ArithmeticError("The sharp source has a different diagonal or residual")
    return {"source": list(SOURCE), "endpoint_degrees": [int(board.degree[p]) for p in SOURCE],
            "energy": str(Fraction(degree-2, 2)-2*strength),
            "residual_norm_squared": str(Fraction(sum(v*v for v in moves.values()), 4)),
            "residual_destinations": [list(pair) for pair in sorted(moves)]}


def floor_certificate(energy, variance, cross_squared, half_degree):
    """Rational verification of the d=1 optimizer and radical enclosure."""
    a, variance = -Fraction(energy), Fraction(variance)
    b2, c = Fraction(cross_squared), Fraction(half_degree)
    if a <= 0 or variance < 0 or b2 <= 0 or c <= 0:
        raise ValueError("The source needs negative energy and positive boundary estimates")
    # Here sigma=3/2 exactly; keep the optimizer calculation explicit.
    sn, sd = isqrt(variance.numerator), isqrt(variance.denominator)
    if sn*sn != variance.numerator or sd*sd != variance.denominator:
        raise ValueError("This exact certificate requires a rational residual norm")
    sigma, crossover = Fraction(sn, sd), b2/c
    derivative_margin = c*(a-crossover)**2 - 4*variance*crossover
    if a <= crossover or derivative_margin <= 0:
        raise ValueError("The boundary interpolation maximum is not at the claimed crossover")
    scale = 10**12
    radical_square = b2 * scale**2
    lower_integer = isqrt(radical_square.numerator // radical_square.denominator)
    lower, upper = Fraction(lower_integer, scale), Fraction(lower_integer+1, scale)
    if not lower**2 <= b2 < upper**2:
        raise ArithmeticError("Invalid integer enclosure of the boundary radical")
    height = (a+crossover)**2
    floor = height/(height+(sigma+upper)**2)
    if floor <= Fraction(9, 10):
        raise ArithmeticError("The source-specific all-time adjacency floor misses 0.90")
    return {"energy_magnitude": str(a), "residual_norm": str(sigma),
            "optimizer_d": str(crossover), "derivative_margin_squared": str(derivative_margin),
            "boundary_radical_lower": str(lower), "boundary_radical_upper": str(upper),
            "exact_floor": "625/(654+12*sqrt(5))",
            "rational_floor": str(floor), "floor_decimal": float(floor)}


def finite_control():
    """Same source and reader at an exposed comparison, with all contacts active."""
    verify_local_reader()
    sector = PairSector(federation(0))
    launch = np.zeros(sector.P)
    launch[sector.index(*SOURCE)] = 1
    home_edges = [edge for edge in sector.car.seams if edge[0]//12 == edge[1]//12 == 13]
    detector = [sector.index(*edge) for edge in home_edges]
    inside = (sector.aa//12 == 13) & (sector.bb//12 == 13)
    internal = sector.aa//12 == sector.bb//12
    rows = []
    for strength in (8, 0):
        state = expm_multiply(-1j * sector.hamiltonian(strength), launch)
        probability = abs(state)**2
        bins = [*probability[detector], probability[inside & ~sector.adjacent].sum(),
                probability[~internal & sector.adjacent].sum(),
                probability[internal & ~inside & sector.adjacent].sum(),
                probability[~inside & ~sector.adjacent].sum()]
        if not np.isfinite(bins).all() or min(bins) < 0 or abs(sum(bins)-1) > TOLERANCE:
            raise ArithmeticError("Complete blank-memory census loses a terminal outcome")
        home = float(sum(bins[:30]))
        adjacent = float(home+bins[31]+bins[32])
        rows.append({"kappa": strength, "exposure": 1,
                     "census34": [float(p) for p in bins],
                     "local31": [*[float(p) for p in bins[:30]], float(sum(bins[30:]))],
                     "home_adjacent": home, "adjacent_anywhere": adjacent})
    if rows[0]["adjacent_anywhere"] <= .90 or rows[1]["adjacent_anywhere"] >= .50:
        raise ArithmeticError("The unchanged-reader contrast control fails its stated comparison")
    return {"scope": "Ideal finite-record probabilities reproducing an exposed appointment; no sampled trials or prospective credit.",
            "home_seams": [list(edge) for edge in home_edges],
            "census_complement": ["home nonadjacent", "external adjacent seam", "remote internal adjacent seam", "nonadjacent outside home"],
            "arms": rows}


def experiment():
    if candidate.KAPPA != 8:
        raise ValueError("The bound is for the adopted strength 8")
    native_blocks()
    boards = []
    for level in range(3):
        board = federation(level)
        source = source_certificate(board)
        boundary = boundary_certificate(board)
        floor = floor_certificate(source["energy"], source["residual_norm_squared"],
                                  boundary["cross_norm_squared_upper"], boundary["P_boundary_half_degree_upper"])
        boards.append({"level": level, "cells": board.cells, "ports": board.n,
                       "full_binary_sector_dimension": board.n*(board.n-1)//2,
                       "board_sha256": sha256((CARRIER_HOME/"carriers"/f"FEDERATION_LEVEL{level}.npz").read_bytes()).hexdigest(),
                       "source": source, "boundary": boundary, "bound": floor})
    return {"schema": "oph.interval_echo.sharp_pair_all_time.v1",
            "proof_sha256": sha256((HERE/"PROOF.md").read_bytes()).hexdigest(),
            "experiment_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "scope": "All-time adjacency anywhere for the original sharp blank-memory source on the three complete registered boards; no all-time home containment, physical particle identity or independent selection of kappa.",
            "law": "Every seam and degree reference active; binary H=K-16P after a common energy shift.",
            "boards": boards, "finite_control": finite_control()}


def compare(actual, archived):
    if type(actual) is not type(archived):
        raise ArithmeticError("Archive type differs from the certificate")
    if isinstance(actual, dict):
        if actual.keys() != archived.keys():
            raise ArithmeticError("Archive fields differ from the certificate")
        for key in actual:
            compare(actual[key], archived[key])
    elif isinstance(actual, list):
        if len(actual) != len(archived):
            raise ArithmeticError("Archive outcome partition differs from the certificate")
        for left, right in zip(actual, archived, strict=True):
            compare(left, right)
    elif isinstance(actual, float):
        if not np.isfinite(archived) or abs(actual-archived) > TOLERANCE:
            raise ArithmeticError("Numerical replay differs from the archive")
    elif actual != archived:
        raise ArithmeticError("Exact certificate metadata differs from the archive")


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record", type=Path, help="explicit, exclusive first-record creation")
    args = cli.parse_args()
    result = experiment()
    if args.record:
        with args.record.open("x", encoding="utf-8") as destination:
            json.dump(result, destination, indent=2, allow_nan=False)
            destination.write("\n")
    else:
        compare(result, json.loads(ARCHIVE.read_text(encoding="utf-8")))
    print("All three complete boards: adjacency at every exposure >=",
          result["boards"][0]["bound"]["floor_decimal"])


if __name__ == "__main__":
    main()
