"""Verify the carrier constraints in FLAT_BAND.md using integer arithmetic."""

import ast
import json
import struct
from pathlib import Path
from zipfile import ZipFile


PORTS = 12
PRIME = 101
SQRT5_MOD_PRIME = 45


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_cell_seams(path):
    syntax = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    assignments = [
        statement.value
        for statement in syntax.body
        if isinstance(statement, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "CELL_SEAMS"
                for target in statement.targets)
    ]
    require(len(assignments) == 1, "Expected one literal CELL_SEAMS assignment")
    seams = ast.literal_eval(assignments[0])
    require(len(seams) == len(set(seams)) == 30, "Incorrect cell seam count")
    require(all(0 <= a < b < PORTS for a, b in seams), "Invalid cell seam endpoints")
    return set(seams)


def read_i64_array(archive, name, shape, count):
    """Read the fixed little-endian integer format of the two carrier arrays."""
    with archive.open(name) as stream:
        require(stream.read(8) == b"\x93NUMPY\x01\x00", f"{name}: expected NPY 1.0")
        (header_size,) = struct.unpack("<H", stream.read(2))
        require(header_size <= 1024, f"{name}: header too long")
        header = ast.literal_eval(stream.read(header_size).decode("ascii").strip())
        require(header == {"descr": "<i8", "fortran_order": False, "shape": shape},
                f"{name}: wrong dtype, order or shape")
        body = stream.read(count * 8 + 1)
        require(len(body) == count * 8, f"{name}: incorrect payload length")
        return struct.unpack(f"<{count}q", body)


def multiply(left, right):
    return [[sum(a * b for a, b in zip(row, col)) for col in zip(*right)] for row in left]


def subtract_identity(matrix, scalar):
    return [[value - scalar * (i == j) for j, value in enumerate(row)]
            for i, row in enumerate(matrix)]


def cell_band_rows(seams):
    q = [[0] * PORTS for _ in range(PORTS)]
    for a, b in seams:
        q[a][a] += 1
        q[b][b] += 1
        q[a][b] = q[b][a] = 1
    require(all(q[i][i] == 5 for i in range(PORTS)), "Cell is not 5-regular")
    other = multiply(subtract_identity(q, 4), subtract_identity(q, 10))
    shifted = subtract_identity(q, 5)
    require(not any(map(any, multiply(other, subtract_identity(multiply(shifted, shifted), 5)))),
            "Cell matrix fails the exact spectral polynomial")
    f0 = multiply(other, shifted)
    f1 = [[-entry for entry in row] for row in other]
    # F0 + sqrt(5) F1 equals -40 times the projector onto the lowest eigenspace.
    require(sum(f0[i][i] for i in range(PORTS)) == -120
            and sum(f1[i][i] for i in range(PORTS)) == 0, "Wrong lowest-band dimension")
    return f0, f1


def rank_mod_101(rows):
    residual = [[x % PRIME for x in row] for row in rows]
    rank = 0
    for column in range(PORTS):
        pivot = next((i for i in range(rank, len(rows)) if residual[i][column]), None)
        if pivot is None:
            continue
        residual[rank], residual[pivot] = residual[pivot], residual[rank]
        row = residual[rank]
        inverse = pow(row[column], -1, PRIME)
        for i in range(rank + 1, len(rows)):
            factor = (residual[i][column] * inverse) % PRIME
            residual[i] = [(a - factor * b) % PRIME for a, b in zip(residual[i], row)]
        rank += 1
    return rank


def carrier_tripods(path, cells, seams):
    ports = PORTS * cells
    edge_count = 30 * cells + 3 * cells // 2
    with ZipFile(path) as archive:
        entries = read_i64_array(archive, "edges.npy", (edge_count, 2), 2 * edge_count)
        expected_degrees = read_i64_array(archive, "degree.npy", (ports,), ports)
    edges = list(zip(entries[::2], entries[1::2]))
    require(len(set(edges)) == edge_count, f"{path.name}: duplicate seams")
    inside = [set() for _ in range(cells)]
    tripods = [[] for _ in range(cells)]
    external = []
    degrees = [0] * ports
    used = set()
    for a, b in edges:
        require(0 <= a < b < ports, f"{path.name}: invalid seam endpoints")
        degrees[a] += 1
        degrees[b] += 1
        ca, pa = divmod(a, PORTS)
        cb, pb = divmod(b, PORTS)
        if ca == cb:
            inside[ca].add((pa, pb))
            continue
        require(a not in used and b not in used, f"{path.name}: external seams are not a matching")
        used.update((a, b))
        seam_id = len(external)
        external.append((a, b))
        tripods[ca].append((seam_id, pa))
        tripods[cb].append((seam_id, pb))
    require(all(internal == seams for internal in inside), f"{path.name}: incomplete cell")
    require(tuple(degrees) == expected_degrees, f"{path.name}: stored degrees disagree")
    require(all(len(tripod) == 3 for tripod in tripods), f"{path.name}: incorrect external tripod")
    return tripods, external


def certify_federation(tripods, external, band_rows):
    f0, f1 = band_rows
    reduced = [[(a + SQRT5_MOD_PRIME * b) % PRIME for a, b in zip(row0, row1)]
               for row0, row1 in zip(f0, f1)]
    ties = [set() for _ in external]
    killed = set()
    dependent = 0
    for cell, tripod in enumerate(tripods):
        ids, ports = zip(*tripod)
        rank = rank_mod_101([reduced[port] for port in ports])
        if rank == 3:
            killed.update(ids)
            continue
        pairs = [(i, j) for i in range(3) for j in range(i) if all(
            all(a == -b for a, b in zip(matrix[ports[i]], matrix[ports[j]]))
            for matrix in band_rows)]
        require(rank == 2 and len(pairs) == 1, f"Cell {cell}: unclassified tripod")
        i, j = pairs[0]
        ties[ids[i]].add(ids[j])
        ties[ids[j]].add(ids[i])
        killed.add(ids[3 - i - j])
        dependent += 1
    reached = set(killed)
    pending = list(killed)
    while pending:
        fresh = ties[pending.pop()] - reached
        reached.update(fresh)
        pending.extend(fresh)
    require(len(reached) == len(external), "A tie component has no forced-zero seam")
    cells = len(tripods)
    band_dimension = 3 * cells - len(external)
    symmetric_dimension = band_dimension * (band_dimension + 1) // 2
    return {
        "cells": cells,
        "external_seams": len(external),
        "independent_tripods": cells - dependent,
        "antipodal_tripods": dependent,
        "seams_reached_from_forced_zeros": len(reached),
        "lowest_band_dimension": band_dimension,
        "symmetric_dark_lower_bound_local_dimension": symmetric_dimension - 6 * cells,
        "symmetric_dark_lower_bound_seam_constraints": symmetric_dimension - 6 * cells + len(external),
        "disabled_seam_control_lowest_band_dimension": band_dimension + 1,
        "full_rank_seam_control_dimension_at_original_gap": dependent,
    }


def certificates():
    require(pow(SQRT5_MOD_PRIME, 2, PRIME) == 5, "Invalid modular square root of 5")
    package = Path(__file__).resolve().parents[2] / "contrast_memory_spherical_v1"
    seams = read_cell_seams(package / "carrier.py")
    band_rows = cell_band_rows(seams)
    results = []
    for level, cells in enumerate((20, 80, 320)):
        path = package / "carriers" / f"FEDERATION_LEVEL{level}.npz"
        tripods, external = carrier_tripods(path, cells, seams)
        results.append(certify_federation(tripods, external, band_rows))
    return {"one_cell_gap_in_J": "5 - sqrt(5)", "federations": results}


if __name__ == "__main__":
    print(json.dumps(certificates(), indent=2))
