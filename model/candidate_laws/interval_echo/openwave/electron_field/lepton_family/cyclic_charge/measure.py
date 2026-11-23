"""Test physical threefold rotations on the written single-charge sector."""
from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

FIELD = Path(__file__).resolve().parents[2]


def archived_modules(temp):
    work = Path(temp).resolve()
    (work / "reverse-engineering-reality").symlink_to(FIELD.parents[4], target_is_directory=True)
    files = {"charge": "experiments/interval_echo_electron_flux_binding/search.py",
             "symmetry": "experiments/interval_echo_electron_symmetry/measure.py"}
    manifest = dict(line.split("  ", 1)[::-1] for line in (FIELD / "records.sha256").read_text().splitlines())
    hashes = {}
    modules = {}
    with ZipFile(FIELD / "records.zip") as archive:
        for alias, name in files.items():
            data = archive.read(name)
            hashes[name] = sha256(data).hexdigest()
            assert hashes[name] == manifest[name]
            path = work / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            spec = importlib.util.spec_from_file_location(alias, path)
            modules[alias] = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(modules[alias])
    return modules["charge"], modules["symmetry"], hashes


def cells(symmetry):
    symmetric = broken = None
    for level in range(3):
        board = symmetry.federation(level)
        for cell_id, degree in enumerate(board.degree.reshape(-1, 12)):
            data = {"federation_level": level, "federation_cells": board.cells, "cell_id": cell_id,
                    "physical_port_degrees": degree.tolist()}
            group = sorted(symmetry.stabilizer(degree))
            order_three = [g for g in group if sorted(map(len, symmetry.cycles(g))) == [3] * 4]
            if order_three and symmetric is None:
                symmetric = data | {"degree_preserving_proper_rotations": len(group), "generator": list(order_three[0])}
            if len(group) == 1 and broken is None:
                broken = data | {"degree_preserving_proper_rotations": 1}
            if symmetric and broken:
                return symmetric, broken
    raise ValueError("Missing physical symmetric or broken-degree cell")


def backgrounds(patch, sign, copies):
    result = []
    head_vectors = {}
    for copy in copies:
        keys = set()
        for p, word in copy:
            if word not in head_vectors:
                heads = [0] * 12
                for (a, b), bit in patch.bit.items():
                    heads[a if word & bit else b] += 1
                head_vectors[word] = heads
            head = head_vectors[word].copy()
            head[p] -= sign
            keys.add(tuple(head))
        assert len(keys) == 1
        result.append(keys.pop())
    assert len(set(result)) == 3
    return result


def contacts(patch, sign, copies):
    n = len(copies[0])
    states = [state for orientation in copies for state in orientation]
    ids = {state: i for i, state in enumerate(states)}
    assert len(ids) == len(states)
    rows, cols, values = [], [], []
    transfers = np.zeros((3, 3), dtype=int)
    outside = 0
    for j, (p, word) in enumerate(states):
        moves = [((p, w), 1) for _, w in patch.ring_moves(word)]
        moves.extend((state, 8) for state in patch.charge_moves(p, word, sign))
        for state, weight in moves:
            if state not in ids:
                outside += 1
                continue
            k = ids[state]
            transfers[j // n, k // n] += 1
            rows.append(j)
            cols.append(k)
            values.append(weight)
    assert outside == 0 and not np.any(transfers - np.diag(np.diag(transfers)))
    adjacency = sparse.csr_matrix((values, (rows, cols)), shape=(len(states), len(states)), dtype=np.int64)
    assert (adjacency != adjacency.T).nnz == 0
    kin = sparse.csgraph.laplacian(adjacency).tocsr()
    return kin, np.array([p for p, _ in states]), transfers, outside


def ground(h):
    rng = np.random.default_rng(179)
    energies, v = eigsh(h.astype(float) * 0.5, k=2, which="SA", tol=1e-12, v0=rng.normal(size=h.shape[0]))
    order = np.argsort(energies)
    energies, v = energies[order], v[:, order]
    residual = float(np.max(np.linalg.norm(h @ v * 0.5 - v * energies, axis=0)))
    assert residual < 1e-8 and 0 < energies[0] < energies[1]
    return {"lowest_energy_native": float(energies[0]), "next_energy_native": float(energies[1]),
            "maximum_eigenpair_residual": residual}


def read_cell(degree, kin, ports, n, archived, original_degree, symmetric):
    v = 16 * (6 - np.asarray(degree)[ports])
    h = (kin + sparse.diags(v, dtype=np.int64)).tocsr()
    # Rotation sends the r-th n-state block to the next with identical transported basis order.
    rotated = np.roll(np.arange(3 * n), -n)
    delta = h[rotated][:, rotated] - h
    blocks = []
    for r in range(3):
        sl = slice(r * n, (r + 1) * n)
        blocks.append(h[sl, sl])
    off = h - sparse.block_diag(blocks, format="csr")
    assert off.nnz == 0
    expected = archived + sparse.diags(16 * (np.asarray(original_degree)[ports[:n]] - np.asarray(degree)[ports[:n]]),
                                        dtype=np.int64)
    assert (blocks[0] != expected).nnz == 0
    mismatches = [int((b != blocks[0]).nnz) for b in blocks]
    assert ((delta.nnz == 0) == symmetric) and (all(m == 0 for m in mismatches) == symmetric)
    spectra = [ground(b) for b in blocks]
    report = {"cross_orientation_nonzero_integer_columns": int(off.nnz),
              "integer_rotation_commutator_nonzeros": int(delta.nnz),
              "max_abs_integer_rotation_commutator_entry": int(max(abs(delta.data), default=0)),
              "integer_block_mismatches_relative_to_original": mismatches, "oriented_component_spectra": spectra}
    if symmetric:
        # H J_r = J_r H_0 over integers for each injection of a copy. Any Fourier
        # linear combination of the J_r therefore intertwines the whole H_0 exactly.
        ident = sparse.eye(n, dtype=np.int64, format="csr")
        injections = [sparse.vstack([ident if s == r else sparse.csr_matrix((n, n), dtype=np.int64)
                                     for s in range(3)], format="csr") for r in range(3)]
        differences = [int((h @ j != j @ blocks[0]).nnz) for j in injections]
        assert differences == [0, 0, 0]
        values = [s["lowest_energy_native"] for s in spectra]
        assert max(values) - min(values) < 1e-11
        report["fourier_character_blocks"] = {
            "integer_embedding_intertwiner_mismatches": differences,
            "identity": "H J_r = J_r H_0 for each r; U_k=(sum_r exp(-2 pi i k r / 3) J_r)/sqrt(3)",
            "character_exponents": [0, 1, 2], "dimension_each": n,
            "lowest_energies_native": values, "energy_ratios_same_intrinsic_level": [1, 1, 1],
            "balance_same_positive_level_exact": "1/3",
            "balance_ground_level_numerical": float(sum(values) / np.square(np.sqrt(values).sum()))}
    else:
        spread = np.ptp([s["lowest_energy_native"] for s in spectra])
        report["orientation_ground_energy_spread_native"] = float(spread)
        report["fourier_characters_are_invariant_sectors"] = False
    return report


def run(temp):
    charge, symmetry, hashes = archived_modules(temp)
    selected, broken = cells(symmetry)
    patch = charge.Patch()
    gen = selected["generator"]
    rotate = symmetry.head_rotation(patch, gen)
    records = []
    for sign in (-1, 1):
        h, states, _ = charge.charged_sector(patch, sign)
        copies = [states]
        for _ in range(2):
            copies.append([(gen[p], rotate(w)) for p, w in copies[-1]])
        assert all((gen[p], rotate(w)) == orig for (p, w), orig in zip(copies[-1], states, strict=True))
        n = len(states)
        sectors = backgrounds(patch, sign, copies)
        kin, ports, transfers, outside = contacts(patch, sign, copies)
        original_degree = charge.federation(0).degree[:12]
        sharp = (patch.faces[0][0], patch.seed)
        sid = states.index(sharp)
        records.append({"sign": sign, "component_dimension": n, "whole_dimension": 3 * n,
                        "rotated_sharp_sources_port_and_head_word": [copy[sid] for copy in copies],
                        "gauss_head_counts_minus_signed_charge_at_its_port": sectors, "background_sum": sum(sectors[0]),
                        "native_directed_move_counts_between_components": transfers.tolist(),
                        "destinations_outside_enumerated_components": outside,
                        "symmetric_cell": read_cell(selected["physical_port_degrees"], kin, ports, n, h, original_degree, True),
                        "broken_cell": read_cell(broken["physical_port_degrees"], kin, ports, n, h, original_degree, False)})
    return {"construction": "three rotated sharp reference components of the archived written single-charge extension",
            "symmetric_cell": selected, "broken_degree_control_cell": broken,
            "verified_archived_dependencies": hashes, "sectors": records}


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    args = cli.parse_args()
    with TemporaryDirectory() as temp:
        result = run(temp)
    with args.output.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
