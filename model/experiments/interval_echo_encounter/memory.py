"""A scheduled native seam write retained during free evolution, including all memory outcomes."""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier

from .__main__ import HERE, bound_project, verify

PROTOCOL = "175803436d7e7cd08d33ccd0f1972e20de18afb791d84c3901a5150e57266b34"
BITS = ((0, 0), (0, 1), (1, 0), (1, 1))


def local_components(kappa):
    blank = np.eye(11)[:, :1]
    bits = np.eye(9)[:, [law.SRC.index(loads) for loads in BITS]]
    embed_blank = np.kron(blank, bits)
    u = law.echo()
    g0, g1 = law.generator_parts()
    tagged_h = u @ (g0 + kappa*g1) @ u.T
    written = u @ embed_blank
    erase = u.T @ written
    closure = np.linalg.norm(erase - embed_blank)
    # The write uses contrast zero or magnitude one on the binary inputs.
    base = np.eye(11)
    mem = {(c, d): base[:, 1 + law.LABELS.index((c, d))] for c, d in law.LABELS}
    sectors = [(mem[0, 0] + mem[1, 0])[:, None] / np.sqrt(2)]
    for parity in (1, -1):
        sectors.append(np.column_stack([(mem[c, 1] + parity*mem[c, -1]) / np.sqrt(2)
                                        for c in (0, 1)]))
    components = []
    reconstructed = np.zeros_like(written)
    for memory in sectors:
        embed = np.kron(memory, bits)
        block = embed.T @ tagged_h @ embed
        initial = embed.T @ written
        closure = max(closure, np.linalg.norm(tagged_h @ embed - embed @ block))
        reconstructed += embed @ initial
        components.append((memory, block, initial))
    closure = max(closure, np.linalg.norm(reconstructed - written))
    if closure > 1e-10:
        raise ArithmeticError("Write, inverse or generator leaves the declared tagged binary subspace")
    other = np.asarray([law.SRC.index(pair) for pair in BITS])
    blank_term = law.seam_term(kappa)[np.ix_(other, other)]
    return components, blank_term, embed_blank.T @ erase, float(closure)


def supports(sector, edge):
    a, b = edge
    remote = [p for p in range(sector.n) if p not in edge]
    outside = np.flatnonzero(~np.isin(sector.aa, edge) & ~np.isin(sector.bb, edge))
    return [outside, np.asarray([sector.index(b, p) for p in remote]),
            np.asarray([sector.index(a, p) for p in remote]), np.asarray([sector.index(a, b)])]


def lift(matrix, positions, size):
    """Lift the local binary matrix, matching unchanged spectator ports on both sides."""
    rows, cols, values = [], [], []
    for i, j in zip(*np.nonzero(abs(matrix) > 1e-13), strict=True):
        out_mem, out_bits = divmod(int(i), 4)
        in_mem, in_bits = divmod(int(j), 4)
        if sum(BITS[out_bits]) != sum(BITS[in_bits]):
            raise ArithmeticError("Local gate changes the total load")
        rows.append(out_mem*size + positions[out_bits])
        cols.append(in_mem*size + positions[in_bits])
        values.append(np.full(len(rows[-1]), matrix[i, j]))
    shape = (matrix.shape[0]//4 * size, matrix.shape[1]//4 * size)
    if not rows:
        return sparse.csr_matrix(shape, dtype=matrix.dtype)
    return sparse.csr_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))), shape=shape)


def run(car=None):
    if hashlib.sha256((HERE / "MEMORY_PROTOCOL.md").read_bytes()).hexdigest() != PROTOCOL:
        raise ValueError("The retained-write protocol differs from its frozen edition")
    if candidate.native_law is not law or candidate.KAPPA != 8:
        raise ValueError("Different candidate seam or coupling")
    car = carrier.federation(0) if car is None else car
    sector = PairSector(car)
    seams = sorted(e for e in car.seams if max(e) < 12)
    target = seams[0]
    internal = np.asarray([sector.index(*e) for e in seams])
    starts = [sector.index(a, b) for a, b in combinations(range(12), 2) if not sector.adjacent[sector.index(a, b)]]
    source = np.zeros((sector.P, len(starts)))
    source[starts, np.arange(len(starts))] = 1.
    h = sector.hamiltonian(8)
    components, seam_term, erase, closure = local_components(8)
    pos = supports(sector, target)
    other_terms = h - lift(seam_term, pos, sector.P)
    first = expm_multiply(-1j*h, source)
    erased = lift(erase, pos, sector.P) @ first
    if np.max(abs(erased - first)) > 1e-11:
        raise ArithmeticError("Material write and inverse fail to restore the input")
    late_erase = expm_multiply(-4j*h, erased)
    # The material inverse has been checked on the whole input basis. Verify
    # the free bypass against the independent full-observer readout as well.
    erase_probs = np.mean(abs(late_erase[internal])**2, axis=1)
    erase_probs = np.append(erase_probs, 1-erase_probs.sum())
    if car.cells == 20:
        bypass_probs = np.array(json.loads((HERE / "result.json").read_text())
                                ["arms"]["contrast"]["complete_patch_distributions"]["bypass"])
    else:
        late_free = expm_multiply(-4j*h, first)
        bypass_probs = np.mean(abs(late_free[internal])**2, axis=1)
        bypass_probs = np.append(bypass_probs, 1 - bypass_probs.sum())
    if np.max(abs(erase_probs - bypass_probs)) > 1e-10:
        raise ArithmeticError("Write-erase differs from the free patch control")

    branch_records = []
    norm = np.zeros(len(starts))
    retained = np.zeros(len(seams))
    w_after = 0.
    for memory, local, write in components:
        width = memory.shape[1]
        size = sector.P
        state = lift(write, pos, size) @ first
        joint_h = sparse.kron(sparse.eye(width), other_terms, format="csr") + lift(local, pos, size)
        evolved = expm_multiply(-4j * joint_h, state)
        norm += np.sum(abs(evolved)**2, axis=0)
        branch_weight = float(np.sum(abs(state)**2) / len(starts))
        branch_records.append(branch_weight)
        for register in range(width):
            branch = state[register*size:(register+1)*size]
            projected, error = bound_project(h, branch)
            w_after += np.vdot(branch, projected).real / len(starts)
            retained += np.mean(abs(evolved[register*size + internal])**2, axis=1)
    if np.max(abs(norm - 1.)) > 1e-10:
        raise ArithmeticError("Retained state lost norm")
    projected, error = bound_project(h, source)
    w_erase = float(np.vdot(source, projected).real / len(starts))
    retained = np.append(retained, 1-retained.sum())
    disjoint = [i for i, e in enumerate(seams) if set(e).isdisjoint(target)]
    disjoint_ret, disjoint_erase = float(sum(retained[disjoint])), float(sum(erase_probs[disjoint]))
    change = disjoint_ret - disjoint_erase
    return {
        "protocol_sha256": PROTOCOL, "carrier_seam_sha256": carrier.seam_digest(car),
        "cells": car.cells, "separated_starts": len(starts), "kappa": 8, "seam": list(map(int, target)),
        "disjoint_seam_count": len(disjoint), "invariant_memory_dimensions": [1, 2, 2],
        "local_basis_and_inverse_error_below_1e_minus_10": closure < 1e-10,
        "all_start_norms_and_erase_control_within_1e_minus_10": True,
        "weight_by_zero_even_odd_contrast_memory": branch_records,
        "final_readout_order": [list(map(int, e)) for e in seams] + ["quiet"],
        "final_distribution": {"retained": retained.tolist(), "erased": erase_probs.tolist(),
                               "bypass": bypass_probs.tolist()},
        "disjoint_seam_measurement": {"retained": disjoint_ret, "erased": disjoint_erase,
                                      "retained_minus_erased": change,
                                      "absolute_change_exceeds_one_per_million": abs(change) > 1e-6},
        "load_band_immediately_after_service": {"retained": float(w_after), "erased": w_erase,
                                                "retained_minus_erased": float(w_after - w_erase),
                                                "analytic_error_bound_per_weight": error,
                                                "increase_exceeds_one_per_million":
                                                    bool(w_after - w_erase - 2*error > 1e-6)},
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path, help="write measurements to a fresh file")
    args = parser.parse_args()
    receipt = run()
    if args.record:
        with args.record.open("x") as f:
            f.write(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    else:
        verify(receipt, json.loads((HERE / "memory.json").read_text()))
    print(receipt["disjoint_seam_measurement"], receipt["load_band_immediately_after_service"])
