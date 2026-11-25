"""Independent finite-flight enclosure for the published three-triplet read."""

import argparse
from itertools import product
import json
from math import factorial
from pathlib import Path

import numpy as np
from scipy import sparse

from model.candidate_laws.contrast_memory_spherical_v1.sectors import TripleSector
from model.executable_model.carrier import Carrier, cell, federation
from ..frames import EPS
from .measure import source_data

HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "flight_check.json"
STEPS, ORDER = 256, 10


def source():
    """Rebuild the accepted finite filter from its sharp launch and native spectrum."""
    frame = json.loads((HERE.parent / "frames.json").read_text())
    sector = TripleSector(cell())
    values, vectors = np.linalg.eigh(sector.hamiltonian(8).toarray())
    target = values[np.abs(values + 144 - frame["energy"]) < 1e-8].mean()
    launch = sector.rank(*frame["face_for_energy_filter"])
    x = (values - target) / (64 * np.pi)
    count = 2**20
    filt = np.exp(-1j * np.pi * (count - 1) * x) * np.sinc(count*x) / np.sinc(x)
    initial = vectors @ (filt * vectors[launch])
    admission = float(np.linalg.norm(initial)**2)
    return frame, initial / np.sqrt(admission), admission


def neighbourhood(level):
    board = federation(level)
    owners = {0} | {v // 12 for u in range(12) for v in board.neighbours[u]}
    ports = sorted(12*c + i for c in owners for i in range(12))
    indices = {v:i for i, v in enumerate(ports)}
    inside = [(indices[a], indices[b]) for a,b in board.seams if a in indices and b in indices]
    local = Carrier("owner zero and its neighbours", len(ports), inside)
    lost = board.degree[ports] - local.degree
    outside = [v for a,b in board.seams for u,v in ((a,b), (b,a))
               if u in indices and v not in indices]
    if max(lost) > 1 or len(set(outside)) != len(outside):
        raise ValueError("boundary fails the partial-matching condition")
    sector = TripleSector(local)
    return sector, lost[sector.tri].sum(axis=1)


def panel(level, kappa, frame, initial, source_error):
    sector, missing = neighbourhood(level)
    # Principal compression retains the exit cost on its diagonal. The off-block
    # channels each have one parent here, because the external seams are a matching.
    h = sector.hamiltonian(kappa) + sparse.diags(missing/2)
    home = np.flatnonzero(np.all(sector.tri < 12, axis=1))
    state = np.zeros(sector.dim, complex)
    state[home] = initial
    # Gershgorin for the three-occupant exchange Laplacian and diagonal potential.
    low, high = -6*kappa, 18.
    centre = (low+high)/2
    radius = (high-low)/2
    step = (h - centre*sparse.eye(sector.dim, format="csr")) * (-1j / STEPS)
    groups = frame["groups"]
    idx = np.asarray([sector.rank(*sorted(t)) for t in product(*groups)]).reshape(tuple(map(len,groups)))
    readers = []
    for stored_frame in frame["sharp_energy_filtered"]["frames"]:
        left, singular_values, right = np.linalg.svd(stored_frame, full_matrices=False)
        if max(abs(singular_values - 1)) > 1e-9:
            raise ValueError("published receiver frame is outside its rounding allowance")
        readers.append(left @ right)
    boundary = missing > 0
    weights = np.sqrt(missing[boundary])/2
    norm_bound = np.sqrt(max(missing))/2
    integrated = 0.
    residual = float(np.linalg.norm(weights * state[boundary]))
    panels = {}
    for tick in range(1, STEPS + 1):
        term = state
        evolved = state.copy()
        for order in range(1, ORDER+1):
            term = (step @ term) / order
            evolved += term
        state = evolved
        next_residual = float(np.linalg.norm(weights * state[boundary]))
        integrated += (residual + next_residual) / (2*STEPS)
        residual = next_residual
        if tick not in (STEPS//16, STEPS):
            continue
        time = tick/STEPS
        # Use the nearer endpoint on each half of a time step.
        # Integrating each nearest cone costs L*dt²/4 per interval.
        truncation = np.exp(radius/STEPS)*(radius/STEPS)**(ORDER+1)/factorial(ORDER+1)
        polynomial_error = np.expm1(tick * truncation)
        radius_state = (integrated + time * norm_bound * radius/(4*STEPS)
                        + (1+time*norm_bound)*polynomial_error + source_error + 1e-8)
        if abs(np.linalg.norm(state)-1) > 1e-9:
            raise ArithmeticError("flight integrator lost norm")
        code = np.einsum("ijk,ia,jb,kc->abc", state[idx], *readers, optimize=True)
        triplet = float(np.linalg.norm(code)**2)
        singlet = float(abs(np.vdot(EPS, code))**2)
        witness = singlet - triplet/3
        # For normalized states the trace distance is at most their vector
        # distance; both effects have diameter 1, as does their signed witness.
        panels[str(time)] = dict(local_triplets=triplet, local_singlet=singlet,
                                 local_witness=witness, state_distance_upper=radius_state,
                                 full_witness_lower=witness-radius_state)
    return panels


def report():
    frame, initial, admission = source()
    _, _, audited, metadata = source_data()
    overlap = np.vdot(initial, audited)
    discrepancy = np.linalg.norm(initial - audited * (overlap.conjugate()/abs(overlap)))
    source_error = float(metadata["finite_source_state_error"] + discrepancy)
    boards = {str(20*4**level): {str(k):panel(level,k,frame,initial,source_error) for k in (8,0)}
              for level in range(3)}
    if min(b["8"]["1.0"]["full_witness_lower"] for b in boards.values()) <= 0:
        raise ArithmeticError("open-board colour witness unresolved")
    return dict(source_filter_admission=admission, source_state_error_upper=source_error,
                time_step=1/STEPS,
                polynomial_degree=ORDER, numerical_allowance=1e-8, boards=boards)


def compare(read, expected):
    if isinstance(expected, dict):
        if read.keys() != expected.keys():
            raise ValueError("different receipt keys")
        for key in expected:
            compare(read[key], expected[key])
    else:
        np.testing.assert_allclose(read, expected, atol=4e-9, rtol=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = report()
    if args.write:
        with RECEIPT.open("x") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False)+"\n")
    else:
        compare(json.loads(RECEIPT.read_text()), result)
        print("Three open-board colour witnesses replayed for native and exchange-only flight.")
