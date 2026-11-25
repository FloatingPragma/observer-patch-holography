"""Same-source open-board binding with exact integer boundary certificates."""
import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import isqrt
import json
from pathlib import Path

import numpy as np
from scipy import linalg, sparse

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.executable_model import cell, federation
from model.executable_model.carrier import HERE as CARRIER_HOME
from ..frames import EPS


HERE = Path(__file__).parent
ARCHIVE = HERE/"results.json"
FRAMES = HERE.parent/"frames.json"
TOLERANCE = 2e-9


def root_upper(value, scale=1 << 100):
    value = Fraction(value)
    return Fraction(isqrt(value.numerator*scale**2//value.denominator)+1, scale)


def dyadic_spectral_audit(h, energies, vectors, target):
    """Exact integer residual and Gram checks on the numerical spectral frame."""
    n, scale = len(h), 1 << 80
    twice = np.rint(2*h).astype(np.int64)
    if not np.array_equal(2*h, twice):
        raise ArithmeticError("The source matrix is not exactly half-integral")
    frame = np.array([[round(Fraction.from_float(float(v))*scale) for v in row]
                      for row in vectors], dtype=object)
    levels = np.array([round(Fraction.from_float(float(v))*scale) for v in energies], dtype=object)
    gram = frame.T@frame-scale**2*np.eye(n, dtype=object)
    eta = Fraction(max(sum(abs(int(v)) for v in row) for row in gram), scale**2)
    residual = (twice.astype(object)@frame)*scale-2*frame*levels
    rho = root_upper(Fraction(sum(int(v)**2 for v in residual.flat), 4*scale**4))
    hnorm = Fraction(int(np.max(abs(twice).sum(axis=1))), 2)
    lnorm = Fraction(max(abs(int(v)) for v in levels), scale)
    if eta >= Fraction(1, 2):
        raise ArithmeticError("Dyadic spectral frame is not close to an orthogonal matrix")
    delta = rho+(hnorm+lnorm)*eta
    selected = np.flatnonzero(target)
    other = np.flatnonzero(~target)
    gap = min(Fraction(abs(int(levels[i])-int(levels[j])), scale) for i in selected for j in other)
    width = Fraction(max(levels[selected])-min(levels[selected]), scale)+2*delta
    selected_upper = Fraction(max(levels[selected]), scale)+delta
    if gap <= 2*delta or width >= Fraction(1, 10**6) or selected_upper > Fraction(201, 2):
        raise ArithmeticError("The isolated source cluster is not resolved")
    frame_error = eta+Fraction(n, 2*scale)
    projector_error = 2*delta/(gap-2*delta)+frame_error*(2+frame_error)
    return {"dyadic_bits": 80, "exact_Gram_operator_error_upper": str(eta),
            "exact_residual_frobenius_upper": str(rho),
            "nearby_orthogonal_spectral_operator_error": str(delta),
            "selected_energy_upper": str(selected_upper),
            "target_separation_lower": str(gap-2*delta), "target_cluster_width_upper": str(width),
            "target_projector_error_upper": str(projector_error)}, float(delta), float(projector_error), float(frame_error)


def source_data():
    saved = json.loads(FRAMES.read_text())
    sector = BinarySector(cell(), 3)
    h = sector.hamiltonian(8, False)
    # The proof's dyadic Gram and residual numerators refer to a particular
    # numerical witness. Other valid bases of a degenerate eigenspace give
    # different numerators. Certify the deposited frame against the live h.
    with np.load(HERE / "spectral_frame.npz", allow_pickle=False) as frame:
        energies, vectors = frame["energies"], frame["vectors"]
    if energies.shape != (len(h),) or vectors.shape != h.shape:
        raise ArithmeticError("Spectral frame dimensions differ from the native three-load sector")
    target = abs(energies-saved["energy"]) < 1e-8
    if target.sum() != 4:
        raise ArithmeticError("The published trimer source line differs")
    face = sector.agreed == 3
    if not np.array_equal(2*h[np.ix_(face, face)], 201*np.eye(20)):
        raise ArithmeticError("The exact isolated face compression differs")
    sharp = sector.index[tuple(saved["face_for_energy_filter"])]
    projection = vectors[:, target]@vectors[sharp, target]
    ideal_admission = float(np.vdot(projection, projection).real)
    ideal = projection/np.sqrt(ideal_admission)
    audit, delta, projector_error, frame_error = dyadic_spectral_audit(h, energies, vectors, target)
    arms, step = saved["finite_filter"]["timing_arms"], saved["finite_filter"]["step"]
    energy = float(energies[target].mean())
    z = (energies-energy)*step/(2*np.pi)
    filt = np.exp(-1j*np.pi*(arms-1)*z)*np.sinc(arms*z)/np.sinc(z)
    raw = vectors@(vectors[sharp]*filt)
    admission = float(np.vdot(raw, raw).real)
    finite = raw/np.sqrt(admission)
    # Time-average Kraus error: mean delay times operator residual, plus
    # the spectral-frame replacement and a disclosed scalar/dot tolerance.
    kraus_error = ((arms-1)*step/2)*delta+frame_error*(2+frame_error)+1e-10
    if np.sqrt(admission) <= kraus_error:
        raise ArithmeticError("Finite filter admission is unresolved")
    finite_error = float(2*kraus_error/(np.sqrt(admission)-kraus_error))
    ideal_error = float(2*projector_error/(np.sqrt(ideal_admission)-projector_error))
    distance = float(np.linalg.norm(finite-ideal*np.vdot(ideal, finite)))
    distance_upper = float(distance+finite_error+ideal_error+1e-11)
    groups = saved["groups"]
    frames, corrections = [], []
    for old in saved["sharp_energy_filtered"]["frames"]:
        old = np.array(old)
        u, _, vt = linalg.svd(old, full_matrices=False)
        fixed = u@vt
        corrections.append(float(linalg.norm(fixed-old)))
        frames.append(fixed)
    index = np.array([sector.index[tuple(sorted(w))] for w in product(*groups)]).reshape(tuple(map(len, groups)))
    code = np.einsum("ijk,ia,jb,kc->abc", ideal[index], *frames, optimize=True)
    triplets = float(np.vdot(code, code).real)
    singlet = float(abs(np.vdot(EPS, code))**2)
    witness = singlet-triplets/3
    witness_error = float(ideal_error+8*sum(corrections)+1e-10)
    admission_interval = [(np.sqrt(admission)-kraus_error)**2, (np.sqrt(admission)+kraus_error)**2]
    metadata = {"home_owner": 0, "sharp_face": saved["face_for_energy_filter"],
                "target_energy": energy, "target_dimension": int(target.sum()),
                "ideal_source_admission": ideal_admission, "finite_filter_admission": admission,
                "finite_filter_admission_interval": [float(v) for v in admission_interval],
                "filter_refusal": 1-admission, "timing_arms": arms, "timing_step": step,
                "spectral_audit": audit, "scalar_and_dot_tolerance": 1e-10,
                "finite_source_state_error": finite_error, "ideal_source_state_error": ideal_error,
                "finite_vs_ideal_trace_distance_center": distance,
                "finite_vs_ideal_trace_distance_upper": distance_upper,
                "receiver_groups": groups, "rounded_frame_polar_corrections": corrections,
                "initial_ideal_local_triplet_probability": triplets,
                "initial_ideal_local_singlet_probability": singlet,
                "initial_ideal_witness": witness, "initial_witness_error": witness_error,
                "no_face_or_connected_flag_before_colour_read": True}
    return sector, ideal, finite, metadata


def connected_boundary(board):
    connected = {tuple(sorted((a, middle, b))) for middle in range(board.n)
                 for a, b in combinations(board.neighbours[middle], 2)}
    rows, columns, outside, maximum = [], [], {}, 0
    for j, word in enumerate(sorted(connected)):
        occupied, departures = set(word), []
        for moving in word:
            for arrival in board.neighbours[moving]:
                if arrival not in occupied:
                    target = tuple(sorted((occupied-{moving})|{arrival}))
                    if target not in connected:
                        departures.append(target)
        if len(departures) != len(set(departures)):
            raise ArithmeticError("A connected-to-fragmented hop has duplicate parents")
        maximum = max(maximum, len(departures))
        for target in departures:
            rows.append(outside.setdefault(target, len(outside)))
            columns.append(j)
    m = sparse.csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, columns)),
                           shape=(len(outside), len(connected)))
    gram_row = int(np.asarray((m.T@m).sum(axis=1)).max())
    if maximum > 12 or gram_row > 28:
        raise ArithmeticError("The full-board cross-boundary certificate fails")
    return {"connected_triples": len(connected), "fragmented_boundary_triples": len(outside),
            "max_cross_departures": maximum, "max_integer_Gram_row_sum": gram_row,
            "cross_norm_squared_upper": str(Fraction(gram_row, 4)),
            "P_boundary_half_degree_upper": str(Fraction(maximum, 2)),
            "fragment_compression_energy_floor": -16}


def action(board, source):
    result = defaultdict(complex)
    for word, amplitude in source.items():
        occupied = set(word)
        pairs = sum(q in occupied for p in word for q in board.neighbours[p])//2
        moves = [tuple(sorted((occupied-{p})|{q})) for p in word for q in board.neighbours[p] if q not in occupied]
        result[word] += amplitude*(len(moves)/2-16*pairs)
        for target in moves:
            result[target] -= amplitude/2
    return result


def seam_energies(board, source):
    incident = sorted({tuple(sorted((p, q))) for word in source for p in word for q in board.neighbours[p]})
    records = []
    for edge in incident:
        both, odd, visited = 0., 0., set()
        for word, amplitude in source.items():
            if all(p in word for p in edge):
                both += abs(amplitude)**2
            if len(set(word).intersection(edge)) != 1 or word in visited:
                continue
            target = tuple(sorted(set(word).symmetric_difference(edge)))
            visited.update((word, target))
            odd += abs(amplitude-source.get(target, 0))**2/2
        records.append({"edge": list(edge), "odd_both_other": [float(odd), float(both), float(1-odd-both)]})
    kinetic = float(sum(row["odd_both_other"][0] for row in records))
    pairs = float(sum(row["odd_both_other"][1] for row in records))
    return {"nontrivial_edges": records, "remaining_host_edges": len(board.seams)-len(records),
            "host_odd_both_other": [0., 0., 1.], "same_source_shifted_energy_kappa8": kinetic-16*pairs,
            "same_source_energy_exchange_only": kinetic,
            "native_absolute_energy": kinetic-16*pairs+144,
            "native_disconnected_compression_floor_absolute": 128,
            "exchange_only_disconnected_compression_floor": 0}


def run():
    sector, ideal, finite, source_info = source_data()
    sigma = Fraction(1000001, 1000000)
    a, crossover, c = Fraction(26), Fraction(7, 6), Fraction(6)
    if c*(a-crossover)**2 <= 4*sigma**2*crossover:
        raise ArithmeticError("The boundary-interpolation maximum is unresolved")
    height = (a+crossover)**2
    floor = height/(height+(sigma+root_upper(7))**2)
    prepared_floor = float(floor)-source_info["finite_vs_ideal_trace_distance_upper"]
    witness_bounds = [{"exposure": t, "conditional_witness_lower": source_info["initial_ideal_witness"]
                      -source_info["initial_witness_error"]-float(sigma)*t-source_info["finite_vs_ideal_trace_distance_upper"]}
                     for t in (0., 1/16, .07)]
    if prepared_floor <= .975 or min(row["conditional_witness_lower"] for row in witness_bounds) <= 0:
        raise ArithmeticError("The released source fails a binding or local-colour bound")
    boards = []
    ideal_dict = {word: complex(value) for word, value in zip(sector.states, ideal, strict=True)}
    finite_dict = {word: complex(value) for word, value in zip(sector.states, finite, strict=True)}
    for level in range(3):
        board = federation(level)
        boundary = connected_boundary(board)
        moments = []
        for name, wave in (("ideal_line", ideal_dict), ("actual_finite_filter", finite_dict)):
            applied = action(board, wave)
            energy = float(sum(np.conj(value)*applied[word] for word, value in wave.items()).real)
            variance = float(sum(abs(value-energy*wave.get(word, 0))**2 for word, value in applied.items()))
            moments.append({"source": name, "shifted_energy": energy, "residual_norm_squared": variance})
        boards.append({"level": level, "cells": board.cells, "ports": board.n,
                       "complete_three_load_dimension": board.n*(board.n-1)*(board.n-2)//6,
                       "board_sha256": sha256((CARRIER_HOME/"carriers"/f"FEDERATION_LEVEL{level}.npz").read_bytes()).hexdigest(),
                       "boundary": boundary, "source_moment_diagnostics": moments,
                       "complete_seam_energy_records": seam_energies(board, finite_dict)})
    return {"schema": "oph.interval_echo.released_colour_trimer.v1",
            "measure_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "spectral_frame_sha256": sha256((HERE / "spectral_frame.npz").read_bytes()).hexdigest(),
            "proof_sha256": sha256((HERE/"PROOF.md").read_bytes()).hexdigest(),
            "original_frames_receipt_sha256": sha256(FRAMES.read_bytes()).hexdigest(),
            "source": source_info, "boards": boards,
            "single_terminal_read_bound_uniform_in_exposure": {"ideal_cluster_residual_upper": str(sigma),
                "rational_ideal_connected_lower": str(floor), "finite_filter_connected_lower": prepared_floor,
                "original_start_joint_lower": source_info["finite_filter_admission_interval"][0]*prepared_floor,
                "scope": "For every chosen t, a fresh admitted release has this connected-read probability. No continuous monitoring or joint all-time path probability."},
            "all_board_local_colour_witness_bounds": witness_bounds,
            "scope": "Same sharp trimer, timing filter and frozen local frames after all contacts open. Ordinary integral loads; no quark identity, pure colour singlet or baryon assignment. No face flag before the local colour witness."}


def compare(actual, archived):
    if type(actual) is not type(archived):
        raise ArithmeticError("Receipt type differs")
    if isinstance(actual, dict):
        if actual.keys() != archived.keys():
            raise ArithmeticError("Receipt fields differ")
        for key in actual:
            compare(actual[key], archived[key])
    elif isinstance(actual, list):
        if len(actual) != len(archived):
            raise ArithmeticError("Receipt outcome partition differs")
        for left, right in zip(actual, archived, strict=True):
            compare(left, right)
    elif isinstance(actual, float):
        if not np.isfinite(archived) or abs(actual-archived) > TOLERANCE:
            raise ArithmeticError("Receipt probability or bound differs")
    elif actual != archived:
        raise ArithmeticError("Exact certificate differs")


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record", type=Path, help="explicit exclusive first-record creation")
    args = cli.parse_args()
    result = run()
    payload = json.dumps(result, indent=2, allow_nan=False)+"\n"
    if args.record:
        with args.record.open("x", encoding="utf-8") as destination:
            destination.write(payload)
    else:
        compare(result, json.loads(ARCHIVE.read_text()))
    print(result["single_terminal_read_bound_uniform_in_exposure"])
    print(result["all_board_local_colour_witness_bounds"])


if __name__ == "__main__":
    main()
