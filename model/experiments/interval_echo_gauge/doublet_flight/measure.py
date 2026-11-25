"""Two-load transmission of the intrinsic pair plane across an owner gateway."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
from model.experiments.count_intervals import clopper_pearson

from model.candidate_laws.interval_echo import currents
from model.candidate_laws.contrast_memory_spherical_v1.sectors import (
    PairSector,
    pair_index,
)
from model.executable_model import cell, federation
from model.executable_model.carrier import Carrier

HERE = Path(__file__).parent
SHOTS = 2**40
TIMES = (0, 8, 32)
# Simultaneous Clopper-Pearson coverage for every terminal seam in every bank.
ALPHA = 0.01 / (2 * 2 * 2 * len(TIMES) * 2 * 2 * 30)


def endpoint(port):
    car = cell()
    sector, stars, _, _, bands, _, _, _, h, _ = currents._service(car, True, 1.0)
    band = bands[0]
    opposite = currents._geometry(car)[0][port]
    bright = band @ (band.T @ stars[:, port])
    bright /= np.linalg.norm(bright)
    turn = next(rot for rot in currents._one_edge_turns(car) if rot[port] == port)
    neighbor = min(car.neighbours[port])
    plane = []
    for q in (neighbor, turn[neighbor]):
        v = band @ (band.T @ stars[:, q])
        v -= bright * (bright @ v)
        for old in plane:
            v -= old * (old @ v)
        plane.append(v / np.linalg.norm(v))
    frame = np.column_stack([bright, *plane])
    basis = np.column_stack(bands)
    triplet_current = h[port][:3, :3] + h[opposite][:3, :3]
    rates = np.linalg.eigvalsh(triplet_current)
    duration = 2 * np.pi / rates[-1]
    word = expm(-1j * duration * h[opposite]) @ expm(-1j * duration * h[port])
    tag = np.eye(sector.P) + basis @ (word - np.eye(6)) @ basis.T
    phase = np.exp(-1j * duration * rates[0])
    p = frame[:, 1:] @ frame[:, 1:].T
    np.testing.assert_allclose(tag, np.eye(sector.P) + (phase - 1) * p, atol=3e-13)
    r = currents._permutation(turn, sector)
    averaged = sum(np.linalg.matrix_power(r, j) for j in range(5)) / 5
    np.testing.assert_allclose(averaged @ frame, frame * [1, 0, 0], atol=1e-12)
    return dict(
        frame=frame,
        tag=tag,
        phase=phase,
        stars=stars[:, [port, neighbor]],
        turn=np.asarray(turn),
        neighbor=neighbor,
        duration=duration,
    )


def fixtures():
    full = federation(0)
    a, b = full.external[0]
    assert a // 12 == 0
    guarded = Carrier(
        "two owners and one bridge",
        24,
        [*cell().seams, *((x + 12, y + 12) for x, y in cell().seams), (a, 12 + b % 12)],
    )
    return [(guarded, 1), (full, b // 12)], a, b % 12


def turn_preserves(board, owner, turn):
    rotation = np.arange(board.n)
    rotation[owner * 12 : (owner + 1) * 12] = owner * 12 + turn
    return {tuple(sorted((rotation[i], rotation[j]))) for i, j in board.seams} == set(
        board.seams
    )


def intervals(counts):
    return clopper_pearson(counts, SHOTS, ALPHA / 2)


def trajectory(board, destination, home, away, kappa, geometry):
    sector = PairSector(board)
    small = PairSector(cell())
    ins = pair_index(board.n, small.aa, small.bb)
    outs = pair_index(board.n, 12 * destination + small.aa, 12 * destination + small.bb)
    h = sector.hamiltonian(kappa)
    starts = np.zeros((sector.P, 5), complex)
    starts[ins] = np.column_stack([home["frame"], home["stars"]])
    # The two symmetry columns and their phase-tagged variants follow by linearity.
    states = expm_multiply(-1j * h, starts, start=0, stop=32, num=5)
    results = []
    for time in TIMES:
        local = states[time // 8, outs]
        flying_frame, flying_stars = local[:, :3], local[:, 3:]
        plane = home["frame"][:, 1:]
        tag_update = (
            (home["phase"] - 1) * flying_frame[:, 1:] @ (plane.T @ home["stars"])
        )
        transfer = away["frame"].T @ flying_frame
        plane_sq = np.linalg.svd(transfer[1:, 1:], compute_uv=False) ** 2
        whole_plane_sq = np.linalg.svd(flying_frame[:, 1:], compute_uv=False) ** 2
        rows = []
        for origin in (0, 1):
            for receiver in (0, 1):
                probabilities, counts = [], []
                for source in (0, 1):
                    v = flying_stars[:, origin] + source * tag_update[:, origin]
                    if receiver:
                        v = away["tag"] @ v
                    p = abs(v[small.seam_states]) ** 2
                    total = np.r_[p, 1 - np.sum(p)]
                    rng = np.random.default_rng(
                        [683214, geometry, int(kappa), time, origin, receiver, source]
                    )
                    probabilities.append(total)
                    counts.append(rng.multinomial(SHOTS, total))
                low0, high0 = intervals(counts[0][:30])
                low1, high1 = intervals(counts[1][:30])
                distinguished = (low0 > high1) | (low1 > high0)
                rows.append(
                    dict(
                        source_kind=["on_axis", "cross_axis"][origin],
                        receiver_tag=bool(receiver),
                        law=[p.tolist() for p in probabilities],
                        counts=[c.tolist() for c in counts],
                        resolved_seams=np.flatnonzero(distinguished).tolist(),
                        total_variation=float(
                            np.abs(probabilities[1] - probabilities[0]).sum() / 2
                        ),
                    )
                )
        results.append(
            dict(
                exposure=time,
                plane_to_plane_squared=plane_sq.tolist(),
                plane_to_whole_receiver_squared=whole_plane_sq.tolist(),
                bright_to_bright_squared=float(abs(transfer[0, 0]) ** 2),
                bright_to_whole_receiver_squared=float(
                    np.linalg.norm(flying_frame[:, 0]) ** 2
                ),
                mixed_triplet_endpoint_weights=(abs(transfer) ** 2).tolist(),
                phase_records=rows,
            )
        )
    return results


def portable(x):
    if isinstance(x, dict):
        return {k: portable(v) for k, v in x.items()}
    if isinstance(x, list):
        return [portable(v) for v in x]
    if isinstance(x, float):
        return float(f"{x:.10g}") if abs(x) > 1e-25 else 0.0
    return x


def run():
    setup, p, q = fixtures()
    home, away = endpoint(p), endpoint(q)
    rows = []
    for geometry, (board, destination) in enumerate(setup):
        for kappa in (currents.KAPPA, 0):
            rows.append(
                dict(
                    geometry=board.name,
                    destination_owner=destination,
                    kappa=kappa,
                    independent_turns=[
                        turn_preserves(board, 0, home["turn"]),
                        turn_preserves(board, destination, away["turn"]),
                    ],
                    samples=trajectory(board, destination, home, away, kappa, geometry),
                )
            )
    return portable(
        dict(
            protocol_sha256=hashlib.sha256(
                (HERE / "PROTOCOL.md").read_bytes()
            ).hexdigest(),
            exposures=list(TIMES),
            trials_per_bank=SHOTS,
            source_axis=p,
            receiver_axis=q,
            source_neighbor=home["neighbor"],
            tag_duration=float(home["duration"]),
            receiver_seams=cell().seams,
            arms=rows,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    destination = parser.parse_args().output
    result = run()
    with destination.open("x") as stream:
        json.dump(result, stream, indent=1, allow_nan=False)
        stream.write("\n")
