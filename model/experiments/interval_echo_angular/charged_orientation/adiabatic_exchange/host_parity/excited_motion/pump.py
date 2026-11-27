"""Two ordinary vacancies in an excited band of a native three-spoke guide."""

from functools import lru_cache
from itertools import combinations

import numpy as np
from scipy.linalg import expm

from ..experiment import LAW
from ..passive import commission, local_certificate, sharp_route

SPOKES = (1, 2, 3)
PAIRS = ((1, 2), (2, 3), (1, 3))
PHI = 2 * np.pi * np.arange(3) / 3


def aperture(theta):
    return 2 * (1 + np.cos(theta - PHI))


@lru_cache(None)
def stocks(count):
    return tuple(
        tuple(int(p not in holes) for p in range(4))
        for holes in combinations(range(4), count)
    )


@lru_cache(None)
def source_geometry(count):
    words = stocks(count)
    lookup = {w: i for i, w in enumerate(words)}
    h = np.zeros((len(words), len(words)))
    edge = LAW.seam_term(8)
    for j, w in enumerate(words):
        for b in SPOKES:
            pair = (w[0], w[b])
            idx = LAW.SRC.index(pair)
            h[j, j] += edge[idx, idx]
            if pair[0] != pair[1]:
                moved = list(w)
                moved[0], moved[b] = moved[b], moved[0]
                h[lookup[tuple(moved)], j] += edge[LAW.SRC.index(pair[::-1]), idx]
    h -= np.eye(len(h)) * h[0, 0]  # common phase in both observed arms
    pairs = np.array([[(1 - w[a]) * (1 - w[b]) for a, b in PAIRS] for w in words])
    return h, pairs, np.array([w[0] for w in words])


@lru_cache(None)
def native_flight(count, dt):
    h, _, _ = source_geometry(count)
    return expm(-1j * dt * h)


def clocks(count, detuning, theta, dt):
    h, pairs, central = source_geometry(count)
    d = detuning * central - pairs @ aperture(theta)
    ph = np.exp(-0.5j * dt * d)
    step = ph[:, None] * native_flight(count, dt) * ph
    # Sine is injective on the occupied half circle.
    _, q = np.linalg.eigh(-step.imag)
    phase = np.diag(q.T @ step @ q)
    if np.linalg.norm(step @ q - q * phase) > 1e-8:
        raise ArithmeticError("Native step has unresolved eigenspaces")
    energy = -np.angle(phase) / dt
    assert np.max(abs(energy)) * dt < np.pi / 2
    order = np.argsort(energy)
    return energy[order], q[:, order]


def route_certificate():
    local_certificate()
    native = LAW.seam_term(8)
    pair_count = 0
    for angle in (0.0, 0.073, 0.349, 1.0):
        # Native two-port controlled dwell, including the precise idle phases.
        u = expm(-0.5j * angle * native / 17)
        z = np.diag([(-1) ** (a * a) for a, b in LAW.SRC])
        private = np.diag([np.exp(-0.5j * angle * (a * a + b * b)) for a, b in LAW.SRC])
        full = private @ z @ u @ z @ u
        for a, b in ((0, 0), (0, 1), (1, 0), (1, 1)):
            v = np.eye(9)[:, LAW.SRC.index((a, b))]
            np.testing.assert_allclose(
                full @ v, v * np.exp(-1j * angle * (1 - (1 - a) * (1 - b))), atol=1e-12
            )
            pair_count += 1
    # Same physical Y as the previously commissioned one; no flight apertures
    # beyond these three blank spokes and their four private references act.
    old = commission(0, 3)
    ports = old["ports"]
    stock = (1, 0, 0, 1)
    start = stock
    moves = []
    for a, b in old["quiet_contacts"]:
        chosen = [ports.index(p) for p in (a, b) if p in ports]
        for leaf in sharp_route(stock, lambda w: all(w[k] == 1 for k in chosen)):
            nxt = list(stock)
            nxt[0], nxt[leaf] = nxt[leaf], nxt[0]
            stock = tuple(nxt)
            moves.append(("swap", ports[0], ports[leaf]))
        assert all(stock[k] == 1 for k in chosen)
        moves.append(("write_equal", a, b))
    for leaf in sharp_route(stock, lambda w: w == start):
        nxt = list(stock)
        nxt[0], nxt[leaf] = nxt[leaf], nxt[0]
        stock = tuple(nxt)
        moves.append(("swap", ports[0], ports[leaf]))
    assert stock == start
    return dict(
        pair_gate_columns=pair_count, ports=ports, sharp_stock=start, operations=moves
    )


def calendar(count, detuning, segments, dt, guarded=False):
    h, pairs, central = source_geometry(count)
    angles = (np.arange(segments) + 0.5) * 2 * np.pi / segments
    apertures = 2 * (1 + np.cos(angles[:, None] - PHI))
    d = detuning * central[None, :] - apertures @ pairs.T
    ph = np.exp(-0.5j * dt * d)
    u = ph[:, :, None] * native_flight(count, dt)[None, :, :] * ph[:, None, :]
    if guarded:
        words = stocks(count)
        for leg, b in enumerate((2, 3, 1)):
            section = slice(leg * segments // 3, (leg + 1) * segments // 3)
            restricted = h.copy()
            paid = []
            for i, w in enumerate(words):
                charge = 8.5 * (w[0] - w[b]) ** 2
                paid.append(charge)
                restricted[i, i] -= charge
                if w[0] != w[b]:
                    swapped = list(w)
                    swapped[0], swapped[b] = swapped[b], swapped[0]
                    restricted[words.index(tuple(swapped)), i] = 0
            if guarded == "parked":
                # Compiled unconditional unary controls never operate on b.
                # Retain their true action even in events where b is filled.
                masks = pairs.copy()
                for link, (a, c) in enumerate(PAIRS):
                    if b in (a, c):
                        movable = c if a == b else a
                        masks[:, link] = [1 - w[movable] for w in words]
                controlled = (
                    detuning * central - apertures[section] @ masks.T + 8.5 * central
                )
                diag = np.exp(-0.5j * dt * controlled)
            else:
                diag = ph[section] * np.exp(-0.5j * dt * np.array(paid))[None, :]
            u[section] = (
                diag[:, :, None]
                * expm(-1j * dt * restricted)[None, :, :]
                * diag[:, None, :]
            )
            occupied = np.array([w[b] for w in words])
            if (
                np.max(abs(u[section] * (occupied[:, None] != occupied[None, :])))
                > 1e-13
            ):
                raise ArithmeticError("The guarded vacancy can leave its spoke")
    _, q = np.linalg.eigh(-u.imag)
    diag = np.einsum("tjk,tjk->tk", q, u @ q)
    if np.max(abs(u @ q - q * diag[:, None, :])) > 1e-8:
        raise ArithmeticError("Calendar phase collision")
    e = -np.angle(diag) / dt
    assert np.max(abs(e)) * dt < np.pi / 2
    order = np.argsort(e, axis=1)
    return np.take_along_axis(e, order, axis=1), np.take_along_axis(
        q, order[:, None, :], axis=2
    )


def index_pair(pair):
    return stocks(2).index(tuple(int(p not in pair) for p in range(4)))


def launch_law(row, accept):
    if np.isscalar(row):
        row = [1 - row, row]
    p = np.r_[1 - accept, np.array(row) * accept]
    if min(p) < -1e-8 or abs(sum(p) - 1) > 1e-7:
        raise ArithmeticError("Source/exit normalization")
    return p.clip(0) / p.clip(0).sum()


def run_case(detuning, dt, segments, turns, shots=2**22, guarded=False):
    count = 2
    initial = index_pair(PAIRS[0])
    source = np.eye(6)[:, initial]
    eigenindex = 3 if detuning > 0 else 0
    e0, q0 = clocks(count, detuning, 0.0, dt)
    branch = q0[:, eigenindex]
    admitted = float(branch[initial] ** 2)
    # Project by a period-8192 eigenphase clock with the nearest-bin effect,
    # retaining the probabilities and amplitudes of all six original lines.
    labels = 2**16
    tick = int(np.round(e0[eigenindex] * dt * labels / (2 * np.pi)))
    phi = (tick / labels - e0 * dt / (2 * np.pi) + 0.5) % 1 - 0.5
    response = (
        np.exp(1j * np.pi * (labels - 1) * phi) * np.sinc(labels * phi) / np.sinc(phi)
    )
    retained = q0 @ (response * (q0.T @ source))
    hit = float(np.linalg.norm(retained) ** 2)
    retained /= np.sqrt(hit)
    states = np.column_stack((source, retained))
    # Accumulate each full turn and each reflection-matched out-and-back path.
    energies, bases = calendar(count, detuning, segments, dt, guarded)
    addresses = [index_pair(x) for x in PAIRS]
    intervals = {segments * j // 12 for j in range(13)}
    overlap = 1.0
    last = branch
    gap = float("inf")
    for energy, q in zip(energies, bases):
        current = q[:, eigenindex]
        overlap *= np.sign(last @ current)
        last = current
        gap = min(
            gap, float(min(abs(np.delete(energy, eigenindex) - energy[eigenindex])))
        )
    overlap *= np.sign(last @ branch)
    records = []
    rng = np.random.default_rng(
        [
            9102707,
            detuning + 16,
            round(1 / dt),
            segments,
            2 if guarded == "parked" else int(guarded),
        ]
    )
    for duration in turns:
        m = duration / (segments * dt)
        if m != int(m):
            raise ValueError("Noninteger native microsteps")
        f = np.eye(6, dtype=complex)
        r = f.copy()
        # Separate sharp-address cohorts make sequential QND boundary reads.
        kept = states.copy()
        kept[[(w[2] != 0) for w in stocks(2)]] = 0
        pass_prob = [np.linalg.norm(kept, axis=0) ** 2]
        frames = [abs(states.T) ** 2]
        ref_frames = [frames[0]]
        all_steps = (
            bases * np.exp(-1j * duration / segments * energies)[:, None, :]
        ) @ bases.transpose(0, 2, 1)
        for j, u in enumerate(all_steps):
            f = u @ f
            if guarded == "parked":
                kept = u @ kept
                if (j + 1) % (segments // 3) == 0:
                    reached = (j + 1) // (segments // 3)
                    required = (
                        (3,) if reached == 1 else (1,) if reached == 2 else (1, 2)
                    )
                    kept[[any(w[p] != 0 for p in required) for w in stocks(2)]] = 0
                    pass_prob.append(np.linalg.norm(kept, axis=0) ** 2)
            mirrored = j if j < segments // 2 else segments - 1 - j
            ur = all_steps[mirrored]
            r = ur @ r
            if j + 1 in intervals:
                frames.append(abs((f @ states).T) ** 2)
                ref_frames.append(abs((r @ states).T) ** 2)
        row = dict(duration=duration)
        if guarded == "parked":
            row["surviving_sequential_reads"] = np.array(pass_prob).tolist()
            # Output 0 is source refusal. Then first failed boundary (0,1,2,3)
            # and the final success. Probabilities are nested in a single cohort.
            row["joint_route_counts"] = []
            for i, accept in enumerate((1, hit)):
                surv = np.r_[1, np.array(pass_prob)[:, i], 0]
                row["joint_route_counts"].append(
                    rng.multinomial(shots, launch_law(-np.diff(surv), accept)).tolist()
                )
        for name, a, b in [("once", f, r), ("twice", f @ f, r @ r)]:
            arm0, arm1 = a @ states, b @ states
            bright = np.linalg.norm(arm0 + arm1, axis=0) ** 2 / 4
            dark = np.linalg.norm(arm0 - arm1, axis=0) ** 2 / 4
            np.testing.assert_allclose(bright + dark, 1, atol=5e-9)
            row[name] = dict(
                dark=dark.tolist(),
                joint_launch_counts=[
                    rng.multinomial(
                        shots, launch_law(dark[i], 1 if i == 0 else hit)
                    ).tolist()
                    for i in range(2)
                ],
                signed_overlaps=[
                    float(np.vdot(arm1[:, i], arm0[:, i]).real) for i in range(2)
                ],
            )
        row["forward_word_frames"] = np.array(frames).tolist()
        row["ref_word_frames"] = np.array(ref_frames).tolist()
        for name, bank in [("forward", frames), ("reference", ref_frames)]:
            row[name + "_frame_counts"] = [
                [
                    rng.multinomial(shots, launch_law(p, 1 if i == 0 else hit)).tolist()
                    for i, p in enumerate(frame)
                ]
                for frame in bank
            ]
        row["min_pinned"] = [
            min(frames[4 * j][i, addresses[j % 3]] for j in range(4)) for i in range(2)
        ]
        row["min_crossing_shared_two_pairs"] = [
            min(
                frames[4 * j + 2][i, addresses[j]]
                + frames[4 * j + 2][i, addresses[(j + 1) % 3]]
                for j in range(3)
            )
            for i in range(2)
        ]
        records.append(row)
    return dict(
        detuning=detuning,
        microstep=dt,
        appointments=segments,
        shots=shots,
        guarded=guarded,
        source_energy=e0[eigenindex],
        source_admission=admitted,
        finite_bin=tick,
        finite_acceptance=hit,
        finite_purity=abs(branch @ retained) ** 2,
        min_sampled_gap=gap,
        polygon_sign=int(overlap),
        records=records,
    )
