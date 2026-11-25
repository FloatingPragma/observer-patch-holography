"""Lowest source-visible and lowest extended bands at a fixed charged emitter."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.stats import beta as beta_law

from .. import measure as source
from ..apparatus import LAW, blocks, native_graph
from ..reciprocal import channel
from .guides import detour
from .measure import clock

HERE = Path(__file__).resolve().parent
SHOTS = 2**46
WITNESS_SHOTS = 2**51
DELAYS = 2**36
DT = 1 / 32
SEED = 903825
TAIL = 0.01 / 1024


def ci(k, n):
    return np.array(
        [
            beta_law.ppf(TAIL, k, n - k + 1) if k else 0.0,
            beta_law.isf(TAIL, k + 1, n - k) if k < n else 1.0,
        ]
    )


def draw(probs, rng, n=SHOTS):
    assert abs(sum(probs) - 1) < 2e-12 and min(probs) >= -1e-14
    p = np.maximum(probs, 0.0)
    return rng.multinomial(n, p / p.sum()).tolist()


def positions(z, xx, rng):
    p = np.sum(abs(z.reshape(-1, 64)) ** 2, axis=1)
    mid = (3 * xx >= xx.max()) & (3 * xx <= 2 * xx.max())
    fractions = np.array([p[mid].sum(), p[1], p[~mid].sum() - p[1]]) / 3
    emitted = float(sum(fractions))
    counts = draw([*fractions, 1 - emitted], rng)
    adm = sum(counts[:3])
    return dict(
        original_admission=emitted,
        mid_and_far_conditional=(fractions[:2] / emitted).tolist(),
        counts_mid_far_other_refusal=counts,
        classified_as_visible=bool(ci(adm, SHOTS)[0] * SHOTS >= 100),
        mid99=ci(counts[0], adm).tolist() if adm else [0.0, 1.0],
        remote_original99=ci(counts[1], SHOTS).tolist(),
    )


def witness(K, rng):
    success, minor, J = source.pt(K)
    axes = []
    boxes = []
    for op in [source.X, source.Y, source.Z]:
        corr = np.trace(J @ np.kron(op, op)).real
        pair = np.array([success + corr, success - corr]) / 18
        counts = draw([*pair, 1 - sum(pair)], rng, WITNESS_SHOTS)
        axes.append(counts)
        boxes.append([ci(c, WITNESS_SHOTS) for c in counts[:2]])
    a, b, c = boxes
    return dict(
        yield_per_original=success / 9,
        conditional_pt_minimum=minor / success,
        counts_same_different_all_refusals=axes,
        witness99=[
            (2 * a[1][0] + b[0][0] - b[1][1] - c[0][1] + c[1][0]) / 4,
            (2 * a[1][1] + b[0][1] - b[1][0] - c[0][0] + c[1][1]) / 4,
        ],
    )


def length_run(count, rng):
    _, ports, legs, h, actions, q, T = detour(count)
    ev, v = eigh(T, driver="evd")
    L = len(legs[0]) + 1
    xx = np.zeros(len(ports), dtype=int)
    xx[1] = L
    for leg in legs:
        for i, p in enumerate(leg, 1):
            xx[ports.index(p)] = i
    groups = [[ports.index(leg[j]) for leg in legs] for j in range(L - 1)]
    # Derive the bulk recurrence directly on every memory configuration away
    # from both local triangles, including every nonzero physical matrix entry.
    for j in range(3, L - 4):
        for a, b, c in zip(groups[j - 1], groups[j], groups[j + 1], strict=True):
            for head in range(64):
                r = h.getrow(64 * b + head)
                assert dict(zip(r.indices, r.data, strict=True)) == {
                    64 * a + head: -0.5,
                    64 * b + head: 17.0,
                    64 * c + head: -0.5,
                }
    beta = source.BASE[:, 0] @ (q[source.ONE] @ v)
    modes = q @ v

    def band(omega):
        return beta * clock(ev, dict(center=omega, tick=DT, number_delays=DELAYS))

    first = ev.searchsorted(16.0)
    second = ev.searchsorted(17 - np.cos(1.5 * np.pi / L))
    # Aggregate groups unresolved at the clock's first sinc zero. This also
    # removes the arbitrary choice of eigenvectors within a nearly degenerate pair.
    spacing = 2 * np.pi / (DT * DELAYS)
    breaks = np.r_[0, np.flatnonzero(np.diff(ev[:first]) > spacing) + 1, first]
    below = []
    for i, j in zip(breaks[:-1], breaks[1:], strict=True):
        center = float(np.mean(ev[i:j]))
        field = modes @ band(center)
        projection = beta[i:j]
        below.append(
            dict(
                center=center,
                multiplicity=int(j - i),
                source_weight=np.dot(projection, projection),
                **positions(field, xx, rng),
            )
        )
    cycle = actions[source.CY]
    src = np.zeros(len(xx) * 64)
    src[source.ONE] = source.BASE[:, 0]
    sign = np.sign(
        source.BASE[:, 1] @ ((src[cycle] - src[cycle[cycle]])[source.ONE] / np.sqrt(3))
    )
    assert abs(sign) == 1
    low = []
    packets = []
    for k in [first, second]:
        coeff = band(ev[k])
        packets.append(coeff)
        wave = modes @ coeff
        doublet = np.stack(
            [wave, sign * (wave[cycle] - wave[cycle[cycle]]) / np.sqrt(3)], axis=1
        )
        K = channel(doublet, 0)[0]
        low.append(
            dict(
                energy=ev[k],
                input_weight=beta[k] ** 2,
                filtered_mode_weight=(abs(coeff[k]) / np.linalg.norm(coeff)) ** 2,
                wave_number=np.arccos(17 - ev[k]),
                **positions(wave, xx, rng),
                winding=witness(K, rng),
            )
        )
    # Coherently erase the choice of the two equally weighted clock centers.
    # Then let the unmodified guide run freely and compare fresh local banks.
    coeff = (packets[0] + packets[1]) / 2
    admitted = np.linalg.norm(coeff) ** 2 / 3
    gap = ev[second] - ev[first]
    cut = L // 2
    left, right = groups[cut - 1], groups[cut]
    quarter = 2 * np.pi / gap / 4
    trace = []
    for n in range(5):
        t = n * quarter
        wave = (modes @ (coeff * np.exp(-1j * (ev - ev[first]) * t))).reshape(-1, 64)
        wave /= np.linalg.norm(wave)
        plus = 0.5 * np.linalg.norm(1j * wave[left] + wave[right]) ** 2
        minus = 0.5 * np.linalg.norm(wave[left] + 1j * wave[right]) ** 2
        # H_lr=-1/2, so current = (right-left probability)/2 after the physical half-swap.
        current = (plus - minus) / 2
        vector = wave.ravel()
        Rmask = np.repeat(xx > cut, 64)
        actual_velocity = 2 * np.vdot(vector[Rmask], (h @ vector)[Rmask]).imag
        assert abs(actual_velocity - current) < 1e-12
        stream = draw(
            [
                admitted * plus,
                admitted * minus,
                admitted * (1 - plus - minus),
                1 - admitted,
            ],
            rng,
        )
        nread = sum(stream[:3])
        hi, lo = (ci(c, nread) for c in stream[:2])
        far = np.linalg.norm(wave[xx > cut]) ** 2
        pos = draw([admitted * far, admitted * (1 - far), 1 - admitted], rng)
        trace.append(
            dict(
                time=t,
                current=float(current),
                current99=((hi - lo[::-1]) / 2).tolist(),
                current_counts_forward_backward_else_refused=stream,
                far_half=float(far),
                far_half99=ci(pos[0], sum(pos[:2])).tolist(),
                counts_far_near_refused=pos,
            )
        )
    # Incoherently mixing the two same narrow filters removes the beat while
    # retaining within-window finite-clock tails. Perform that control at the
    # same measurement times on the same freely evolving registered generator.
    mixed = []
    for row in trace:
        weights = np.zeros(len(xx))
        for a in packets:
            w = (modes @ (a * np.exp(-1j * (ev - ev[first]) * row["time"]))).reshape(
                -1, 64
            )
            weights += np.sum(abs(w) ** 2, axis=1) / 12
        p = float(sum(weights))
        far = float(sum(weights[xx > cut]))
        counts = draw([far, p - far, 1 - p], rng)
        mixed.append(
            dict(
                far_half=far / p,
                counts_far_near_refused=counts,
                far_half99=ci(counts[0], sum(counts[:2])).tolist(),
            )
        )
    return dict(
        length=L,
        physical_cuts=[[ports[a], ports[b]] for a, b in zip(left, right, strict=True)],
        below=below,
        lowest_and_next_harmonic=low,
        scaled_first_excess=L**2 * (ev[first] - 16),
        coherent_source_admission=admitted,
        free_beat=trace,
        incoherent_control=mixed,
    )


def run():
    # Both registered charge preparations reduce to the same matrix, including
    # the detours because their blank exchange term is charge-complement invariant.
    assert (native_graph(1) != native_graph(-1)).nnz == 0
    raw = blocks()["blank"]
    for q in [-1, 1]:
        pair = [LAW.SRC.index((1 + q, 1)), LAW.SRC.index((1, 1 + q))]
        np.testing.assert_array_equal(
            raw[np.ix_(pair, pair)] / 8, [[8.5, -0.5], [-0.5, 8.5]]
        )
    # Matched changed-law bulk control. On each actual, uniform interior
    # site these two checked blank bonds comprise every incident nonzero term.
    changed = blocks(4)["blank"] / 8
    rails = [LAW.SRC.index((0, 1)), LAW.SRC.index((1, 0))]
    local = changed[np.ix_(rails, rails)]
    np.testing.assert_array_equal(local, [[4.5, -0.5], [-0.5, 4.5]])
    control_bottom = 2 * (local[0, 0] + local[0, 1])
    rng = np.random.default_rng(SEED)
    return dict(
        changed_kappa4_homogeneous_bottom=control_bottom,
        shots_per_position_or_current=SHOTS,
        shots_per_witness_axis=WITNESS_SHOTS,
        seed=SEED,
        clock_delays=DELAYS,
        clock_tick=DT,
        tail=TAIL,
        lengths=[length_run(n, rng) for n in range(3)],
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    result = run()
    if ap.parse_args().write:
        (HERE / "threshold.json").write_text(json.dumps(result, indent=2) + "\n")
    else:
        for row in result["lengths"]:
            print(
                row["length"],
                row["lowest_and_next_harmonic"][0]["energy"],
                row["free_beat"][0]["far_half99"],
                row["free_beat"][2]["far_half99"],
            )
