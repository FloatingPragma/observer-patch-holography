"""Counter records for the six winding/arm preparations of a released neutral pair."""

import argparse
import hashlib
import json
from functools import cache
from itertools import permutations, product
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm
from scipy.sparse.linalg import expm_multiply

from model.experiments.interval_echo_angular.charged_orientation.live_trimer.measure import (
    LAW,
    blocks,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin import (
    apparatus as arms,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge import (
    instruments,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release import (
    flight,
    geometry,
)

HERE = Path(__file__).resolve().parent
SIZE = 2**52
SEED = 26112610
TIMES = (256, 512, 2048)
OWNERS = (12, 17)
BASES = ("character", "arm")
WORDS = tuple(permutations(range(3)))
PANELS = ("sharp", "left", "right")


def bank(key, probabilities, size=SIZE):
    """Independent, order-invariant counting cohorts. Preserve small bins explicitly."""
    p = np.asarray(probabilities, dtype=float)
    if p.min() < -2e-15 or abs(p.sum() - 1) > 1e-11:
        raise ValueError("Incomplete or negative outcome probabilities")
    digest = hashlib.sha256(json.dumps([SEED, *key]).encode("ascii")).digest()
    rng = np.random.default_rng(int.from_bytes(digest[:16], "big"))
    order = np.argsort(p, kind="stable")
    ordered = np.maximum(p[order], 0)
    draws = rng.multinomial(size, ordered / ordered.sum())
    return draws[np.argsort(order)].tolist()


@cache
def source_weights():
    f = arms.fixture()
    selected = [f["marked"][3:].index(e) for e in geometry.layout()[1]]
    code = instruments.material("hub")
    source = np.column_stack([code[b] for b in BASES])
    h = -32j * arms.native_graph()
    initial = np.zeros((h.shape[0], 6), complex)
    initial[:8] = source
    hub = f["ports"].index(f["tip"])
    field = expm_multiply(h, initial, traceA=h.diagonal().sum())
    received = field[64 * hub : 64 * (hub + 1)].reshape(8, 8, 6)
    sharp = (abs(received) ** 2).sum(axis=1)
    masks = np.arange(8)
    labels = sum(((masks >> j) & 1) << i for i, j in enumerate(selected))
    return np.array([sharp[labels == m].sum(axis=0) for m in range(4)])


@cache
def receiver(owner):
    ports = tuple(p for p in geometry.layout()[0] if p // 12 == owner)
    if len(ports) != 3 or owner not in OWNERS:
        raise ValueError("Receiver needs the three enrolled ring ports")
    return ports


def contact_word_matrix(which):
    h = np.zeros((len(WORDS), len(WORDS)))
    a, b = which, which + 1
    blank = blocks()["blank"] / 8
    for j, w in enumerate(WORDS):
        col = LAW.SRC.index((w[a], w[b]))
        for r in np.flatnonzero(blank[:, col]):
            new = list(w)
            new[a], new[b] = LAW.SRC[r]
            h[WORDS.index(tuple(new)), j] = blank[r, col]
    return h


@cache
def analyzers():
    return dict(
        zip(
            PANELS,
            [
                np.eye(6),
                *(expm(-0.5j * np.pi * contact_word_matrix(i)) for i in range(2)),
            ],
            strict=True,
        )
    )


def enrol_matrix(states, strength, cut):
    """Accumulate host-subtracted registered contact blocks in this source sector."""
    loc = {state: i for i, state in enumerate(states)}
    h8 = np.zeros((len(states), len(states)), dtype=int)
    for a, b, bit in flight.edges(True):
        if cut and (a // 12 == 13) != (b // 12 == 13):
            continue
        recorded = bit is not None
        raw = blocks(strength)[(1, 1) if recorded else "blank"]
        for j, (vac, exc, mask) in enumerate(states):
            if not {a, b}.intersection((vac, exc)):
                continue
            sites = sorted((a, b))
            pair = tuple(0 if s == vac else 2 if s == exc else 1 for s in sites)
            coin = (((mask >> bit) & 1) ^ int(a > b)) if recorded else 0
            col = 9 * coin + LAW.SRC.index(pair)
            written_host = 9 * coin + LAW.SRC.index((1, 1))
            h8[j, j] -= raw[written_host, written_host]
            for r in np.flatnonzero(raw[:, col]):
                if r == col:
                    i = j
                else:
                    if LAW.SRC[r % 9] != pair[::-1] or (recorded and r // 9 == coin):
                        raise ValueError("Contact escaped registered swap sector")
                    swap = {a: b, b: a}
                    target = (
                        swap.get(vac, vac),
                        swap.get(exc, exc),
                        mask ^ (1 << bit) if recorded else mask,
                    )
                    i = loc[target]
                h8[i, j] += raw[r, col]
    if not np.array_equal(h8, h8.T):
        raise ValueError("Free flight generator lost reciprocity")
    return h8 / 8


def sector_fields(strength, cut, times=TIMES):
    """Each local vector is unnormalized. Outside the owner all other records match."""
    by_owner = {owner: [] for owner in OWNERS}
    for mask in range(4):
        states, native = flight.chart((160, 167, mask), True)
        h = enrol_matrix(states, strength, cut)
        if strength == 8 and not cut:
            np.testing.assert_array_equal(h, native)
        vals, vecs = eigh(h, driver="evd")
        evolved = (np.exp(-1j * np.outer(times, vals)) * vecs[0]) @ vecs.T
        index = {s: i for i, s in enumerate(states)}
        for owner in OWNERS:
            ports = receiver(owner)
            ids = [index[ports[w.index(0)], ports[w.index(2)], mask] for w in WORDS]
            if any(m != mask for v, e, m in states if v in ports and e in ports):
                raise ValueError("Source sector cannot be traced independently")
            field = evolved[:, ids]
            if cut:
                np.testing.assert_allclose(field, 0, atol=2e-11, rtol=0)
                field = np.zeros_like(field)  # graph-disconnected components
            by_owner[owner].append(field)
    return {k: np.array(f) for k, f in by_owner.items()}


def density(weights, columns):
    return np.einsum("m,mi,mj->ij", weights, columns, columns.conj())


def acquire():
    weights = source_weights()
    records, ideal, standards = [], [], []
    for panel, u in analyzers().items():
        for j in range(6):
            standards.append(
                dict(
                    panel=panel,
                    word=j,
                    counts=bank(["standard", panel, j], abs(u[:, j]) ** 2),
                )
            )
    for strength, cut in ((8, False), (4, False), (8, True)):
        fields = sector_fields(strength, cut)
        for owner, (i, t) in product(OWNERS, enumerate(TIMES)):
            dms = []
            for a in range(6):
                basis, j = BASES[a // 3], a % 3
                key = [strength, cut, owner, t, basis, j]
                rho = density(weights[:, a], fields[owner][:, i, :])
                prior = float(weights[:, a].sum())
                dms.append(rho)
                for panel, u in analyzers().items():
                    clicks = (u @ rho @ u.conj().T).diagonal().real
                    probs = [*clicks, prior - clicks.sum(), 1 / 3 - prior, 2 / 3]
                    records.append(
                        dict(key=[*key, panel], counts=bank([*key, panel], probs))
                    )
            for jj, family in enumerate(BASES):
                for a, b in ((0, 1), (0, 2), (1, 2)):
                    x, y = jj * 3 + a, jj * 3 + b
                    rx, ry = dms[x], dms[y]
                    px, py = rx.trace().real, ry.trace().real
                    # With equal prior, the optimum conditional success given the pair
                    # is (1+D)/2. For the complete protocol every refused branch counts.
                    trace_norm = np.abs(np.linalg.eigvalsh(rx - ry)).sum()
                    qa, qb = weights[:, x].sum(), weights[:, y].sum()
                    d_all = (trace_norm + abs((qa - px) - (qb - py)) + abs(qa - qb)) / 2
                    d_arr = None
                    if px > 0 and py > 0:
                        d_arr = float(
                            np.abs(np.linalg.eigvalsh(rx / px - ry / py)).sum() / 2
                        )
                    ideal.append(
                        dict(
                            key=[strength, cut, owner, t, family, a, b],
                            pair_given_both_classes_arrive=d_arr,
                            complete_trial_distance=float(d_all),
                        )
                    )
    return dict(
        seed=SEED,
        starts_per_bank=SIZE,
        receivers={k: receiver(k) for k in OWNERS},
        sharp_words=WORDS,
        source_bases=BASES,
        final_appointments=TIMES,
        count_order=[
            "six ordered pair words",
            "neutral pair outside receiver",
            "charged vacancy failed to reach converter",
            "source aperture refused",
        ],
        retained_mask_weights=weights.tolist(),
        local_word_standards=standards,
        banks=records,
        ideal_receiver_bounds=ideal,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    args.output.write_text(
        json.dumps(acquire(), separators=(",", ":"), allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()
