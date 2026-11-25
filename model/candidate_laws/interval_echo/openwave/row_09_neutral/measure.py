"""Release the two wall daughters and measure energy and material charge on the same wing."""

import argparse
import json
from itertools import combinations, product
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix, eye
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation
from model.candidate_laws.interval_echo.openwave.row_08_annihilation import (
    experiment,
    walls,
    wall_flight,
)
from model.candidate_laws.interval_echo.openwave.row_08_annihilation.daughter_spectrum.measure import (
    check,
    compare,
    reference,
)
from model.experiments.count_intervals import clopper_pearson

HERE = Path(__file__).resolve().parent
N, SEED = 2**30, 9061125
COUNTER_ARMS = ((8, False), (8, True), (0, False))
SETTINGS = ("keep", "detach", "release")
CLOCK_BINS = (2**15, 2**28)
DT = 1 / 128
TAIL = 0.01 / 6000


def register_words():
    rows = []
    for filled in combinations(range(6), 3):
        word = np.zeros(6, int)
        word[list(filled)] = 2
        rows.append(word)
    return np.array(rows)


def full_wing(k, sites):
    words = register_words()
    index = {tuple(w): i for i, w in enumerate(words)}
    h = np.diag(k * ((words**2) @ (6 - federation(0).degree[sites])).astype(float))
    local = experiment.term(k)
    for i, w in enumerate(words):
        for a in range(5):
            pair = tuple(w[a : a + 2])
            j = experiment.LAW.SRC.index(pair)
            h[i, i] += local[j, j]
            if pair[0] != pair[1]:
                v = w.copy()
                v[a : a + 2] = pair[::-1]
                h[index[tuple(v)], i] += local[experiment.LAW.SRC.index(pair[::-1]), j]
    return h, words


def wing(k, side, boundary, cut, setting):
    sites = wall_flight.extended_geometry()[2][side]
    guarded_words = np.array([walls.loads(bits) for bits in product((0, 1), repeat=3)])
    if setting == "release":
        h, words = full_wing(k, sites)
    else:
        h = walls.wing_ref(k, boundary, cut or setting == "detach", 3)[0]
        h += np.diag(k * ((guarded_words**2) @ (6 - federation(0).degree[sites])))
        words = guarded_words
    lookup = {tuple(w): j for j, w in enumerate(words)}
    embed = np.eye(len(words))[:, [lookup[tuple(w)] for w in guarded_words]]
    e, v = eigh(h)
    # Equality tolerance groups an entire unresolved near-degeneracy into one level.
    labels = np.cumsum(np.r_[False, np.diff(e) > 1e-8])
    energies = np.array([e[labels == i].mean() for i in range(1 + max(labels))])
    return h, words, embed, e, v, labels, energies


def filters(e, band, levels):
    yield "ideal", (band == 1).astype(float)
    for count in CLOCK_BINS:
        phase = DT * (e - levels[1])
        weights = (
            np.exp(-0.5j * (count - 1) * phase)
            * np.sinc(count * phase / (2 * np.pi))
            / np.sinc(phase / (2 * np.pi))
        )
        yield f"clock_{count}", weights


def finite_number(x):
    return float(np.round(x, 12))


def sample(raw, rng):
    p = np.asarray(raw).ravel().clip(min=0)
    check(abs(p.sum() - 1) < 1e-8, "incomplete outcome distribution")
    p = np.round(p / p.sum(), 14)
    last = int(np.argmax(p))
    order = [j for j in range(len(p)) if j != last] + [last]
    p[last] = 1 - p[order[:-1]].sum()
    counts = np.zeros(len(p), dtype=np.int64)
    counts[order] = rng.multinomial(N, p[order])
    return counts


def interval(n, total):
    if not total:
        return None
    return [float(x) for x in clopper_pearson(n, int(total), TAIL)]


def meter_counts(wave, setups, filter_name, side, rng):
    # For each ordered source phase retain all six-port words, then its own
    # energy refusal. The final entry keeps all other central-source outcomes.
    outcomes = []
    for state, pair in zip(wave, setups, strict=True):
        _, words, embed, e, v, band, levels = pair[side]
        incoming = embed @ (state if side == 0 else state.T)
        selector = dict(filters(e, band, levels))[filter_name]
        after = v @ (selector[:, None] * (v.T @ incoming))
        p = (abs(after) ** 2).sum(axis=1)
        outcomes.extend(
            [*p, max(0.0, float(np.vdot(incoming, incoming).real - p.sum()))]
        )
    counts = sample([*outcomes, max(0.0, 1 - float(np.vdot(wave, wave).real))], rng)
    accepted = counts[:-1].reshape(2, -1)[:, :-1].sum(axis=0)
    admitted = int(accepted.sum())
    q = words[:, -2:].sum(axis=1) - 2
    signed = [int(accepted[q == val].sum()) for val in (-2, 0, 2)]
    return dict(
        counts=counts.tolist(),
        accepted=admitted,
        outer_charge_counts=signed,
        line_per_start_interval=interval(admitted, N),
        line_given_ordered_interval=interval(admitted, N - int(counts[-1])),
        split_given_line_interval=interval(signed[0] + signed[2], admitted),
    )


def observe(wave, setups, keeps, rng):
    spectrum = []
    work = []
    for state, pair, intact in zip(wave, setups, keeps, strict=True):
        vl, vr = [p[4].T @ p[2] for p in pair]
        table = abs(vl @ state @ vr.T) ** 2
        left_levels, right_levels = pair[0][5], pair[1][5]
        joint = np.zeros((len(pair[0][6]), len(pair[1][6])))
        np.add.at(joint, (left_levels[:, None], right_levels[None, :]), table)
        spectrum.extend(joint.ravel())
        phase_work = []
        for side, (op, old) in enumerate(zip(pair, intact, strict=True)):
            psi = state if side == 0 else state.T
            change = op[2].T @ op[0] @ op[2] - old[0]
            phase_work.append(float(np.vdot(psi, change @ psi).real))
        work.append(phase_work)
    admitted = float(np.vdot(wave, wave).real)
    joint = sample([*spectrum, max(0.0, 1 - admitted)], rng)
    return dict(
        joint_energy_counts=joint.tolist(),
        mean_switch_work_given_ordered=(np.sum(work, axis=0) / admitted).tolist(),
        local=[
            {
                name: meter_counts(wave, setups, name, side, rng)
                for name in ("ideal", *(f"clock_{b}" for b in CLOCK_BINS))
            }
            for side in (0, 1)
        ],
    )


def calibration(setups):
    out = []
    for pair in setups:
        sides = []
        for _, words, _, e, v, band, levels in pair:
            ground, first = v[:, band == 0], v[:, band == 1]
            # Maximum ground-to-first signed-port matrix element, well-defined at degeneracy.
            signed = [
                np.linalg.norm(ground.T @ ((words[:, j] - 1)[:, None] * first), 2)
                for j in range(6)
            ]
            readings = {
                name: [
                    finite_number(np.mean(abs(a[band == j]) ** 2))
                    for j in range(len(levels))
                ]
                for name, a in filters(e, band, levels)
            }
            sides.append(
                dict(
                    energies=levels.tolist(),
                    excitation_costs=(levels - levels[0]).tolist(),
                    signed_port_transition_norm=signed,
                    word_labels=words.tolist(),
                    clock_pass_probabilities_by_level=readings,
                )
            )
        out.append(sides)
    return out


def source_scale(k, rng):
    # Independent local contact-energy meter: equal and opposite 0/2 port loads.
    local = experiment.term(k)
    scale = []
    for a, b in [(0, 0), (0, 2)]:
        i = experiment.LAW.SRC.index((a, b))
        mu = float(local[i, i])
        # Resolve squared contrast sharply, then measure exchange-odd parity separately.
        swap = mu - k * (a - b) ** 2
        check(-1e-12 <= swap <= 1 + 1e-12, "invalid swap parity")
        high = int(rng.binomial(N, np.clip(swap, 0, 1)))
        lo, hi = interval(high, N)
        scale.append(
            dict(
                contact=[a, b],
                contrast_energy=k * (a - b) ** 2,
                odd_count=high,
                mean=k * (a - b) ** 2 + high / N,
                low=k * (a - b) ** 2 + lo,
                high=k * (a - b) ** 2 + hi,
            )
        )
    return dict(
        reads=scale,
        two_wall_surplus=2 * (scale[1]["mean"] - scale[0]["mean"]),
        interval=[
            2 * (scale[1]["low"] - scale[0]["high"]),
            2 * (scale[1]["high"] - scale[0]["low"]),
        ],
    )


def run():
    rng = np.random.default_rng(SEED)
    boundary, initial, ordered = reference()
    arms = []
    for k, cut in COUNTER_ARMS:
        h = csr_matrix(walls.operators(k, cut, 3))
        evolution = -32j * (h - eye(len(initial)) * np.mean(h.diagonal()))
        propagated = expm_multiply(evolution, initial).reshape(16, 8, 8, 2)[ordered]
        conditions = {}
        for mode in SETTINGS:
            setup = [
                [
                    wing(k, side, int(b if phase == 0 else 2 - b), cut, mode)
                    for side, b in enumerate(boundary)
                ]
                for phase in range(2)
            ]
            if mode == "keep":
                intact = setup
            conditions[mode] = dict(
                calibration=calibration(setup),
                sources={
                    name: observe(propagated[:, :, :, src], setup, intact, rng)
                    for src, name in enumerate(("sham", "wall"))
                },
            )
        arms.append(
            dict(
                k=k,
                cut_in_flight=cut,
                source_scale=source_scale(k, rng),
                conditions=conditions,
            )
        )
    sizes = {}
    for n in (1, 2, 3):
        hs = walls.wing_ref(8, 0, True, n)[0]
        e, v = eigh(hs)
        words = np.array([walls.loads(bits) for bits in product((0, 1), repeat=n)])
        sizes[n] = dict(
            gap=float(e[1] - e[0]),
            port_transition=(abs(v[:, 0] @ ((words - 1) * v[:, 1, None]))).tolist(),
        )
    return dict(
        shots=N,
        seed=SEED,
        exposure=32,
        wing_ports=wall_flight.extended_geometry()[2],
        charge_labels=[-2, 0, 2],
        release_lengths={str(k): v for k, v in sizes.items()},
        arms=arms,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = run()
    if args.output:
        with args.output.open("x") as f:
            json.dump(report, f, indent=1, allow_nan=False)
            f.write("\n")
    else:
        compare(
            report, json.loads((HERE / "records.json").read_text()), "daughter release"
        )
        print("Replayed the complete released-wing energy and signed-load records.")
