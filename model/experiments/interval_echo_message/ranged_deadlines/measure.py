"""Range three senders, then read their post-setup bits on one material hand."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import clock
from model.experiments.material_inertia.neutral_defect import ordered_index
from ..contact_deadline import measure as pilot

HERE = Path(__file__).resolve().parent
ROUTE = (1, 82, 72, 182)
COMMON = (78, 73, 183)
OTHER = (74, 77, 184)
HAND = (2, 3)
DEADLINES = (1, 2, 3)
HALF_WIDTH = float(np.sqrt(np.log(2 * 72 / 0.01) / (4 * pilot.SHOTS)))


def banks(h, car, y, tick, starts, marks=None):
    """One unnormalized joint register per initial single port and clock history."""
    state = np.zeros((len(y), len(starts)), dtype=complex)
    for col, port in enumerate(starts):
        state[ordered_index(HAND[0], port, car.n), col] = 1
    outputs = {d: {} for d in DEADLINES}
    prop = (-1j * pilot.STEP * h).tocsr()
    trace = prop.diagonal().sum()
    address = np.minimum(y, 12)
    remote = None if marks is None else y[:, None] == np.asarray(marks)[None, :]

    def enter(word, psi, mark):
        emerged = expm_multiply(prop, psi, traceA=trace)
        for digit, mask in [("0", ~tick), ("1", tick)]:
            record = word + digit
            part = mask[:, None] * emerged
            elapsed = len(record)
            hits = record.count("1")
            done = elapsed == pilot.SLOTS
            if digit == "1" and hits in outputs:
                mass = abs(part) ** 2
                hist = np.array(
                    [np.bincount(address, weights=col, minlength=13) for col in mass.T]
                )
                for j, p in enumerate(hist.T):
                    outputs[hits][f"{mark}:{record}:{j}"] = p
            if done:
                mass = (abs(part) ** 2).sum(axis=0)
                for d in DEADLINES:
                    if hits < d:
                        outputs[d][f"{mark}:{record}:timeout"] = mass
            elif hits < max(DEADLINES):
                if elapsed == 1 and remote is not None:
                    enter(record, part * remote, "visited")
                    enter(record, part * ~remote, "missed")
                else:
                    enter(record, part, mark)

    enter("", state, "unmarked")
    for bank in outputs.values():
        np.testing.assert_allclose(sum(bank.values()), 1, rtol=0, atol=2e-11)
    return outputs


def signal_samples(outputs, pair, rng):
    result = pilot.acquire(
        {d: {k: p[pair] for k, p in hist.items()} for d, hist in outputs.items()}, rng
    )
    return result


def passive_samples(outputs, col, rng):
    rows = {}
    for d, hist in outputs.items():
        keys = sorted(hist)
        probs = np.array([hist[k][col] for k in keys])
        counts = pilot.sample(probs, rng)
        returned = np.array(
            [
                k.startswith("visited:") and not k.endswith((":12", ":timeout"))
                for k in keys
            ]
        )
        rows[str(d)] = dict(
            keys=keys,
            counts=counts.tolist(),
            marked_returns=int(counts[returned].sum()),
            calculated_return_probability=float(probs[returned].sum().round(13)),
        )
    return rows


def addressed_ranges(car):
    # Sharp ports in the blank-memory sector. Audit every applied nontrivial pulse
    # directly against the full native local unitary before transporting the word.
    h0, h1 = clock.law.local_term_parts()
    gate = expm(-1j * np.pi * (h0 + clock.KAPPA * h1))
    rows = {}
    for end in range(1, len(ROUTE)):
        target = ROUTE[end]
        edges = list(zip(ROUTE[:end], ROUTE[1 : end + 1]))
        rows[str(target)] = {}
        for sham in (False, True):
            word = np.zeros(car.n, dtype=int)
            word[[HAND[0], ROUTE[0]]] = [2, 1]
            ticks = []
            local_readings = []
            mark = False
            first_return = None
            for i, (a, b) in enumerate(edges + edges[::-1], 1):
                if b not in car.neighbours[a]:
                    raise ValueError("invalid route")
                if not (sham and i == end):
                    m, n = word[[a, b]]
                    column = gate[:, clock.law.SRC.index((m, n))]
                    output = np.eye(clock.law.DIM)[:, clock.law.SRC.index((n, m))]
                    np.testing.assert_allclose(column, output, atol=1e-12, rtol=0)
                    word[[a, b]] = [n, m]
                ticks.append(int(word[list(HAND)].sum() == 2))
                local_readings.append(int(word[ROUTE[0]] == 1))
                if i == end:
                    mark = bool(word[target] == 1)
                if i > end and mark and local_readings[-1] and first_return is None:
                    first_return = sum(ticks)
            rows[str(target)]["outward_cut" if sham else "service"] = dict(
                tick_word=ticks,
                local_gateway=local_readings,
                remote_mark=mark,
                return_tick=first_return,
                deterministic_launches=pilot.SHOTS,
            )
    return rows


def run():
    pilot.native_gates()
    car, x, y, _, kinetic, contrast = pilot.flight()
    for j, zero in enumerate(ROUTE[1:]):
        common, one = COMMON[j], OTHER[j]
        if any(p not in car.neighbours[common] for p in (one, zero)):
            raise ArithmeticError("sender contact missing")
        if not zero // 12 == one // 12 == common // 12:
            raise ArithmeticError("encoder crosses a cell")
    tick = np.isin(x, HAND) & ~np.isin(y, HAND)
    ports = sorted(set(ROUTE + OTHER))
    print("native signal histories", flush=True)
    signal = banks(kinetic + contrast, car, y, tick, ports)
    print("changed-law histories", flush=True)
    changed = banks(0.5 * kinetic + contrast, car, y, tick, ports)
    print("passive ranging histories", flush=True)
    range_outputs = banks(kinetic + contrast, car, y, tick, [ROUTE[0]] * 3, ROUTE[1:])
    result = dict(
        protocol="ranged_deadlines_v1",
        shots_per_bit=pilot.SHOTS,
        threshold=0.52,
        simultaneous_half_width=HALF_WIDTH,
        addressed_ranging=addressed_ranges(car),
        stations={},
    )
    rng = np.random.default_rng(523071)
    for j, zero in enumerate(ROUTE[1:]):
        one = OTHER[j]
        station = {}
        for name, source, pair in [
            ("passive", signal, (zero, one)),
            ("relay", signal, (ROUTE[j], one)),
            ("sham", signal, (zero, zero)),
            ("half_exchange", changed, (zero, one)),
        ]:
            station[name] = signal_samples(source, [ports.index(p) for p in pair], rng)
            # Fresh sharp reset scan, before any free flight; the terminal
            # address is quiet for sources outside cell zero, exact on each bit.
            station[name]["reset"] = {
                "bit_addresses": [min(p, 12) for p in pair],
                "shots_per_bit": pilot.SHOTS,
            }
        station["passive_ranging"] = passive_samples(range_outputs, j, rng)
        result["stations"][str(zero)] = station
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    result = run()
    if args.record:
        with args.record.open("x") as f:
            json.dump(result, f, separators=(",", ":"))
            f.write("\n")
    elif result != json.loads((HERE / "records.json").read_text()):
        raise SystemExit("reacquired count banks differ")
    for station, data in result["stations"].items():
        print(
            station,
            {
                arm: [data[arm][str(d)]["correct"]["address"] for d in DEADLINES]
                for arm in ("passive", "relay", "sham", "half_exchange")
            },
        )
