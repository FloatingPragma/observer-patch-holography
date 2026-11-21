"""Open a native guide detector at ticks posted by a local stationary Echo patch."""

import argparse
import json
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import nbinom

from model.candidate_laws.interval_echo import clock
from model.executable_model import federation
from .experiment import Guide, HERE, SHOTS, portable, select_route

THRESHOLDS = np.arange(33)
PERIODS = np.arange(193)


def layout(car, route):
    occupied = {p // 12 for p in route}

    def neighbouring_cells(c):
        return sorted(
            {p // 12 for a in range(12 * c, 12 * c + 12) for p in car.neighbours[a]}
            - {c}
        )

    for index, rail in enumerate(route[1::2]):
        cell = rail // 12
        vacancies = [
            p for p in car.neighbours[rail] if p // 12 == cell and p not in route
        ]
        if not vacancies:
            continue
        for c in neighbouring_cells(cell):
            if c in occupied:
                continue
            bank = [c]
            for parent in bank:
                for neighbour in neighbouring_cells(parent):
                    if neighbour not in occupied and neighbour not in bank:
                        bank.append(neighbour)
                if len(bank) >= 4:
                    return index, vacancies[0], bank[:4]
    raise ValueError("A physical readout clock could not be attached")


def timing_witness(car, front, bank, b_rail, vacancy):
    reg = clock.Registers(car, seed=306140)
    owner = {front, *bank, b_rail // 12}
    reg.native(b_rail, vacancy, owner)
    reg.prepare_front(front)
    digits = [tuple(12 * c + i for i in d) for c in bank for d in clock.DIGITS]
    publication = tuple(12 * bank[0] + i for i in clock.EXTRA_DIGITS[0])
    reg.initialize_digits(digits + [publication])
    crossings = [0]
    cycle = 0
    while len(crossings) <= 32:
        cycle += 1
        reg.publish_tick(publication, reg.fair(front, owner), owner)
        if reg.read_tick(publication, owner):
            reg.post(front, digits, owner)
            if reg.reading(front, digits, owner) != len(crossings):
                raise AssertionError("local count and hand disagree")
            crossings.append(cycle)
    assert np.allclose(list(reg.rates.values()), 0.5, atol=1e-12)
    return crossings


def decoding_curve(omega, probability):
    z = np.exp(2j * omega)
    transformed = (probability * z / (1 - (1 - probability) * z)) ** THRESHOLDS
    return (1 - transformed.real) / 2


def decode(data, p):
    def score(w):
        v = decoding_curve(w, p)
        a = max(0.0, (v @ data) / (v @ v))
        return np.linalg.norm(a * v - data) ** 2

    grid = np.linspace(0.005, np.pi / 2 - 0.005, 1025)
    i = int(np.argmin([score(w) for w in grid]))
    return minimize_scalar(
        score,
        bounds=(grid[max(0, i - 1)], grid[min(i + 1, len(grid) - 1)]),
        method="bounded",
        options={"xatol": 1.0e-12},
    ).x


def measure():
    car = federation(2)
    route, _ = select_route(car)
    guide = Guide(car, route, 0.2)
    station, vacancy, banks = layout(car, route)
    b_rail = route[2 * station + 1]
    witness = timing_witness(car, banks[0], banks[1:], b_rail, vacancy)
    cells = sorted({b_rail // 12, *banks})
    regs = clock.Registers(car, 0)
    probs = list(regs.rates.values())
    assert abs(probs[0] - probs[1]) < 1.0e-12
    p = float(np.mean(probs))
    waiting = np.zeros((len(THRESHOLDS), len(PERIODS)))
    waiting[0, 0] = 1
    for n in THRESHOLDS[1:]:
        waiting[n, n:] = nbinom.pmf(PERIODS[n:] - n, n, p)
    tail = float(nbinom.sf(PERIODS[-1] - 32, 32, p))
    assert tail < 1.0e-15
    main = json.loads((HERE / "results.json").read_text())
    primary = next(b for b in main["trials"][-1]["bands"] if b["detuning"] == 0.2)
    target = round(0.23 * guide.size / (2 * np.pi))
    momenta = {row["harmonic"]: row["k_measured"] for row in primary["held"]}
    rng = np.random.default_rng(306141)
    records = []
    for harmonic in (0, -target, target):
        ket = guide.uniform(2 * np.pi * harmonic / guide.size)
        local = []
        for t in PERIODS:
            local.append(abs(ket[b_rail]) ** 2)
            if t < PERIODS[-1]:
                guide.step(ket)
        exact = waiting @ local
        counts = rng.binomial(SHOTS, exact)
        freq = decode(counts / SHOTS, p)
        records.append(
            dict(
                harmonic=harmonic,
                k=0 if harmonic == 0 else momenta[harmonic],
                click_probabilities=exact.tolist(),
                clicks=counts.tolist(),
                omega_period=freq,
                omega_mean_home_tick=freq / p,
            )
        )
    for r in records[1:]:
        r["relative_kg_error"] = (
            abs(np.hypot(records[0]["omega_period"], r["k"]) - r["omega_period"])
            / r["omega_period"]
        )
    return dict(
        cells=cells,
        periods=192,
        repairs_per_period=p,
        detector_rail=b_rail,
        vacant_probe=vacancy,
        timing_register_cells=banks,
        witness_first_hitting_periods=witness,
        omitted_probability_upper_bound_1e_minus_15=tail < 1.0e-15,
        log10_omitted_probability=np.log10(tail),
        shots_per_setting=SHOTS,
        records=records,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record", action="store_true")
    data = portable(measure())
    path = HERE / "patch_time.json"
    if cli.parse_args().record:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    elif data != json.loads(path.read_text()):
        raise SystemExit("Clock-triggered archive differs")
    print(
        [
            (r["harmonic"], r["omega_mean_home_tick"], r.get("relative_kg_error"))
            for r in data["records"]
        ]
    )
