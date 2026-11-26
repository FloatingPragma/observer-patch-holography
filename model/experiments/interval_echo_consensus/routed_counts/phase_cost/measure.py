"""Read the material quadratic load through an empty-arm matter interferometer."""

import argparse
import json
from functools import lru_cache
from math import cos, log, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import federation
from model.experiments.interval_echo_gauge.arrival_ladder.measure import read_loads
from ..material import Routing, attempt, certify

HERE = Path(__file__).resolve().parent
SHOTS = 65536
ANGLES = (1.0, 0.25, 0.375)
FORMS = ("native", "linear_private", "exchange_only")


@lru_cache(maxsize=None)
def packet_rates(form):
    """Evolve both endpoints of the blank seam; sum final addresses into ternary bins."""
    if form not in FORMS:
        raise ValueError(form)
    exchange, contrast = native_law.local_term_parts()
    k = 0 if form == "exchange_only" else KAPPA
    full = exchange + k * contrast
    if np.max(abs(full[9:, :9])) > 1e-12:
        raise ArithmeticError("blank contact leaks")
    split = expm(-0.5j * pi * full[:9, :9])
    undo = expm(-1.5j * pi * full[:9, :9])
    private = np.diag(contrast[:9, :9])[
        [native_law.SRC.index((i, 0)) for i in range(3)]
    ]
    if form == "linear_private":
        private = np.arange(3)
    first, second = np.array(native_law.SRC).T
    banks = []
    for fraction in ANGLES:
        diagonal = np.exp(-1j * pi * fraction / KAPPA * k * private[first])
        read = undo @ (diagonal[:, None] * split)
        rates = []
        for q in range(3):
            probs = abs(read[:, native_law.SRC.index((q, 0))]) ** 2
            bins = np.bincount(second, weights=probs, minlength=3)
            np.testing.assert_allclose(bins.sum(), 1, atol=1e-12, rtol=0)
            if q and bins[3 - q] > 1e-12:
                raise ArithmeticError("the detector resolved a third load size")
            rates.append(float(bins[1:].sum()))
        banks.append(rates)
    return banks


def settings(board, routing, source, receiver):
    start = 12 * source
    meter = next(
        (a, b)
        for a, b in board.seams
        if a // 12 == b // 12 == receiver and board.degree[a] == 5
    )
    # The two meter addresses are prepared empty outside the complete source word.
    paths = [routing.transpose(start + i, meter[0]) for i in range(12)]
    if any(tuple(sorted(edge)) not in board.seams for path in paths for edge in path):
        raise ArithmeticError("route leaves the federation")
    return meter, paths


def histories(board, routing, source, reserved):
    start = 12 * source
    ports = list(range(start, start + 12))
    loads = [0] * board.n
    loads[start : start + 12] = [2] * 6 + [0] * 6
    # Both the splitter and the terminal echo require untouched blank heads.
    chords = [c for c in routing.chords if c not in reserved]
    cases = [("initial", loads.copy())]
    events = []
    offers = [(0, 1)] + [(j, j + 6) for j in range(6)]
    for j, (a, b) in enumerate(offers):
        pair = (ports[a], ports[b])
        before = [loads[i] for i in pair]
        chord = chords[len([ev for ev in events if ev["accepted"]])]
        accepted, work = attempt(loads, routing, pair, chord)
        after = [loads[i] for i in pair]
        events.append(
            dict(
                event=j,
                pair=pair,
                chord=chord,
                before=before,
                after=after,
                accepted=accepted,
                routing_contacts=work,
            )
        )
        cases.append(("refusal" if j == 0 else f"repair_{j}", loads.copy()))
    direct = [0] * board.n
    direct[start : start + 12] = [1] * 12
    cases.append(("direct_units", direct))
    return cases, events


def acquire(rng, form, loads, paths, meter):
    data = [[] for _ in ANGLES]
    for path in paths:
        copy = loads.copy()
        for a, b in path:
            copy[a], copy[b] = copy[b], copy[a]
        q, empty = [copy[i] for i in meter]
        if empty:
            raise ArithmeticError("reader bypass is not empty")
        for rows, rates in zip(data, packet_rates(form), strict=True):
            # The exact sector forbids a third word, including for q=0.
            p = round(float(np.clip(rates[q], 0, 1)), 12)
            n = int(rng.binomial(SHOTS, p))
            bins = [SHOTS - n, 0, 0]
            bins[q] += n
            rows.append(bins)
    return data


def rates_from_counts(banks):
    sums = np.asarray(banks, dtype=np.int64).sum(axis=1)
    return (sums[:, 1:].sum(axis=1) / sums.sum(axis=1)).tolist()


def decode(banks, a):
    first, second, held = rates_from_counts(banks)
    double = second - a * first
    unit_held = (1 - cos(ANGLES[2] * pi)) / 2
    double_held = (1 - cos(4 * ANGLES[2] * pi)) / 2
    return dict(
        quadratic=round(12 * (first + 4 * double), 8),
        held_frequency=round(held, 8),
        held_forecast=round(unit_held * first + double_held * double, 8),
    )


def run():
    certify()
    if not native_law.exact_checks()["checks"][
        "reference_pair_gives_n_squared_at_reference_zero"
    ]:
        raise ArithmeticError("unverified private zero reference")
    if read_loads() != [0, 1, 2]:
        raise ArithmeticError("signed terminal echo does not resolve all loads")
    board = federation(0)
    routing = Routing(board)
    rng = np.random.default_rng(112607)
    appointments = []
    # Each pooled rate averages 12*SHOTS independent Bernoulli trials (stratified by address).
    eps = sqrt(log(2 * 512 / 0.01) / (2 * 12 * SHOTS))
    for src, recv in ((0, 13), (13, 0)):
        meter, paths = settings(board, routing, src, recv)
        terminal_ref = next(
            v for v in board.neighbours[meter[1]] if v // 12 == recv and v != meter[0]
        )
        reserved = {meter, tuple(sorted((meter[1], terminal_ref)))}
        cases, events = histories(board, routing, src, reserved)
        standards = {}
        for q in (0, 1):
            word = [0] * board.n
            word[12 * src : 12 * (src + 1)] = [q] * 12
            standards[str(q)] = acquire(rng, "native", word, paths, meter)
        a = rates_from_counts(standards["1"])[1]
        observations = []
        for label, loads in cases:
            banks = {form: acquire(rng, form, loads, paths, meter) for form in FORMS}
            observations.append(
                dict(
                    case=label,
                    counts=banks,
                    decoded={key: decode(bank, a) for key, bank in banks.items()},
                )
            )
        appointments.append(
            dict(
                source=src,
                receiver=recv,
                meter=meter,
                terminal_empty_port=terminal_ref,
                splitter_contacts=2,
                terminal_echo_calls=2,
                routing_contacts=[len(path) for path in paths],
                events=events,
                standards=standards,
                single_B=round(a, 8),
                cases=observations,
            )
        )
    # For a product of empirical rates in [0,1], the error is at most the sum of the two radii.
    # The exact a is below 1; using that looser bound avoids any data-dependent coverage radius.
    return dict(
        address_shots=SHOTS,
        appointments_in_pi_over_kappa=ANGLES,
        simultaneous_family=512,
        pooled_radius=round(eps, 8),
        quadratic_radius=round(12 * 13 * eps, 8),
        held_difference_radius=round(
            (2 + 3 * (1 - cos(4 * pi * ANGLES[2])) / 2) * eps, 8
        ),
        readers=appointments,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    contents = json.dumps(run(), indent=1) + "\n"
    if parser.parse_args().write:
        (HERE / "records.json").write_text(contents)
    else:
        print(contents, end="")
