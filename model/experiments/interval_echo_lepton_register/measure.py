"""Charge-port interferometry on a marked seam of a twelve-port cell."""

import argparse
from functools import cache
import json
from math import acos, atan, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model.carrier import CELL_SEAMS, federation
from model.experiments.interval_echo_gauge.flag_measure.measure import (
    meter,
    path,
    routes,
)

N = 12
PORT_BASIS = np.eye(N, dtype=complex)
HERE = Path(__file__).resolve().parent
SHOTS = 2**22


@cache
def seam_block(speed):
    exchange, contrast = native_law.local_term_parts()
    h = speed * exchange + KAPPA * contrast
    indices = [native_law.SRC.index(word) for word in [(0, 1), (1, 0)]]
    others = np.delete(h[:, indices], indices, axis=0)
    np.testing.assert_allclose(others, 0, atol=6e-14)
    return h[np.ix_(indices, indices)]


@cache
def route(a, b):
    """Deterministic shortest path inside the first guarded cell."""
    pending = [(a,)]
    seen = {a}
    for trail in pending:
        if trail[-1] == b:
            return trail
        for port in federation(0).neighbours[trail[-1]]:
            if port < N and port not in seen:
                seen.add(port)
                pending.append((*trail, port))
    raise ValueError((a, b))


def contact(state, a, b, exposure, speed):
    if (min(a, b), max(a, b)) not in CELL_SEAMS:
        raise ValueError("An addressed contact requires a material seam")
    state[[a, b]] = expm(-1j * exposure * seam_block(speed)) @ state[[a, b]]


def phase(state, port, angle, speed, owner=0):
    """A vacancy gains +kappa*t relative to the occupied private reference."""
    board = federation(0)
    if not 0 <= owner < board.cells:
        raise ValueError("Phase owner outside the twenty-cell board")
    trail = min(
        (route(port, q) for q in range(N) if board.degree[owner * N + q] == 5), key=len
    )
    swaps = list(zip(trail, trail[1:]))
    for a, b in swaps:
        contact(state, a, b, pi, speed)
    exposure = (angle % (2 * pi)) / KAPPA
    state[trail[-1]] *= np.exp(1j * KAPPA * exposure)
    for a, b in reversed(swaps):
        contact(state, a, b, pi, speed)


def rotate(state, a, b, angle, speed, owner=0):
    """Native contact, flanked by private pulses, implements a real rail rotation."""
    trail = route(a, b)
    swaps = list(zip(trail[:-2], trail[1:-1]))
    for x, y in swaps:
        contact(state, x, y, pi, speed)
    x, y = trail[-2:]
    exposure = 2 * (angle % (2 * pi))
    phase(state, y, pi / 2, speed, owner)
    contact(state, x, y, exposure, speed)
    phase(state, x, (KAPPA + 0.5) * exposure, speed, owner)
    phase(state, y, (KAPPA + 0.5) * exposure - pi / 2, speed, owner)
    for x, y in reversed(swaps):
        contact(state, x, y, pi, speed)


def prepare(speed):
    trail = path()
    state = PORT_BASIS[:, trail[0]].copy()
    for i, (a, b) in enumerate(zip(trail, trail[1:])):
        rotate(state, a, b, acos(1 / sqrt(N - i)), speed)
    return state


def addresses(edge):
    a, b = edge
    neighbours = federation(0).neighbours
    tips = sorted(set(neighbours[a]) & set(neighbours[b]))
    if len(tips) != 2 or max(tips) >= N:
        raise ArithmeticError("An internal edge must have two internal face tips")
    return (a, b, *tips)


def reader(ports, angle, speed, owner=0):
    a, b, plus, minus = ports
    unitary = PORT_BASIS.copy()
    phase(unitary, plus, angle, speed, owner)
    phase(unitary, minus, -angle, speed, owner)
    rotate(unitary, a, b, -pi / 4, speed, owner)
    rotate(unitary, plus, minus, -pi / 4, speed, owner)
    rotate(unitary, a, plus, -atan(sqrt(2)), speed, owner)
    phase(unitary, plus, pi, speed, owner)
    phase(unitary, minus, pi / 2, speed, owner)
    rotate(unitary, plus, minus, -pi / 4, speed, owner)
    return unitary


def terminal_labels():
    """Test both external echo references, including every intermediate swap."""
    tags = meter()
    references = routes(0)
    tables = []
    for reference in references:
        labels = []
        for vacancy in range(N):
            records = []
            for port, trail in enumerate(reference["paths"]):
                word = [1] * N
                word[vacancy] = 0
                walk = list(zip(trail, trail[1:]))
                swaps = walk + walk[:-1][::-1]
                for x, y in swaps:
                    word[x], word[y] = word[y], word[x]
                # A double at the external reference lights the zero-load pointer.
                records.append(tags[word[trail[-1]]] + 2 == 2)
                for x, y in reversed(swaps):
                    word[x], word[y] = word[y], word[x]
                if word.index(0) != vacancy:
                    raise ArithmeticError("Reader disturbed the sharp source")
            if sum(records) != 1:
                raise ArithmeticError("Incomplete vacancy record")
            labels.append(records.index(True))
        tables.append(labels)
    np.testing.assert_array_equal(tables, [range(N)] * 2)
    return tables


def count_bank(rng, probabilities):
    """Simulate a complete multinomial bank on a fixed numerical probability grid."""
    weights = np.round(probabilities, 12)
    weights /= weights.sum()
    remaining_probability = 1.0
    remaining_starts = SHOTS
    counts = []
    for weight in weights[:-1]:
        if not remaining_starts:
            count = 0
        elif weight == 0:
            # Reserve this draw so numerical leakage into an ideal dark port
            # cannot change the random stream for subsequent detector bins.
            rng.random()
            count = 0
        else:
            conditional = min(1.0, weight / remaining_probability)
            count = rng.binomial(remaining_starts, conditional)
        counts.append(count)
        remaining_starts -= count
        remaining_probability -= weight
    return np.asarray([*counts, remaining_starts])


def free_release(ports, analyzer):
    board = federation(0)
    # Relative to the uniformly filled host, a vacancy costs kappa per
    # incident contrast and saves kappa for each vacant private contact.
    h = board.laplacian() / 2 + diags(KAPPA * (2 * board.degree - 6))
    a, _, plus, minus = ports
    modes = np.zeros((board.n, 3), complex)
    modes[:N] = analyzer.conj().T[:, [a, minus, plus]]
    columns = h @ modes
    leakage = columns - modes @ (modes.conj().T @ columns)
    read_times = [0, 1 / 8, 1 / 4, 1]
    rows = []
    rng = np.random.default_rng(15504)
    for time in read_times:
        evolved = expm_multiply((-1j * time) * h, modes)
        inside = analyzer @ evolved[:N]
        returned = (abs(inside[[a, minus, plus]]) ** 2).sum(axis=0)
        home = (abs(inside) ** 2).T
        outside = (abs(evolved[N:]) ** 2).sum(axis=0)
        complete = np.column_stack((home, outside))
        counts = [count_bank(rng, p).tolist() for p in complete]
        rows.append(
            dict(
                time=time,
                character_return=returned.tolist(),
                home_port_probabilities=home.tolist(),
                outside_home=outside.tolist(),
                home_and_no_vacancy_counts=counts,
            )
        )
    return dict(
        leakage_squared_eigenvalues=np.linalg.eigvalsh(
            leakage.conj().T @ leakage
        ).tolist(),
        records=rows,
    )


def measure():
    tables = terminal_labels()
    ports = addresses(CELL_SEAMS[0])
    a, b, plus, minus = ports
    state = prepare(1)
    density = abs(state) ** 2
    endpoint = density[a] + density[b]
    tip = density[plus]
    joint = np.outer(density, density)
    flag = np.arange(N) == plus
    agree = np.isin(np.arange(N), [a, b])
    xor = flag[:, None] ^ agree[None, :]
    angle = float(joint[xor].sum())
    rng = np.random.default_rng(15503)
    raw_pairs = count_bank(rng, joint.ravel()).reshape(N, N)
    measured_xor = int(raw_pairs[xor].sum()) / SHOTS
    rows = []
    for name, speed, theta in [
        ("count_phase", 1, angle),
        ("zero_phase", 1, 0),
        ("reversed", 1, -angle),
        ("incoherent", 1, angle),
        ("exchange_times_three", 3, angle),
        ("measured_count_phase", 1, measured_xor),
    ]:
        unitary = reader(ports, theta, speed)
        initial = prepare(speed)
        probabilities = (
            abs(unitary) ** 2 @ abs(initial) ** 2
            if name == "incoherent"
            else abs(unitary @ initial) ** 2
        )
        counts = count_bank(rng, probabilities)
        click = probabilities[[a, minus, plus]]
        rows.append(
            dict(
                setting=name,
                tip_phase_radians=theta,
                character_probabilities=click.tolist(),
                all_port_probabilities=probabilities.tolist(),
                counts=counts.tolist(),
                character_ratios=(click / min(click)).tolist(),
                intensity_quotient=float(sum(click) / np.sqrt(click).sum() ** 2),
            )
        )
    all_seams = []
    for edge in CELL_SEAMS:
        local = addresses(edge)
        ua, _, up, um = local
        matrix = reader(local, angle, 1)
        prob = abs(matrix @ state) ** 2
        all_seams.append(
            dict(seam=edge, character_probabilities=prob[[ua, um, up]].tolist())
        )
    return dict(
        schema=1,
        carrier_cells=20,
        charge_relative_to_filled_host=-1,
        ports=ports,
        character_ports=[a, minus, plus],
        external_readbacks=tables,
        branch_probabilities=[
            1 - endpoint - tip - density[minus],
            endpoint,
            tip,
            density[minus],
        ],
        independent_xor=angle,
        measured_xor=measured_xor,
        shared_address_xor=float(density[flag ^ agree].sum()),
        exposure_convention="one radian per unit of prior exclusive probability",
        paired_port_counts=raw_pairs.tolist(),
        shots_per_bank=SHOTS,
        banks=rows,
        all_internal_seams=all_seams,
        release=free_release(ports, reader(ports, angle, 1)),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    write = parser.parse_args().write
    result = measure()
    if write:
        (HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    else:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
