"""Two-matching matter-wave spectroscopy on actual Interval Echo seams."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar
from scipy.sparse import eye

from model.candidate_laws.interval_echo import candidate, clock
from model.candidate_laws.contrast_memory_spherical_v1.sectors import (
    doubled_hamiltonian,
)
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
SHOTS = 2**25
TIMES = np.arange(65)


def select_route(car):
    """Canonical DFS back-edge cycles, with no scores from the dynamics."""
    path, positions, seen = [0], {0: 0}, {0}
    stack = [(0, iter(car.neighbours[0]))]
    best = []
    while stack:
        _, neighbours = stack[-1]
        target = next(neighbours, None)
        if target is None:
            stack.pop()
            positions.pop(path.pop())
        elif target not in seen:
            seen.add(target)
            positions[target] = len(path)
            path.append(target)
            stack.append((target, iter(car.neighbours[target])))
        elif target in positions:
            length = len(path) - positions[target]
            if length <= len(best) or not 96 <= length <= 256 or length % 2:
                continue
            route = path[positions[target] :]
            if (
                sum(a // 12 != b // 12 for a, b in zip(route, route[1:] + route[:1]))
                >= 3
            ):
                best = route
    if not best:
        raise ValueError("No qualifying native route")
    for shift in range(len(best)):
        route = best[shift:] + best[:shift]
        # Two terminal rails need actual vacancies for their sharp detectors.
        read_vacancies = [
            sorted(set(car.neighbours[p]) - set(route)) for p in route[:2]
        ]
        if all(read_vacancies):
            assert all(
                b in car.neighbours[a] for a, b in zip(route, route[1:] + route[:1])
            )
            return route, [vs[0] for vs in read_vacancies]
    raise ValueError("No two adjacent receiver rails with spare native readers")


def native_matrix(duration, kappa=candidate.KAPPA):
    law = candidate.native_law
    a, b = law.local_term_parts()
    h = a + kappa * b
    basis = np.eye(99)[:, [law.SRC.index((2, 0)), law.SRC.index((0, 2))]]
    local = basis.T @ h @ basis
    assert np.linalg.norm(h @ basis - basis @ local) < 1.0e-11
    assert np.linalg.norm(h[:, law.SRC.index((0, 0))]) < 1.0e-12
    # A neighbouring vacant rail and the registered echo resolve this double.
    ok = np.arange(99) % 9 == law.SRC.index((1, 1))
    for load in (0, 2):
        p = np.linalg.norm(law.echo()[ok, law.SRC.index((load, 0))]) ** 2
        assert abs(p - int(load == 2)) < 1.0e-12
    return expm(-1j * duration * local)


class Guide:
    def __init__(self, car, route, delta, kappa=candidate.KAPPA):
        self.route = np.asarray(route)
        self.size = len(route) // 2
        self.car = car
        self.layers = [
            (self.route[::2], self.route[1::2], native_matrix(np.pi, kappa)),
            (
                self.route[1::2],
                np.roll(self.route[::2], -1),
                native_matrix(np.pi - 2 * delta, kappa),
            ),
        ]

    def step(self, ket):
        for a, b, u in self.layers:
            ket[a], ket[b] = u @ np.stack([ket[a], ket[b]])

    def uniform(self, k, polarization=(1.0, 0.0)):
        ket = np.zeros(self.car.n, complex)
        wave = np.exp(1j * np.arange(self.size) * k) / np.sqrt(self.size)
        ket[self.route[::2]] = polarization[0] * wave
        ket[self.route[1::2]] = polarization[1] * wave
        return ket

    def character_matrix(self, k):
        phase = np.exp(1j * k)
        _, _, first = self.layers[0]
        _, _, second = self.layers[1]
        # The second local gate acts on (B_n,A_{n+1}).
        return (
            np.array(
                [
                    [second[1, 1], second[1, 0] / phase],
                    [second[0, 1] * phase, second[0, 0]],
                ]
            )
            @ first
        )


def spatial_reader(guide, state):
    ket = state.copy()
    a, b, c = guide.route[:3]
    ket[[b, c]] = native_matrix(np.pi) @ ket[[b, c]]
    ket[[a, b]] = native_matrix(np.pi / 2) @ ket[[a, b]]
    pa, pb = abs(ket[[a, b]]) ** 2
    return np.array([pa, pb, 1 - pa - pb])


def spatial_inverse(outputs):
    left, right = outputs[:2]
    return np.arcsin(np.clip((right - left) / (right + left), -1.0, 1.0))


def temporal_inverse(data):
    """Free-amplitude sinusoidal frequency, exhaustive coarse grid + local refinement."""
    grid = np.linspace(0.003, np.pi / 2 - 0.003, 2049)
    waves = np.sin(grid[:, None] * TIMES) ** 2
    loss = np.sum(data**2) - (waves @ data) ** 2 / np.sum(waves**2, axis=1)
    index = int(np.argmin(loss))

    def objective(w):
        y = np.sin(TIMES * w) ** 2
        return np.sum((data - (data @ y) / (y @ y) * y) ** 2)

    return minimize_scalar(
        objective,
        bounds=(grid[max(0, index - 1)], grid[min(len(grid) - 1, index + 1)]),
        method="bounded",
        options={"xatol": 1.0e-13},
    ).x


def measure(guide, harmonic, stream):
    k = 2 * np.pi * harmonic / guide.size
    ket = guide.uniform(k)
    spatial = spatial_reader(guide, ket)
    counts = stream.multinomial(SHOTS, spatial)
    signal = []
    for t in TIMES:
        signal.append(abs(ket[guide.route[1]]) ** 2)
        if t != TIMES[-1]:
            guide.step(ket)
    if abs(np.vdot(ket, ket).real - 1) > 1.0e-10:
        raise ArithmeticError("guide lost the excitation")
    clicks = stream.binomial(SHOTS, signal)
    return dict(
        harmonic=harmonic,
        input_k=k,
        k_exact=spatial_inverse(spatial),
        k_measured=spatial_inverse(counts),
        omega_exact=temporal_inverse(np.array(signal)),
        omega_measured=temporal_inverse(clicks / SHOTS),
        spatial_probabilities=spatial.tolist(),
        spatial_counts=counts.tolist(),
        temporal_probabilities=signal,
        temporal_clicks=clicks.tolist(),
    )


def speed(guide):
    state = np.zeros(guide.car.n, complex)
    state[guide.route[0]] = 1
    for _ in range(6):
        guide.step(state)
    position = int(np.argmax(abs(state)))
    pair = list(guide.route[::2]).index(position)
    prob = abs(state[position]) ** 2
    assert prob > 1 - 1.0e-11
    return dict(
        source_port=int(guide.route[0]),
        endpoint_port=position,
        receipts=6,
        probability=prob,
        speed=pair / 6,
    )


def stationary_clock(car, home):
    # Separate timing preparation, on a connected four-cell observer patch.
    cells = [home]
    for c in cells:
        cells.extend(
            sorted(
                {v // 12 for p in range(12 * c, 12 * c + 12) for v in car.neighbours[p]}
                - set(cells)
            )
        )
        if len(cells) >= 4:
            break
    cells = cells[:4]
    owner = set(cells)
    reg = clock.Registers(car, seed=83)
    reg.prepare_front(home)
    digits = [tuple(12 * c + i for i in d) for c in cells[1:] for d in clock.DIGITS]
    reg.initialize_digits(digits)
    # Register in count-bank cell uses two rails left over from the three digits.
    publication = tuple(12 * cells[1] + i for i in clock.EXTRA_DIGITS[0])
    reg.initialize_digits([publication])
    readings = [0]
    for _ in range(64):
        mark = reg.fair(home, owner)
        reg.publish_tick(publication, mark, owner)
        # Read the publication locally; the controller coordinates the bank and ring in its patch.
        if reg.read_tick(publication, owner):
            reg.post(home, digits, owner)
        readings.append(reg.reading(home, digits, owner))
    assert np.allclose(list(reg.rates.values()), 0.5, atol=1.0e-12)
    return dict(
        cells=cells,
        rates=[float(rate) for rate in reg.rates.values()],
        material_elapsed_words=readings,
        expected_home_ticks_per_period=0.5,
    )


def packet(guide, centre, mass, c):
    ks = 2 * np.pi * np.fft.fftfreq(guide.size)
    target = ks[np.argmin(abs(ks - centre))]
    src_index = guide.size // 2
    state = np.zeros(guide.car.n, complex)
    for k in ks:
        distance = np.angle(np.exp(1j * (k - target)))
        if abs(distance) > 0.6:
            continue
        # The carrier phase is common to all momentum modes. The two-band splitting
        # chooses the time-directed branch; no KG frequencies enter its preparation.
        uk = guide.character_matrix(k)
        ref = guide.character_matrix(0)
        phase = np.trace(ref) / abs(np.trace(ref))
        values, vectors = np.linalg.eig(uk / phase)
        vec = vectors[:, np.argmin(np.angle(values))]
        dominant = 0 if target > 0 else 1
        vec *= np.exp(-1j * np.angle(vec[dominant]))
        weight = np.exp(-(distance**2) / (4 * 0.09**2) - 1j * k * src_index)
        state += weight * guide.uniform(k, vec)
    state /= np.linalg.norm(state)

    def distribution():
        return abs(state[guide.route[::2]]) ** 2 + abs(state[guide.route[1::2]]) ** 2

    zones = (
        np.angle(np.exp(2j * np.pi * (np.arange(guide.size) - src_index) / guide.size))
        * guide.size
        / (2 * np.pi)
    )
    before = distribution()
    for _ in range(8):
        guide.step(state)
    after = distribution()
    observed = (after - before) @ zones / 8
    predicted = c * c * target / np.hypot(mass, c * target)
    return dict(
        k=target,
        before=before.tolist(),
        after=after.tolist(),
        v_measured=observed,
        kg_group_prediction=predicted,
        absolute_error=abs(predicted - observed),
    )


def lattice_check(guide, delta):
    """Independent check on localized, broadband input (no fitted spatial or temporal mode)."""
    phase = np.trace(guide.character_matrix(0))
    phase /= abs(phase)
    state = np.zeros(guide.car.n, complex)
    state[guide.route[0]] = 1
    history = []
    errors = []
    for _ in range(17):
        history.append(state[guide.route].copy())
        if len(history) >= 3:
            previous, middle, last = history[-3:]
            translated = np.roll(middle, 2) + np.roll(middle, -2)
            errors.append(np.linalg.norm(last + previous - np.cos(delta) * translated))
        guide.step(state)
        state /= phase
    control = Guide(guide.car, guide.route, delta, kappa=0.0)
    differences = []
    for (_, _, original), (_, _, changed) in zip(guide.layers, control.layers):
        relative_phase = np.trace(changed.conj().T @ original) / 2
        differences.append(np.linalg.norm(original - relative_phase * changed))
    assert max(errors + differences) < 1.0e-10
    return dict(
        checked_localized_triples=len(errors),
        maximum_recurrence_residual=max(errors),
        exchange_only_maximum_layer_difference_modulo_phase=max(differences),
    )


def free_control(car):
    # The independent spatial operator is built from the fixed unweighted seam incidence.
    lap = car.laplacian()
    h = candidate.generator(car)
    offset = 6 * candidate.KAPPA
    diff = h - offset * eye(car.n) - lap / 2
    maximum = float(np.max(abs(diff.data), initial=0.0))
    assert maximum < 1.0e-11
    # Spectral moments give the squared-H vs affine prediction on the entire board.
    square_residual = h @ h - offset**2 * eye(car.n) - offset * lap - lap @ lap / 4
    pair_offset = 24 * candidate.KAPPA
    # The archived sector subtracts 12*kappa for two free constituents.
    same_double = doubled_hamiltonian(
        car, candidate.KAPPA
    ) + 12 * candidate.KAPPA * eye(car.n)
    double_residual = same_double - pair_offset * eye(car.n) - lap / 2
    assert np.max(abs(double_residual.data), initial=0.0) < 1.0e-11
    return dict(
        rest_offset=offset,
        matrix_max_error=maximum,
        squared_matrix_max_error=float(np.max(abs(square_residual.data), initial=0.0)),
        original_double_unshifted_rest_offset=pair_offset,
        original_double_matrix_max_error=float(
            np.max(abs(double_residual.data), initial=0.0)
        ),
        lambda_squared_coefficient=0.25,
        laplacian_squared_frobenius=float(np.linalg.norm((lap @ lap).data)),
    )


def evaluate():
    result = dict(
        candidate=candidate.NAME,
        shots_per_setting=SHOTS,
        trials=[],
        calibration="one dimer displacement and one rest beat per installed guide",
    )
    rng = np.random.default_rng(307091)
    for level in range(3):
        car = federation(level)
        route, vacancies = select_route(car)
        guide = Guide(car, route, 0)
        reference = speed(guide)
        size = guide.size
        board = dict(
            cells=car.cells,
            route=route,
            neighbouring_read_vacancies=vacancies,
            periods=len(TIMES) - 1,
            distinct_route_cells=len({p // 12 for p in route}),
            external_seams=sum(
                a // 12 != b // 12 for a, b in zip(route, route[1:] + route[:1])
            ),
            speed_calibration=reference,
            stationary_clock=stationary_clock(car, route[0] // 12),
            unfiltered_single_load_control=free_control(car),
            bands=[],
        )
        for delta in (0.1, 0.2, 0.8):
            guide = Guide(car, route, delta)
            targets = (0.115, 0.23, 0.35, 0.46) if delta < 0.8 else (1.15,)
            positive = {max(1, int(round(k * size / (2 * np.pi)))) for k in targets}
            harmonics = sorted(positive | {-i for i in positive})
            rest = measure(guide, 0, rng)
            measured = []
            for j in harmonics:
                row = measure(guide, j, rng)
                for suffix in ("exact", "measured"):
                    m = rest["omega_" + suffix]
                    k = row["k_" + suffix]
                    w = row["omega_" + suffix]
                    ck = reference["speed"] * k
                    row["kg_relative_error_" + suffix] = abs(np.hypot(m, ck) - w) / w
                    row["parabolic_relative_error_" + suffix] = (
                        abs(m + ck**2 / (2 * m) - w) / w
                    )
                measured.append(row)
            band = dict(
                detuning=delta,
                rest=rest,
                held=measured,
                analytic_controls=lattice_check(guide, delta),
                worst_kg_exact=max(x["kg_relative_error_exact"] for x in measured),
                worst_kg_sampled=max(x["kg_relative_error_measured"] for x in measured),
                worst_parabolic_exact=max(
                    x["parabolic_relative_error_exact"] for x in measured
                ),
            )
            if level == 2 and delta == 0.2:
                band["packets"] = [
                    packet(guide, k, rest["omega_exact"], reference["speed"])
                    for k in (-0.35, 0.35)
                ]
            board["bands"].append(band)
        result["trials"].append(board)
    return result


def portable(obj):
    if isinstance(obj, dict):
        return {k: portable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [portable(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        return round(float(obj), 9) + 0.0
    return obj


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = portable(evaluate())
    if args.output:
        args.output.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    elif json.loads((HERE / "results.json").read_text()) != data:
        raise SystemExit("row 30 record differs from archive")
    for board in data["trials"]:
        print(
            board["cells"],
            [
                (
                    b["detuning"],
                    b["worst_kg_exact"],
                    b["worst_kg_sampled"],
                    b["worst_parabolic_exact"],
                )
                for b in board["bands"]
            ],
        )
        for band in board["bands"]:
            for p in band.get("packets", []):
                print("packet", p["k"], p["v_measured"], p["kg_group_prediction"])


if __name__ == "__main__":
    main()
