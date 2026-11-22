"""Gateway-monitored preparation and terminal one-load records on the frozen boards."""

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from numbers import Integral, Real

import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from .experiment import dark_bands, free, gateway_table, port_reader


ROUNDS = 1024
LAUNCHES = tuple(range(12)) * 4
EXPOSURES = (0, 1, 7, 37)
PULSES = (0.0, 0.5, 1.0, 2.0)
ATOL = 2e-10


def _check_close(actual, expected, label):
    if not np.isfinite(actual).all() or not np.allclose(actual, expected, rtol=0, atol=ATOL):
        raise ArithmeticError(label)


def _index(value, low, high):
    return isinstance(value, Integral) and not isinstance(value, (bool, np.bool_)) and low <= value <= high


def _chart(gateways):
    tips = tuple(gateways)
    if len(tips) != 3 or any(not _index(p, 0, 11) for p in tips) or tuple(sorted(set(tips))) != tips:
        raise ValueError("expected three distinct, ascending home gateway ports")
    return tuple(map(int, tips))


def gateway_wiring(car, home):
    """(local port, exterior port, written sign) in the actual external seam orientation."""
    if not _index(home, 0, car.cells - 1):
        raise ValueError("home cell outside the board")
    wires = []
    for a, b in car.external:
        if a // 12 == home:
            wires.append((a % 12, b, 1))
        if b // 12 == home:
            wires.append((b % 12, a, -1))
    wires.sort()
    _chart(p for p, _, _ in wires)
    return wires


@lru_cache(maxsize=1)
def signed_seam_check():
    """Check the vacancy and both click Kraus columns, including the restored collar."""
    law = candidate.native_law
    echo = law.echo()
    loads = ((0, 0), (1, 0), (0, 1))
    blank = np.eye(law.DIM)[:, [law.SRC.index(pair) for pair in loads]]
    written = echo @ blank
    worst = 0.0
    for sign in (-2, -1, 0, 1, 2):
        labels = [False, *(contrast == sign for _, contrast in law.LABELS)]
        selected = np.repeat(labels, law.N_SRC)[:, None] * written
        recovered = echo.T @ selected
        reference = blank * [a - b == sign for a, b in loads]
        _check_close(recovered, reference, "signed seam read changes a blank-code column")
        worst = max(worst, float(np.max(abs(recovered - reference))))
    return worst


@lru_cache(maxsize=1)
def _internal_generator():
    """Thirty isolated native contacts, with external contacts and references idle."""
    local = cell()
    if len(local.seams) != 30 or local.external:
        raise ValueError("incorrect preparation cell")
    law = candidate.native_law
    one = [law.SRC.index(pair) for pair in ((1, 0), (0, 1))]
    exchange, native = np.zeros((2, 12, 12))
    for operator, kappa in ((exchange, 0.0), (native, candidate.KAPPA)):
        block = law.seam_term(kappa)[np.ix_(one, one)]
        for edge in local.seams:
            operator[np.ix_(edge, edge)] += block
    offset = 5 * candidate.KAPPA
    _check_close(native - exchange, offset * np.eye(12), "nonuniform isolated contrast")
    return exchange, offset


@lru_cache(maxsize=4)
def _exchange_pulse(duration):
    exchange, _ = _internal_generator()
    unitary = expm(-1j * duration * exchange)
    _check_close(unitary.conj().T @ unitary, np.eye(12), "preparation lost unitarity")
    return unitary


@lru_cache(maxsize=1)
def _reset_paths():
    """Route a post-click sharp load by full internal exchanges on the home graph."""
    local = cell()
    distances, parents = shortest_path(local.adjacency(), directed=False, unweighted=True,
                                       return_predecessors=True)
    law = candidate.native_law
    one = [law.SRC.index(pair) for pair in ((1, 0), (0, 1))]
    pulse = expm(-1j * np.pi * law.seam_term(candidate.KAPPA)[np.ix_(one, one)])
    paths = {}
    for start in range(12):
        for target in range(12):
            reverse = [target]
            while reverse[-1] != start:
                reverse.append(int(parents[start, reverse[-1]]))
            path = tuple(reversed(reverse))
            if len(path) - 1 != distances[start, target]:
                raise ArithmeticError("route length disagrees with the home graph")
            state = np.eye(12, dtype=complex)[:, start]
            for a, b in zip(path, path[1:]):
                if b not in local.neighbours[a]:
                    raise ArithmeticError("reset pulse outside the home graph")
                state[[a, b]] = pulse @ state[[a, b]]
            _check_close(abs(state) ** 2, np.eye(12)[target], "retry lost its single load")
            paths[start, target] = path
    return paths


@dataclass
class Preparation:
    gateways: tuple
    pulse: float
    launches: tuple
    rounds: int
    first_click: np.ndarray
    accepted_vectors: np.ndarray
    success_by_launch: np.ndarray
    attempt_reach: np.ndarray
    attempt_accepted: np.ndarray
    accepted_density: np.ndarray
    exhaustion: float
    routes: dict
    max_round_balance_error: float

    def attempt_joint(self, attempt):
        """Unconditional [round, gateway] first-click masses and acceptance at a 1-based attempt."""
        if not _index(attempt, 1, len(self.launches)):
            raise ValueError("attempt outside the source schedule")
        i = attempt - 1
        return self.attempt_reach[i] * self.first_click[:, :, self.launches[i]], self.attempt_accepted[i]

    def history_probability(self, clicks, *, accepted):
        """Mass of one full history. Each failed attempt specifies (1-based round, local gateway)."""
        if not isinstance(accepted, (bool, np.bool_)):
            raise ValueError("acceptance flag must be boolean")
        if len(clicks) > len(self.launches) or accepted != (len(clicks) < len(self.launches)):
            raise ValueError("history must end at first acceptance or at 48 failed attempts")
        weight = 1.0
        for attempt, (round_number, port) in enumerate(clicks):
            if not _index(round_number, 1, self.rounds):
                raise ValueError("first-click round outside the source schedule")
            if not _index(port, 0, 11) or port not in self.gateways:
                raise ValueError("click outside the monitored ports")
            index = self.gateways.index(port)
            weight *= self.first_click[round_number - 1, index, self.launches[attempt]]
        if accepted:
            weight *= self.success_by_launch[self.launches[len(clicks)]]
        return float(weight)


def prepare_layout(gateways, pulse=1.0, *, start=0, direction=1, rounds=ROUNDS):
    """Acquire [round, sorted gateway, sharp launch] records and the unnormalized release."""
    tips = _chart(gateways)
    if not isinstance(pulse, Real) or isinstance(pulse, (bool, np.bool_)) or pulse not in PULSES:
        raise ValueError("duration outside the frozen settings and the exchange-off control")
    if not _index(start, 0, 11) or not _index(direction, -1, 1) or direction == 0:
        raise ValueError("cyclic sweeps require a valid start and direction +1 or -1")
    if not _index(rounds, 1, ROUNDS):
        raise ValueError("clear rounds outside the preparation budget")
    launches = tuple(int((start + direction * i) % 12) for i in range(12)) * 4
    signed_seam_check()
    unitary = _exchange_pulse(float(pulse))
    state = np.eye(12, dtype=complex)
    clicks = np.empty((rounds, len(tips), 12))
    previous = np.ones(12)
    worst = 0.0
    for i in range(rounds):
        state = unitary @ state
        # Different sharp gateways are orthogonal. A vacancy removes only its own row.
        for j, port in enumerate(tips):
            clicks[i, j] = abs(state[port]) ** 2
            state[port] = 0
        clear = (abs(state) ** 2).sum(axis=0)
        error = abs(previous - clear - clicks[i].sum(axis=0))
        if not np.isfinite(error).all():
            raise ArithmeticError("nonfinite gateway branch")
        worst = max(worst, float(error.max()))
        previous = clear
    _check_close(worst, 0.0, "gateway round violates branch normalization")
    # Each clear round has the same scalar 5κ phase. Zero exchange holds the contacts for one unit.
    _, offset = _internal_generator()
    state *= np.exp(-1j * rounds * (pulse or 1.0) * offset)
    failure = clicks.sum(axis=(0, 1))
    success = clear
    _check_close(success + failure, np.ones(12), "one attempt has an incomplete record")
    reach = np.empty(len(launches))
    acceptance = np.empty_like(reach)
    survival = 1.0
    for i, launch in enumerate(launches):
        reach[i] = survival
        acceptance[i] = survival * success[launch]
        survival *= failure[launch]
    _check_close(acceptance.sum() + survival, 1.0, "incomplete source outcome partition")
    weights = np.bincount(launches, weights=reach, minlength=12)
    rho = (state * weights) @ state.conj().T
    _check_close(rho, rho.conj().T, "non-Hermitian accepted density")
    _check_close(np.trace(rho), acceptance.sum(), "released mixture lost attempt weights")
    routes = {(p, q): path for (p, q), path in _reset_paths().items() if p in tips}
    return Preparation(tips, float(pulse), launches, int(rounds), clicks, state, success, reach, acceptance,
                       rho, float(survival), routes, worst)


def _summary(prep):
    """Round and attempt factors in Preparation reconstruct every addressed history."""
    protection = sum(p for _, p in dark_bands(candidate.generator(cell()).toarray(),
                                              list(prep.gateways)))
    # At any exposure the projection onto this reducing subspace stays inside home.
    joint_floor = float(np.trace(protection @ prep.accepted_density).real)
    _check_close(protection @ protection, protection, "invalid local protected projector")
    attempt_click, _ = prep.attempt_joint(len(prep.launches))
    fates = attempt_click.sum(axis=0)
    _check_close(fates.sum(), prep.exhaustion, "last click missing from exhaustion")
    max_route = max(len(path) - 1 for path in prep.routes.values())
    max_contacts = sum(max(len(prep.routes[p, q]) - 1 for p in prep.gateways)
                       for q in prep.launches[1:])
    return dict(
        success_by_sharp_launch=prep.success_by_launch.tolist(),
        first_click_by_gateway_and_launch=prep.first_click.sum(axis=0).tolist(),
        attempt_reach=prep.attempt_reach.tolist(),
        attempt_acceptance=prep.attempt_accepted.tolist(),
        acceptance=float(prep.attempt_accepted.sum()),
        exhaustion=prep.exhaustion,
        exhausted_at_gateway=fates.tolist(),
        maximum_round_balance_error=prep.max_round_balance_error,
        all_time_joint_home_floor=joint_floor,
        floor_exceeds_094=bool(joint_floor > 0.94 + ATOL),
        max_contacts_one_reset=max_route,
        max_reset_contacts_over_47_retries=max_contacts,
    )


def _detector_record(occupations, indices, port_effects, exhaustion):
    if not np.isfinite(occupations).all() or occupations.min() < -ATOL:
        raise ArithmeticError("invalid free-flight occupation")
    occupations = np.maximum(occupations, 0)
    ports = port_effects @ occupations[indices]
    exterior = float(occupations[:indices.start].sum() + occupations[indices.stop:].sum())
    joint_home = float(ports.sum())
    _check_close(joint_home + exterior + exhaustion, 1.0, "incomplete source and terminal record")
    return dict(ports=ports.tolist(), exterior_quiet=exterior,
                source_exhaustion=exhaustion, joint_home=joint_home)


@lru_cache(maxsize=1)
def terminal_port_effects():
    """An internal unsigned seam read has one common record word on the entire exterior."""
    local = cell()
    # Each home port writes its incident seam bits. All exterior addresses write thirty zeros.
    words = {tuple(port in edge for edge in local.seams) for port in range(12)}
    if len(words) != 12 or (False,) * 30 in words:
        raise ArithmeticError("terminal seam words fail to resolve home and exterior")
    effects = np.array([port_reader(local, p)[1] for p in range(12)])
    _check_close(effects, np.eye(12), "terminal reader fails to distinguish the home ports")
    return effects


def scored_appointments(preparation, level=0, home=13, *, exposures=EXPOSURES):
    """Free federation flight from accepted, port-dephased, and uniform sharp ensembles."""
    if not _index(level, 0, 2):
        raise ValueError("level outside the registered boards")
    exposures = tuple(exposures)
    if not exposures or any(not isinstance(t, Real) or isinstance(t, (bool, np.bool_))
                            or not np.isfinite(t) or t < 0 for t in exposures):
        raise ValueError("appointments require finite nonnegative exposures")
    car = federation(level)
    wiring = gateway_wiring(car, home)
    if tuple(p for p, _, _ in wiring) != preparation.gateways:
        raise ValueError("source gateway layout differs from the actual cell")
    h = candidate.generator(car)
    port_effects = terminal_port_effects()
    indices = slice(12 * home, 12 * (home + 1))
    sharp = np.zeros((car.n, 12), dtype=complex)
    sharp[indices] = np.eye(12)
    weights = np.bincount(preparation.launches, weights=preparation.attempt_reach, minlength=12)
    accepted = preparation.accepted_vectors * np.sqrt(weights)
    diagonal = np.diag(preparation.accepted_density).real
    _check_close(accepted @ accepted.conj().T, preparation.accepted_density,
                 "free-flight inputs differ from the accepted mixture")
    # Independent Poisson SWAPs on every edge, each at rate |H_ab| = 1/2.
    heat = -0.5 * car.laplacian()
    seed = np.zeros(car.n)
    seed[indices] = diagonal
    readings = []
    for t in exposures:
        transfer = free(h, sharp, t)
        _check_close((abs(transfer) ** 2).sum(axis=0), np.ones(12),
                     "free flight lost a sharp column")
        sharp_prob = abs(transfer) ** 2
        flight = {
            "coherent": (abs(transfer @ accepted) ** 2).sum(axis=1),
            "dephased": sharp_prob @ diagonal,
            "unfiltered_uniform": sharp_prob.mean(axis=1),
            "stochastic_swap": expm_multiply(t * heat, seed) if t else seed,
        }
        row = {"t": t}
        for arm, occupation in flight.items():
            exhaustion = 0.0 if arm == "unfiltered_uniform" else preparation.exhaustion
            row[arm] = _detector_record(occupation, indices, port_effects, exhaustion)
        readings.append(row)
    return dict(cells=car.cells, home_cell=home, pulse=preparation.pulse,
                gateway_wiring=wiring, readings=readings,
                crosses_below_094={arm: any(row[arm]["joint_home"] < 0.94 - ATOL
                                             for row in readings) for arm in flight})


def gateway_acquisition(levels=(0, 1, 2), scored_level=0):
    levels = tuple(levels)
    if (any(not _index(level, 0, 2) for level in levels)
            or tuple(sorted(levels)) != (0, 1, 2) or not _index(scored_level, 0, 2)):
        raise ValueError("survey requires the three distinct registered levels and a valid scored level")
    boards = []
    layouts = set()
    for level in sorted(levels):
        car = federation(level)
        charts = Counter(tuple(p) for p in gateway_table(car))
        layouts.update(charts)
        boards.append(dict(cells=car.cells, occurrences=[
            dict(gateways=gw, count=n) for gw, n in sorted(charts.items())]))
    scored_car = federation(scored_level)
    scored = tuple(gateway_table(scored_car)[13])
    if scored not in layouts:
        raise ValueError("scored layout absent from the preparation survey")
    summaries = []
    terminal = {}
    for tips in sorted(layouts):
        settings = {}
        for pulse in PULSES:
            preparation = prepare_layout(tips, pulse)
            settings[str(pulse)] = _summary(preparation)
            if tips == scored:
                terminal[str(pulse)] = scored_appointments(preparation, scored_level)
        summaries.append(dict(gateways=tips, settings=settings))
    return dict(rounds=ROUNDS, launches=LAUNCHES, exposures=EXPOSURES,
                isolated_contrast_offset=_internal_generator()[1],
                maximum_preparation_pulses=len(LAUNCHES)*ROUNDS,
                maximum_gateway_reads=len(LAUNCHES)*ROUNDS*3,
                signed_seam_column_error=signed_seam_check(), boards=boards,
                layouts=summaries, scored_appointments=terminal)
