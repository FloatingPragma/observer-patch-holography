"""Enumerate the two-seam strict repair relation and a locally arbitrated control on sharp doubles."""
from itertools import combinations, product

import numpy as np

from model.executable_model import cell
from ..contrast_memory_spherical_v1 import reader
from .candidate import native_law as law

# Order in the shared tensor is (collar A, collar B, left, middle, right).
PORTS = (1, 0, 10)
AXES = ((0, 2, 3), (1, 3, 4))
TOL = 1e-12


def _different(first, second):
    return float(np.linalg.norm(first - second))


def _apply(state, echo, axes):
    target = tuple(range(3))
    block = np.moveaxis(state, axes, target)
    evolved = (echo @ block.reshape(99, -1)).reshape(block.shape)
    return np.moveaxis(evolved, target, axes)


def _attempt(state, echo, seam, probe=False):
    """Write, project sharply on 11, and undo a miss. A readiness probe also undoes a hit."""
    axes = AXES[seam]
    written = _apply(state, echo, axes)
    region = [slice(None)] * state.ndim
    region[axes[1]] = region[axes[2]] = 1
    hit = np.zeros_like(written)
    hit[tuple(region)] = written[tuple(region)]
    miss = written - hit
    weights = [float(np.vdot(v, v).real) for v in (miss, hit)]
    outcome = int(weights[1] > .5)
    if abs(weights[outcome] - 1) > TOL or weights[1 - outcome] > TOL:
        raise ArithmeticError('the declared sharp preparation has a non-deterministic pointer')
    branch = (hit if outcome else miss) / np.sqrt(weights[outcome])
    if probe or not outcome:
        branch = _apply(branch, echo, axes)
        if _different(branch, state) > TOL:
            raise ArithmeticError('uncomputing the sharp pointer changed the source')
    return branch, outcome


def _read(state):
    weight = (abs(state)**2).sum(axis=(0, 1))
    loads = tuple(int(n) for n in np.unravel_index(weight.argmax(), weight.shape))
    if abs(weight[loads] - 1) > TOL:
        raise ArithmeticError('terminal loads are not sharp')
    # B_e is the 11 projector on the incident ports, and gate produces redundant local agreement words.
    flags = [int(loads[s:s + 2] == (1, 1)) for s in range(2)]
    words = [reader.gate(flag, 0) for flag in flags]
    if any((word >> 6) & 3 != 1 or (word >> 8) & 7 != 7 * flag for flag, word in zip(flags, words)):
        raise ArithmeticError('addressed endpoint reader did not agree')
    return dict(loads=list(loads), pair_flags=flags, public_words=words)


def _start(doubled):
    state = np.zeros((11, 11, 3, 3, 3))
    loads = [0, 0, 0]
    loads[doubled] = 2
    state[(0, 0, *loads)] = 1
    return state


def _priority(echo, car, first):
    """Sharp dual rail on 2,3. Port 4 is the adjacent vacant reader rail for port 2."""
    if not {(2, 3), (2, 4)} <= set(car.seams) or set(PORTS) & {2, 3, 4}:
        raise ArithmeticError('priority register or its reader overlaps the tested wedge')
    state = np.zeros((11, 3, 3, 3))  # pointer collar, port 2, port 3, vacancy 4
    state[(0, 2 * (1 - first), 2 * first, 0)] = 1
    axes = (0, 1, 3)
    written = _apply(state, echo, axes)
    hit = np.zeros_like(written)
    hit[:, 1, :, 1] = written[:, 1, :, 1]
    yes = float(np.vdot(hit, hit).real)
    choice = int(yes < .5)
    sharp = written - hit if choice else hit
    if abs(yes - (1 - first)) > TOL or _different(_apply(sharp, echo, axes), state) > TOL:
        raise ArithmeticError('priority read failed')
    return choice


def _run_schedule(source, echo, order, priority):
    state = source.copy()
    trace = []
    for seam in order:
        before = _read(state)['loads']
        ready = []
        if priority is not None:
            # The host requests one seam. Probe the two available addresses in the host's order.
            readiness = [0, 0]
            for address in order:
                state, readiness[address] = _attempt(state, echo, address, probe=True)
            ready = readiness
        if priority is not None and all(ready) and seam != priority:
            trace.append(dict(seam=seam, readiness=ready, outcome='deferred'))
            continue
        state, flag = _attempt(state, echo, seam)
        after = _read(state)['loads']
        energy = [sum(n*n for n in loads) for loads in (before, after)]
        if sum(before) != sum(after) or energy[0] - energy[1] != 2*flag:
            raise ArithmeticError('accepted step did not strictly lower the squared loads at fixed total')
        trace.append(dict(seam=seam, readiness=ready, outcome='accepted' if flag else 'quiet',
                          probability_quiet_hit=[1 - flag, flag], squared_load_before_after=energy))
    # Examine every possible one-step suffix. A quiet attempt also leaves the state vector unchanged.
    selectable = [_attempt(state, echo, s)[1] for s in (0, 1)]
    if any(selectable) or sum(e['outcome'] == 'accepted' for e in trace) != 1:
        raise ArithmeticError('reachable strict-repair relation differs from its terminal certificate')
    return state, dict(order=list(order), attempts=trace, terminal=_read(state), terminal_retries=selectable)


def run():
    car = cell()
    a, middle, b = PORTS
    if a not in car.neighbours[middle] or b not in car.neighbours[middle] or b in car.neighbours[a]:
        raise ArithmeticError('the tested triple is not a native open wedge')
    echo = law.echo()
    rows = []
    for doubled, priority in product(range(3), (None, 0, 1)):
        token = None if priority is None else _priority(echo, car, priority)
        source = _start(doubled)
        paths = [_run_schedule(source, echo, order, token) for order in ((0, 1), (1, 0))]
        states, records = zip(*paths)
        fidelity = abs(np.vdot(*states))**2
        agree = bool(records[0]['terminal']['public_words'] == records[1]['terminal']['public_words'])
        if abs(fidelity - int(agree)) > TOL or agree != (doubled != 1 or priority is not None):
            raise ArithmeticError('endpoint comparison failed')
        rows.append(dict(double_at=PORTS[doubled], priority=priority, initial=_read(source), schedules=list(records),
                         same_addressed_endpoint=agree, terminal_state_fidelity=round(float(fidelity), 12)))
    pairs = [(u, v) for n in car.neighbours for u, v in combinations(n, 2)]
    return dict(schema='oph.interval_echo.two_seam_consensus.v1', ports=list(PORTS), seams=[[a, middle], [middle, b]],
                priority_rails=[2, 3], priority_reader=[2, 4], native_wedges=len(pairs),
                native_open_wedges=sum(v not in car.neighbours[u] for u, v in pairs), rows=rows)
