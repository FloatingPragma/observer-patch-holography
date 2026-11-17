"""Conclusive configuration readout and basis preparation by native repairs.

This is an active, current-adaptive service on the fixed twelve-port,
total-two carrier. No random seam schedule, disclosed initial bit, blank,
stationary preparation, or independent auxiliary is used. Fair placement
coins belong to the native C/E instruments. The controller reads only signed
currents and selects graph seams. Its covering sweep, when requested, also
retains currents: forgetting those would destroy retrodiction of collisions.

A nonzero current seeds two opposite known bits. At a known/unknown boundary,
try the edge r times. If unresolved, route the nearest opposite known bit to
that boundary and try again. Exactly one probe is mixed; a failed mixed probe
is never interpreted as an equal bit. Failed routing and unresolved probes
halt with an explicit timeout. A conclusive tape determines all final bits
and, by continuity, one original configuration. All source Kraus maps are
partial injections, so every conclusive fine word is a scalar |y><x|.

Binary restoration therefore implements a computational-basis measurement;
reset implements its measure-and-prepare version. For these two modes, summing
fine tapes gives success maps s_x |y_x><x| rho |x><y_x|, with y_x=x for restore
and one fixed target for reset. Readout alone retains a path-dependent mixture
of final work configurations; only its reported original label is fixed, and
it is not a Lu\u0308ders instrument on the retained work. The exact s_x below are
derived from the actual finite controller, not sampled. For arbitrary untouched
references, the complete restored/reset success/failure instrument is within
max_x(1-s_x) in half-diamond distance of its ideal success instrument. The same
bound holds for readout after discarding final work. Failure work and its flag
are kept in the reset/restoration comparison.
To compare reset with a blank tensor the ORIGINAL reference marginal, discard
the private resolved label and tape at that boundary. Conditioning on success
can reweight the reference. Full78 restoration is unavailable: binary native
repairs cannot recreate collisions. Full78 readout and binary reset do work.

Resource counts are native attempts, including unsuccessful probes, not
seconds or local incident ticks. This assumes joint access to the selected
seams' current records; it does not turn one passive local observer into that
controller. Coherent state synthesis and choice of conserved total remain
separate services.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from ..observer import native, quantum_mixing as quantum

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('configuration_results.json')
SOURCE_PATHS = (
    'model/native_clock/configuration.py', 'model/observer/native.py',
    'model/observer/quantum_mixing.py', 'model/common_source/slow_geometry.py',
    'evidence/native_geometric_source_20260925/spec.json',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
)
SEED_SCANS, RETRIES = 4, 32


def _configuration(value, *, binary):
    if (type(value) is not tuple or len(value) != 12
            or any(type(x) is not int or x < 0 or x > (1 if binary else 2)
                   for x in value) or sum(value) != 2):
        raise ValueError('an exact twelve-port total-two configuration is required')
    return value


def _parameters(seed_scans, retries, full, mode, target):
    if any(type(x) is not int or not 1 <= x <= 128 for x in (seed_scans, retries)):
        raise ValueError('seed scans and retries must be integers from 1 through 128')
    if type(full) is not bool:
        raise ValueError('full must be an original Boolean')
    if type(mode) is not str or mode not in ('measure', 'restore', 'reset'):
        raise ValueError('mode must be measure, restore, or reset')
    if full and mode == 'restore':
        raise ValueError('native repairs cannot restore full-sector collisions')
    if mode == 'reset':
        _configuration(target, binary=True)
    elif target is not None:
        raise ValueError('a target is used only for reset')


@lru_cache(maxsize=1)
def _current_maps():
    """Coarse current support, obtained from both actual placement outcomes."""
    result = {}
    for seam in quantum.seams():
        for current in (-1, 0, 1):
            mapping = {}
            for state in native.states():
                for coin in (0, 1):
                    out, (_, mark) = native.native(state, seam, coin)
                    if mark == current:
                        if state in mapping and mapping[state] != out:
                            raise ArithmeticError('signed current lost its unique state update')
                        mapping[state] = out
            result[seam, current] = mapping
    return result


def source_audit():
    """Exact Gram coverage and inverse-current identities on all 78 columns."""
    rows = {}
    for law in ('C', 'E'):
        branches = 0
        for seam in quantum.seams():
            groups = defaultdict(dict)
            marks = {}
            coverage = Counter()
            for state in native.states():
                for coin in (0, 1):
                    out, (_, current) = native.native(state, seam, coin)
                    label = current if law == 'C' else state[seam[0]]-state[seam[1]]
                    groups[coin, label][state] = out
                    if (coin, label) in marks and marks[coin, label] != current:
                        raise ArithmeticError('one fine label carries conflicting currents')
                    marks[coin, label] = current
                    rebuilt = list(out)
                    rebuilt[seam[0]] += current
                    rebuilt[seam[1]] -= current
                    if tuple(rebuilt) != state:
                        raise ArithmeticError('native signed-current inverse failed')
                    coverage[state] += 1
            if set(coverage.values()) != {2} or len(coverage) != 78:
                raise ArithmeticError('source instrument lost a Gram column')
            if any(len(set(m.values())) != len(m) for m in groups.values()):
                raise ArithmeticError('source branch is not a partial injection')
            branches += len(groups)
        rows[law] = {'fine_branches': branches, 'columns_per_seam': 78,
                     'exact_gram': 'I', 'current_inverse_verified': True}
    return rows


@lru_cache(maxsize=66)
def _routes(target):
    """Reverse breadth-first search using only actual source mixed-edge swaps."""
    routes = {target: ()}
    todo = deque([target])
    while todo:
        state = todo.popleft()
        for seam in quantum.seams():
            if state[seam[0]] == state[seam[1]]:
                continue
            moves = [(out, current) for coin in (0, 1)
                     for out, (_, current) in (native.native(state, seam, coin),)
                     if current]
            if len(moves) != 1:
                raise ArithmeticError('a binary mixed seam must have one moving coin')
            out = moves[0][0]
            if out not in routes:
                routes[out] = (seam,)+routes[state]
                todo.append(out)
    if len(routes) != 66:
        raise ArithmeticError('binary source routing graph is disconnected')
    return routes


def route(initial, target):
    """A deterministic shortest native swap route between two known states."""
    return _routes(_configuration(target, binary=True))[
        _configuration(initial, binary=True)]


@lru_cache(maxsize=1)
def routing_diameter():
    return max(len(path) for target in native.states(binary=True)
               for path in _routes(target).values())


@dataclass(frozen=True)
class Request:
    seam: tuple[int, int]
    limit: int
    phase: str


@dataclass(frozen=True)
class Outcome:
    success: bool
    reason: str
    initial: tuple[int, ...] | None
    final: tuple[int, ...] | None
    attempts: int
    currents: tuple[tuple[tuple[int, int], int], ...]


class _Tracker:
    """Public model support and known bits; never receives the actual loads."""
    def __init__(self, full):
        self.possible = {x: x for x in native.states(binary=not full)}
        self.known = {}
        self.currents = []

    def accept(self, request, currents):
        if (type(currents) is not tuple or not 1 <= len(currents) <= request.limit
                or any(type(j) is not int or j not in (-1, 0, 1) for j in currents)
                or any(currents[:-1])
                or (not currents[-1] and len(currents) != request.limit)):
            raise ValueError('a complete zero prefix ending in movement or cap is required')
        for current in currents:
            mapping = _current_maps()[request.seam, current]
            possible = {x: mapping[y] for x, y in self.possible.items() if y in mapping}
            if not possible:
                raise ValueError('current word is impossible on the declared source support')
            self.possible = possible
            self.currents.append((request.seam, current))
            a, b = request.seam
            if a in self.known:
                self.known[a] -= current
            if b in self.known:
                self.known[b] += current
            if any(bit not in (0, 1) for bit in self.known.values()):
                raise ValueError('current word violates a known binary bit')
        return currents[-1]

    def opposite(self, vertex):
        todo = deque([(vertex,)])
        seen = {vertex}
        while todo:
            path = todo.popleft()
            if self.known[path[-1]] != self.known[vertex]:
                return path
            for a, b in quantum.seams():
                other = b if a == path[-1] else a if b == path[-1] else None
                if other in self.known and other not in seen:
                    seen.add(other)
                    todo.append(path+(other,))
        raise ArithmeticError('connected known region lost its opposite bit')

    def finish(self, reason, initial=None):
        success = reason == 'complete'
        final = tuple(self.known[i] for i in range(12)) if success else None
        return Outcome(success, reason, initial, final, len(self.currents),
                       tuple(self.currents))


def policy(seed_scans=SEED_SCANS, retries=RETRIES, *, full=False,
           mode='measure', target=None):
    """Current-only generator. Send a tuple of currents for each Request.

    A request stops at its first movement or its finite cap. The generator's
    StopIteration.value is an Outcome; failure has no guessed initial label.
    All current records, including absorption, remain in that private outcome.
    """
    _parameters(seed_scans, retries, full, mode, target)
    tracker = _Tracker(full)

    def probe(seam, limit, phase):
        request = Request(seam, limit, phase)
        return tracker.accept(request, (yield request))

    if full:
        for seam in quantum.seams():
            yield from probe(seam, 1, 'absorb')
        if any(max(y) > 1 for y in tracker.possible.values()):
            raise ArithmeticError('one covering sweep did not absorb collisions')
    seeded = False
    for _ in range(seed_scans):
        for seam in quantum.seams():
            current = yield from probe(seam, 1, 'seed')
            if current:
                a, b = seam
                tracker.known[a], tracker.known[b] = int(current < 0), int(current > 0)
                seeded = True
                break
        if seeded:
            break
    if not seeded:
        return tracker.finish('seed_timeout')

    while len(tracker.known) < 12:
        seam = next(edge for edge in quantum.seams()
                    if sum(v in tracker.known for v in edge) == 1)
        vertex = next(v for v in seam if v in tracker.known)
        current = yield from probe(seam, retries, 'probe')
        if not current:
            path = tracker.opposite(vertex)
            for left, right in reversed(tuple(zip(path[:-1], path[1:], strict=True))):
                if tracker.known[left] == tracker.known[right]:
                    raise ArithmeticError('opposite-bit route contains an unmixed swap')
                if not (yield from probe(tuple(sorted((left, right))), retries, 'flip')):
                    return tracker.finish('flip_timeout')
            current = yield from probe(seam, retries, 'probe_after_flip')
        if not current:
            return tracker.finish('probe_timeout')
        a, b = seam
        tracker.known[a], tracker.known[b] = int(current < 0), int(current > 0)

    final = tuple(tracker.known[i] for i in range(12))
    if len(tracker.possible) != 1 or set(tracker.possible.values()) != {final}:
        raise ArithmeticError('conclusive record did not isolate one original column')
    initial = next(iter(tracker.possible))
    # Independent continuity reconstruction retains the covering-sweep current.
    rebuilt = list(final)
    for (a, b), current in tracker.currents:
        rebuilt[a] += current
        rebuilt[b] -= current
    if tuple(rebuilt) != initial:
        raise ArithmeticError('configuration retrodiction lost a signed current')
    destination = initial if mode == 'restore' else target if mode == 'reset' else final
    for seam in route(final, destination):
        if not (yield from probe(seam, retries, 'route')):
            return tracker.finish('route_timeout')
    if tuple(tracker.known[i] for i in range(12)) != destination:
        raise ArithmeticError('terminal native route missed its target')
    return tracker.finish('complete', initial)


@dataclass(frozen=True)
class Execution:
    outcome: Outcome
    output: tuple[int, ...]
    # Fine native environment records; label is current for C, precontrast for E.
    records: tuple[tuple[tuple[int, int], int, int], ...]
    law: str


def execute(initial, coins, seed_scans=SEED_SCANS, retries=RETRIES, *, full=False,
            mode='measure', target=None, law='C'):
    """Run the actual native source. The hidden initial state is not a policy input."""
    _parameters(seed_scans, retries, full, mode, target)
    state = _configuration(initial, binary=not full)
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('native law C or E required')
    stream = iter(coins)
    controller = policy(seed_scans, retries, full=full, mode=mode, target=target)
    request, records = next(controller), []
    while True:
        currents = []
        for _ in range(request.limit):
            try:
                coin = next(stream)
            except StopIteration as exc:
                raise ValueError('placement tape ended before the complete outcome') from exc
            if type(coin) is not int or coin not in (0, 1):
                raise ValueError('placement tape requires original integer bits')
            out, (_, current) = native.native(state, request.seam, coin)
            label = current if law == 'C' else state[request.seam[0]]-state[request.seam[1]]
            records.append((request.seam, coin, label))
            state = out
            currents.append(current)
            if current:
                break
        try:
            request = controller.send(tuple(currents))
        except StopIteration as done:
            return Execution(done.value, state, tuple(records), law)


def word_support(execution, *, full=False):
    """Exact actual C/E fine-word support; its squared amplitude is 2^-attempts."""
    if type(execution) is not Execution or type(full) is not bool:
        raise ValueError('an Execution and original Boolean domain are required')
    if type(execution.law) is not str or execution.law not in ('C', 'E'):
        raise ValueError('an original C/E execution law is required')
    if (type(execution.records) is not tuple
            or any(type(record) is not tuple or len(record) != 3
                   or type(record[0]) is not tuple
                   or any(type(v) is not int for v in record[0])
                   or record[0] not in quantum.seams()
                   or type(record[1]) is not int or record[1] not in (0, 1)
                   or type(record[2]) is not int
                   or record[2] not in ((-1, 0, 1) if execution.law == 'C'
                                        else (-2, -1, 0, 1, 2))
                   for record in execution.records)):
        raise ValueError('well-formed actual native fine records are required')
    mapping = {}
    for initial in native.states(binary=not full):
        state = initial
        for seam, coin, label in execution.records:
            out, (_, current) = native.native(state, seam, coin)
            actual = current if execution.law == 'C' else state[seam[0]]-state[seam[1]]
            if actual != label:
                break
            state = out
        else:
            mapping[initial] = state
    return mapping


@lru_cache(maxsize=1)
def _discovery_paths():
    """Exhaust all 600 possible binary input/first-moving-seam cases exactly.

    A successful mixed block has probability 1-2^-r and unique output. Equal
    blocks are certain zeros. Waiting before movement cannot affect this
    controller's later action. A failed mixed first probe cannot later succeed:
    after flipping its known endpoint the second probe is equal. Thus each seed
    has ONE successful macro path; omitted paths are precisely failure paths.
    """
    result = {}
    for initial in native.states(binary=True):
        rows, earlier = [], 0
        for seed_index, seed in enumerate(quantum.seams()):
            if initial[seed[0]] == initial[seed[1]]:
                continue
            controller = policy(1, 1)
            request, state, position, exponent = next(controller), initial, 0, 0
            while True:
                options = [native.native(state, request.seam, coin) for coin in (0, 1)]
                moving = [outcome for outcome in options if outcome[1][1]]
                use_move = (position == seed_index if request.phase == 'seed' else bool(moving))
                out, (_, current) = next(item for item in options
                                         if bool(item[1][1]) == use_move)
                if request.phase == 'seed':
                    position += 1
                elif current:
                    exponent += 1
                state = out
                try:
                    request = controller.send((current,))
                except StopIteration as done:
                    if not done.value.success or done.value.initial != initial:
                        raise ArithmeticError('a successful source macro path failed decoding')
                    rows.append((earlier, exponent, state))
                    break
            earlier += 1
        if earlier not in (8, 10):
            raise ArithmeticError('native total-two boundary cut changed')
        result[initial] = tuple(rows)
    if sum(map(len, result.values())) != 600:
        raise ArithmeticError('configuration policy audit did not cover all seed cases')
    return result


@lru_cache(maxsize=1)
def absorption_counts():
    """Exact complete covering-sweep law, including binary motion after absorption."""
    rows = []
    for initial in native.states():
        counts = Counter({initial: 1})
        for seam in quantum.seams():
            following = Counter()
            for state, weight in counts.items():
                for coin in (0, 1):
                    out, _ = native.native(state, seam, coin)
                    following[out] += weight
            counts = following
        if sum(counts.values()) != 2**30 or any(max(x) > 1 for x in counts):
            raise ArithmeticError('covering sweep lost mass or retained a collision')
        rows.append(tuple(counts.get(x, 0) for x in native.states(binary=True)))
    return tuple(rows)


@lru_cache(maxsize=128)
def _success_probabilities(seed_scans, retries, full, mode, target):
    survival = 1-F(1, 2**retries)
    values = []
    for initial, paths in _discovery_paths().items():
        cut = len(paths)
        repeat = sum((F(1, 2**(cut*j)) for j in range(seed_scans)), F())
        value = F()
        for earlier, exponent, final in paths:
            destination = initial if mode == 'restore' else target if mode == 'reset' else final
            exponent += len(route(final, destination))
            value += F(1, 2**(earlier+1))*survival**exponent
        values.append(repeat*value)
    if full:
        values = [sum((F(count, 2**30)*value for count, value in zip(row, values, strict=True)), F())
                  for row in absorption_counts()]
    return tuple(values)


def success_probabilities(seed_scans=SEED_SCANS, retries=RETRIES, *, full=False,
                          mode='measure', target=None):
    """Exact conclusive mass for every native.states(binary=not full) column."""
    _parameters(seed_scans, retries, full, mode, target)
    # Fraction's protected slots are writable; public values must not expose
    # the cache used by reset certificates and downstream preparation laws.
    return tuple(F(v.numerator, v.denominator)
                 for v in _success_probabilities(seed_scans, retries, full, mode, target))


def budget(seed_scans=SEED_SCANS, retries=RETRIES, *, full=False,
           mode='measure', target=None):
    _parameters(seed_scans, retries, full, mode, target)
    distance = 0 if mode == 'measure' else routing_diameter()
    return {'native_attempts_max': (30 if full else 0)+30*seed_scans+(75+distance)*retries,
            'failure_upper': min(F(1), F(1, 2**(8*seed_scans))+F(65+distance, 2**retries)),
            'terminal_swaps_max': distance}


def certificate(seed_scans=SEED_SCANS, retries=RETRIES, *, full=False,
                mode='measure', target=None):
    values = success_probabilities(seed_scans, retries, full=full, mode=mode, target=target)
    bounds = budget(seed_scans, retries, full=full, mode=mode, target=target)
    failure = 1-min(values)
    if not 0 <= failure <= bounds['failure_upper']:
        raise ArithmeticError('exact policy exceeds its uniform failure certificate')
    return {'seed_scans': seed_scans, 'retries': retries, 'full_sector': full,
            'mode': mode, 'target': list(target) if target is not None else None,
            'input_columns': len(values), 'successful_seed_cases_checked': 600,
            'native_attempts_max': bounds['native_attempts_max'],
            'terminal_swaps_max': bounds['terminal_swaps_max'],
            'failure_upper': str(bounds['failure_upper']),
            'exact_max_failure': str(failure),
            'exact_tracial_success_probability': str(sum(values, F())/len(values)),
            'success_table_sha256': hashlib.sha256('\n'.join(map(str, values)).encode()).hexdigest()}


def run():
    target = native.states(binary=True)[0]
    return {'schema': 'oph.native_configuration.v1', 'all_pass': True,
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                              for p in SOURCE_PATHS},
            'source_kraus_audit': source_audit(), 'routing_diameter': routing_diameter(),
            'binary_measurement': certificate(),
            'binary_restoration': certificate(mode='restore'),
            'binary_reset': certificate(mode='reset', target=target),
            'full78_reset': certificate(full=True, mode='reset', target=target),
            'contract': {
                'scheduler': 'deterministic current-adaptive graph-seam controller',
                'quantum_domain': 'native C or E, arbitrary total-two state and untouched reference',
                'full_sector': 'covering sweep currents retained for collision retrodiction',
                'success': 'rank-one fine words, exact configuration label, native target output',
                'failure': 'all timeouts retain their quantum output and public failure flag',
                'comparison': 'complete success/failure half-diamond error <= uniform failure mass',
                'blank_boundary': 'discard private original label and preparation tape',
                'conditional_reference_unchanged': False,
                'initial_bit_disclosure_or_independent_blank_required': False,
                'joint_current_access_and_seam_control_required': True,
                'restoration_of_collision_inputs': False,
                'physical_time_identified': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--check', nargs='?', const=str(OUTPUT))
    action.add_argument('--write', action='store_true')
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    payload = json.dumps(run(), indent=2, sort_keys=True)+'\n'
    if args.check:
        if Path(args.check).read_text(encoding='utf-8') != payload:
            raise SystemExit('native configuration receipt differs')
        print('native configuration receipt verified')
    elif args.write:
        with args.output.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(payload)
        print(f'wrote {args.output}')
    else:
        print(payload, end='')


if __name__ == '__main__':
    main()
