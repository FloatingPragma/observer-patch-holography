"""Continuing integer source and a current-only terminal observer.

No coherent memory law is used. The body control changes only body events;
the terminal native probing service is identical in both experiments.
"""
import json
import random
from pathlib import Path
from native_pin import native

HERE = Path(__file__).resolve().parent
SEAMS = tuple(map(tuple, json.loads((HERE / 'seams.json').read_text())))
N = 12
CANDIDATE = 'A-native-multiload-unit-selection-v1'


def tree_edges():
    seen, queue, tree = {0}, [0], []
    for a in queue:
        for u, v in SEAMS:
            b = v if u == a else u if v == a else None
            if b is not None and b not in seen:
                seen.add(b)
                queue.append(b)
                tree.append((a, b))
    if len(seen) != N:
        raise ValueError('Disconnected observation graph')
    return tuple(tree)


TREE = tree_edges()


def body_event(x, edge, coin, law='native'):
    if law == 'native':
        return native(x, edge, coin)
    if law != 'swap_hold':
        raise ValueError(law)
    a, b = edge
    y = list(x)
    current = coin * (x[a] - x[b])
    y[a] -= current
    y[b] += current
    return tuple(y), (coin, current)


def terminal_probe(x, rng, cap=32, forced_coin=None):
    """The controller sees records; x is confined to the source service."""
    tape = []
    for a, b in TREE:
        first = None
        for attempt in range(cap):
            coin = rng.randrange(2) if forced_coin is None else forced_coin
            x, (_, current) = native(x, (a, b), coin)
            tape.append([a, b, coin, current])
            if first is None:
                first = coin
            elif coin != first:
                break
        else:
            return x, tape, 'inconclusive'
    return x, tape, 'conclusive'


def decode_probe(tape, cap=32):
    """Return pre-probe loads modulo uniform shift, or explicit inconclusive.

    No load, total, preparation label, or simulator checkpoint is an input.
    Incidence corrections refer only to the terminal probe's current tape.
    """
    delta, relative, cursor = [0] * N, [None] * N, 0
    relative[0] = 0
    for a, b in TREE:
        start_delta = delta[a] - delta[b]
        first_coin = first_current = None
        for attempt in range(cap):
            if cursor >= len(tape):
                return None
            aa, bb, coin, current = tape[cursor]
            cursor += 1
            if (aa, bb) != (a, b) or coin not in (0, 1):
                raise ValueError('Malformed probe address/coin')
            if type(current) is not int:
                raise ValueError('Malformed integer-source record')
            delta[a] -= current
            delta[b] += current
            if first_coin is None:
                first_coin, first_current = coin, current
            elif coin != first_coin:
                difference = 2 * first_current + current - start_delta
                relative[b] = relative[a] - difference
                break
        else:
            if cursor != len(tape):
                raise ValueError('Records after timeout')
            return None
    if cursor != len(tape):
        raise ValueError('Extra terminal records')
    return relative


def decode_body_path(body_tape, terminal_relative):
    if terminal_relative is None:
        return None
    x, reverse = list(terminal_relative), [list(terminal_relative)]
    for a, b, coin, current in reversed(body_tape):
        x[a] += current
        x[b] -= current
        reverse.append(list(x))
    return list(reversed(reverse))


def run_history(initial, seed, law='native', horizon=384, cap=32,
                forced_probe_coin=None):
    rng, x, tape = random.Random(seed), tuple(initial), []
    for _ in range(horizon):
        a, b = SEAMS[rng.randrange(len(SEAMS))]
        coin = rng.randrange(2)
        x, (_, current) = body_event(x, (a, b), coin, law)
        tape.append([a, b, coin, current])
    x, probe, outcome = terminal_probe(x, rng, cap, forced_probe_coin)
    decoded = decode_probe(probe, cap)
    path = decode_body_path(tape, decoded)
    spans = None if path is None else [max(y) - min(y) for y in path]
    return {'initial_audit_only': list(initial), 'seed': seed, 'law': law,
            'horizon': horizon, 'probe_cap': cap, 'body': tape, 'probe': probe,
            'probe_outcome': outcome, 'decoded_preprobe_relative': decoded,
            'body_spans_from_records': spans,
            'first_balanced_body_event': None if spans is None else
                next((i for i, span in enumerate(spans) if span <= 1), None)}
