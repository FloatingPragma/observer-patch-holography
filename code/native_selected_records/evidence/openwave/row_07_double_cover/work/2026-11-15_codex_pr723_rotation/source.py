"""Record source; all state reconstruction uses public preparation and currents."""
from pathlib import Path
import json
from native_endpoint import native

HERE = Path(__file__).resolve().parent
SEAMS = tuple(map(tuple, json.loads((HERE / 'seams.json').read_text())))
STATES = ((0,) * 12,) + tuple(tuple(int(i == p) - int(i == r) for i in range(12))
                            for p in range(12) for r in range(12) if p != r)
INDEX = {q: i for i, q in enumerate(STATES)}

def step(q, seam, coin, law):
    if law == 'native':
        return native(q, seam, coin)
    a, b = seam
    j = q[a] - q[b] if law == 'move_always' else 0
    out = list(q)
    out[a] -= j
    out[b] += j
    return tuple(out), (coin, j)

def read_next(preparation_record, seam, current):
    """Incidence reconstruction; no hidden source argument."""
    out = list(preparation_record)
    a, b = seam
    out[a] -= current
    out[b] += current
    return tuple(out)

def all_records(law):
    result = []
    for label, prepared in enumerate(STATES):
        for edge, seam in enumerate(SEAMS):
            for coin in (0, 1):
                actual, (_, current) = step(prepared, seam, coin, law)
                observed = read_next(prepared, seam, current)
                assert actual == observed
                result.append([label, edge, coin, current, INDEX[observed]])
    return result
