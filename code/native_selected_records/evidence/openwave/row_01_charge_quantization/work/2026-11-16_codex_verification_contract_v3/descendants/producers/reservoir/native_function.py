
# V3 verification contract: these checks are unconditional, also under python -O.
def _require(condition, message="Scientific check failed"):
    if not condition:
        raise AssertionError(message)

def native(loads, seam, coin):
    a, b = seam
    d = loads[a] - loads[b]
    current = (d - (2 * coin - 1) * (abs(d) % 2)) // 2
    out = list(loads)
    out[a] -= current
    out[b] += current
    return tuple(out), (coin, current)
