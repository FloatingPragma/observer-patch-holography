"""Exact native() function excerpt from RER model/observer/native.py; see manifest and MIT license."""
def native(loads, seam, coin):
    a, b = seam
    d = loads[a] - loads[b]
    current = (d - (2 * coin - 1) * (abs(d) % 2)) // 2
    out = list(loads)
    out[a] -= current
    out[b] += current
    return tuple(out), (coin, current)
