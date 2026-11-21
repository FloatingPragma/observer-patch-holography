"""Isolated data and material clock patches on the same seam federation."""
from collections import Counter
from copy import copy
from fractions import Fraction

import numpy as np

from model.candidate_laws.interval_echo import clock
from model.executable_model import carrier
from .gate import require


def patches(level):
    host = carrier.federation(level)
    def neighbors(cell):
        return sorted({p // 12 for a, b in host.external for p, q in ((a,b),(b,a)) if q // 12 == cell})
    front = neighbors(0)[0]
    bank = next(c for c in neighbors(front) if c != 0)
    ports = [p for p in range(host.n) if p // 12 not in (front,bank)]
    local = {physical:i for i,physical in enumerate(ports)}
    data = carrier.Carrier("protected data patch", len(ports),
                           [(local[a],local[b]) for a,b in host.seams if a in local and b in local])
    require(ports[:2] == [0,1] and (0,1) in data.seams, "home fuel must stay on its physical seam")
    return host, data, ports, front, bank


class Branch:
    """Select each positive-probability sharp outcome for exact finite tree traversal."""
    def __init__(self, bit):
        self.u = 0.0 if bit else np.nextafter(1.0, 0.0)

    def random(self):
        return self.u


def clock_laws(host, front, bank, spec):
    owner = {front,bank}
    digits = [tuple(12*bank+i for i in layout) for layout in clock.DIGITS]
    public = tuple(12*bank+i for i in clock.EXTRA_DIGITS[0])
    source = tuple(12*front+i for i in (3,2,4))
    reg = clock.Registers(host, spec["clock_seed"])
    require(reg.error < 1e-11 and all(abs(rate-.5) < 1e-12 for rate in reg.rates.values()),
            "local primitive is not a fair instrument in both input orientations")
    reg.prepare_front(front)
    reg.initialize_digits([*digits,public])
    def tick(register):
        bit = register.fair(front, owner)
        delivered = register.read_tick(source, owner, at_hand=bank)
        require(bit == delivered, "native route lost the source bit")
        register.publish_tick(public, delivered, owner)
        posted = register.post_tick(public, front, digits, owner)
        require(posted == delivered, "hand disagrees with published bit")
        return posted, register.reading(front, digits, owner)

    cap = spec["alarm_after_periods"]
    last = max(spec["thresholds"])
    require(0 < last < 2**len(digits), "clock threshold exceeds its memory")
    witness = []
    sampled = copy(reg)
    sampled.loads = reg.loads.copy()
    sampled.gates = Counter()
    for _ in range(cap):
        bit, count = tick(sampled)
        witness.append([bit,count])
        if count == last:
            break

    # Each distinct rail configuration is a sufficient state. Preserve the weight of every prefix.
    live = {reg.loads.tobytes():(reg, 0, Fraction(1))}
    hits = [[Fraction(0) for _ in range(cap+1)] for _ in range(last+1)]
    hits[0][0] = Fraction(1)
    max_cost = Counter()
    for period in range(1,cap+1):
        following = {}
        for prior, count, weight in live.values():
            for bit in (0,1):
                nxt = copy(prior)
                nxt.loads = prior.loads.copy()
                nxt.gates = Counter()
                nxt.rng = Branch(bit)
                posted, observed = tick(nxt)
                require(posted == bit and observed == count+bit, "bad material clock transition")
                max_cost |= nxt.gates
                if posted:
                    hits[observed][period] += weight/2
                if observed < last:
                    key = nxt.loads.tobytes()
                    if key in following:
                        _, previous, prob = following[key]
                        require(previous == observed, "aliased material counts")
                        following[key] = nxt, observed, prob + weight/2
                    else:
                        following[key] = nxt, observed, weight/2
        live = following
    return np.array(hits,dtype=float), dict(front=front, bank=bank, rail_bits=3,
        sampled_bits_and_counts=witness, sampled_gates=dict(sampled.gates),
        maximum_gates_per_opportunity=dict(max_cost), owner_cells=sorted(owner))
