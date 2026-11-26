"""Exact public registers, acceptance supports and conflict-component assembly."""

from dataclasses import dataclass
from itertools import combinations, product

PAIRS = tuple(combinations(range(3), 2))
NAMES = ("b", "g", "x", "y", "z", "p", "s")
DATA = NAMES[2:5]
READ = frozenset(NAMES)
EVENTS = (1, 2, 3, 4)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def mismatch(u):
    return sum(bit != u[0] for bit in u[2:5])


def normal(u):
    return not mismatch(u) and u[5:] == (1, 1)


def rank(u):
    return (mismatch(u), 1 - u[6])


def enabled(u, pair):
    return not u[5] and not u[6] and any(u[2 + i] != u[0] for i in pair)


def seal_ready(u):
    return u[5] and not u[6] and not mismatch(u)


def writer_vector(side, u):
    # The protected boundary witness has its own external stress-test label 5.
    return (
        (5 if u[1] else 0),
        *((1 + side if u[5] else 0),) * 3,
        (3 + side if u[6] else 0),
    )


def parents(side, u):
    return sorted(set(writer_vector(side, u)) - {0})


def functionals(u):
    b, gen, *_, pub, seal = u
    return {
        **{f"mismatch_{i}": int(u[2 + i] != b) for i in range(3)},
        **{f"pair_{i}{j}": int(enabled(u, (i, j))) for i, j in PAIRS},
        "seal_ready": int(bool(seal_ready(u))),
        "unsealed": 1 - seal,
        "boundary": (b, gen),
        "writers": (gen, pub, seal),
        **{f"checkpoint_{k}": (u[2 + k], pub) for k in range(3)},
        "seal_checkpoint": (seal, pub),
    }


def acceptance_supports():
    """Derive exact variable dependence of each registered Boolean functional."""
    inputs = list(product((0, 1), repeat=len(NAMES)))
    table = {u: functionals(u) for u in inputs}
    out = {key: set() for key in table[inputs[0]]}
    for u, row in table.items():
        for i, name in enumerate(NAMES):
            v = (*u[:i], 1 - u[i], *u[i + 1 :])
            for key, val in row.items():
                if val != table[v][key]:
                    out[key].add(name)
    return {k: sorted(v) for k, v in out.items()}


@dataclass(frozen=True)
class Primitive:
    side: int
    tag: str
    writes: frozenset
    payload: tuple

    @property
    def reads(self):
        return READ


def proposals(side, u):
    items = []
    for i, j in PAIRS:
        if enabled(u, (i, j)):
            writes = (DATA[i], DATA[j])
            items.append(
                Primitive(
                    side, f"{i}{j}", frozenset(writes), tuple((k, u[0]) for k in writes)
                )
            )
    if seal_ready(u):
        items.append(Primitive(side, "seal", frozenset({"s"}), (("s", 1),)))
    return items


def conflict(a, b):
    return a.side == b.side and bool(
        a.writes & (b.reads | b.writes) or b.writes & (a.reads | a.writes)
    )


def groups(items):
    groups = []
    for item in items:
        component = [item]
        for g in groups[:]:
            if any(conflict(item, p) for p in g):
                component.extend(g)
                groups.remove(g)
        groups.append(component)
    return groups


@dataclass(frozen=True)
class Prepared:
    side: int
    members: tuple
    snapshot: tuple
    payload: tuple
    writes: frozenset

    @property
    def reads(self):
        return READ

    @property
    def is_seal(self):
        return self.members == ("seal",)

    @property
    def identifier(self):
        return self.side + (3 if self.is_seal else 1)

    def result(self):
        values = dict(zip(NAMES, self.snapshot, strict=True))
        values.update(self.payload)
        return tuple(values[k] for k in NAMES)


def assemble(component, u):
    side = component[0].side
    words = {}
    for p in component:
        require(p.side == side, "crossed independent collars")
        for name, value in p.payload:
            require(name not in words or words[name] == value, "incompatible overlap")
            words[name] = value
    if component[0].tag != "seal":
        # Each bad address belongs to two pairs; their union covers the triangle.
        require(set(words) == set(DATA), "incomplete union")
        words["p"] = 1
    return Prepared(
        side,
        tuple(sorted(p.tag for p in component)),
        u,
        tuple(sorted(words.items())),
        frozenset(words),
    )


def batch(states):
    all_proposals = [p for i, u in enumerate(states) for p in proposals(i, u)]
    primitive_groups = groups(all_proposals)
    transactions = [assemble(g, states[g[0].side]) for g in primitive_groups]
    # Closure adds the publication selector to data writes. Recheck conflicts
    # on these enlarged supports, folding their original members if necessary.
    while True:
        recomputed = groups(transactions)
        if all(len(g) == 1 for g in recomputed):
            return transactions
        new_groups = []
        for g in recomputed:
            members = [
                p
                for c in g
                for p in all_proposals
                if p.side == c.side and p.tag in c.members
            ]
            new_groups.append(assemble(members, states[g[0].side]))
        transactions = new_groups


def validate(tx, current):
    if current != tx.snapshot:
        return "stale"
    members = proposals(tx.side, current)
    if not members or assemble(members, current) != tx:
        return "component"
    after = tx.result()
    if after[:2] != current[:2] or rank(after) >= rank(current):
        return "invariant"
    return "ok"


def check_genesis(states, ids=EVENTS):
    require(len(set(ids)) == 4 and tuple(ids) == EVENTS, "event namespace")
    for i, u in enumerate(states):
        require(not (set(writer_vector(i, u)) & set(ids)), "genesis collision")
        require(u[1] == 0 and mismatch(u) > 0, "inadmissible genesis")


def launch(a, b, x, y):
    return [(a, 0, *x, 0, 0), (b, 0, *y, 0, 0)]
