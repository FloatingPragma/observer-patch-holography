"""Two physical audits select and retire an offer before its fresh noise draw."""

from functools import reduce
from operator import or_, xor
from itertools import product

from .....seam_word_meter.finite_program.machine import template_matching
from ..instrument import Bench, BOUNDARY, COMMIT, FAR, NOMINAL, packed

PAIRS = tuple(
    (12 * o + a, 12 * o + b)
    for o in (18, 19)
    for a, b in template_matching(set(range(12)))
)
SHELVES = PAIRS[:9]


def encode(columns):
    return [format(v, "x") for v in columns]


def choice(remaining, odd, full):
    """Partitions the measured copies by their selected material label."""
    groups = [0] * 8
    for allowed, order in ((full ^ odd, range(8)), (odd, reversed(range(8)))):
        for label in order:
            groups[label] |= allowed & remaining[label]
            allowed &= full ^ remaining[label]
    return groups


def start(leaves, flags, mask=255, scanner=None):
    words = [
        [b ^ ((flags >> (3 * s + i)) & 1) for i in range(3)]
        for s, b in enumerate(BOUNDARY)
    ]
    bench = Bench(leaves, words)
    assert len({p for edge in SHELVES for p in edge}) == 18
    for i, edge in enumerate(SHELVES):
        assert edge[1] in bench.device.car.neighbours[edge[0]]
        bench.set_pair(edge, (mask >> i) & 1)
    if scanner is not None:
        for edge in PAIRS[9:]:
            bench.set_pair(edge, 0)
    return bench


def shelf(bench):
    def read(at):
        return [bench.read_bit(edge, at, bench.full) for edge in SHELVES]

    first, second = read(COMMIT), read(FAR)
    if first != second:
        raise ArithmeticError("Two shelf readers disagree")
    return first


def flags(bench):
    def read(at):
        observed = []
        for side in (0, 1):
            boundary, word, _ = bench.state(side, at)
            observed.extend([v ^ boundary for v in word])
        return observed

    first, second = read(COMMIT), read(FAR)
    if first != second:
        raise ArithmeticError("Two flag readers disagree")
    return first


def proposal(bench, parked=False, scanner=None):
    # Quiescent empty stock has no active group. Nominal runs exclude it.
    at = len(bench.reads)
    mask = shelf(bench)
    state = flags(bench)
    extra = {}
    if scanner is None:
        groups = choice(mask[:8], reduce(xor, state), bench.full)
        for edge, group in zip(SHELVES, groups):
            if not parked:
                bench.swap([edge], group)
    else:
        groups, extra = scanner(bench, state, parked)
    after = shelf(bench)
    expected = [v ^ x for v, x in zip(mask, groups)] + [mask[8]]
    return (
        {
            **extra,
            "before": encode(mask),
            "flags": encode(state),
            "groups": encode(groups),
            "retired": encode(after),
            "read_span": [at, len(bench.reads)],
        },
        groups,
        after,
        after == expected,
    )


def grouped_offer(bench, groups, mode):
    # Each pointwise native call gets only the actually selected copies.
    original = bench.loads
    output = original.copy()
    all_branches = bench.full
    for label, selected in enumerate(groups):
        if not selected:
            continue
        bench.full = selected
        bench.loads = [v & selected for v in original]
        bench.offer(*NOMINAL[label], tuple(v & selected for v in mode))
        output = [(old & ~selected) | new for old, new in zip(output, bench.loads)]
    bench.loads = output
    bench.full = all_branches


def step(bench, k, park=False, scanner=None):
    receipt, groups, after, valid = proposal(bench, parked=park, scanner=scanner)
    if not valid or reduce(or_, groups) != bench.full:
        receipt["status"] = "refused"
        return receipt
    bench.draw(k)
    mode = bench.latched()
    grouped_offer(bench, groups, mode)
    receipt["fault"] = encode(mode)
    receipt["state"] = bench.snapshot()
    wrap = bench.full ^ reduce(or_, after[:8])
    if wrap:
        for edge in SHELVES:
            bench.swap([edge], wrap)
    final = shelf(bench)
    if final != [v ^ wrap for v in after]:
        raise ArithmeticError("Mask refill or wrap marker failed")
    receipt.update(final=encode(final), wrap=format(wrap, "x"), status="committed")
    return receipt


def run(leaves, q, park=False, scanner=None):
    bench = start(leaves, q, scanner=scanner)
    initial = bench.snapshot()
    steps = []
    for k in range(len(leaves[0])):
        record = step(bench, k, park=park, scanner=scanner)
        steps.append(record)
        if record["status"] == "refused":
            break
    return dict(
        initial=initial,
        steps=steps,
        ending=bench.snapshot(),
        reads=bench.reads,
        leaves=len(leaves),
    )


def mask_table(q, scanner=None):
    # Each bit column represents all 256 sharp initial register words.
    bench = start([(0,)] * 256, q, scanner=scanner)
    for i, edge in enumerate(SHELVES[:8]):
        bits = packed([(m >> i) & 1 for m in range(256)])
        bench.swap([edge], bench.full ^ bits)
    obs, groups, _, valid = proposal(bench, scanner=scanner)
    assert valid and reduce(or_, groups) == bench.full ^ 1
    return dict(proposal=obs, reads=bench.reads)


def report():
    return dict(
        shelves=SHELVES,
        tree=run(list(product(range(3), repeat=8)), 8),
        two_blocks=run([tuple(i % 3 for i in range(16))], 8),
        parked=run([(0,)], 8, park=True),
        mask_probes=[mask_table(q) for q in (0, 8)],
    )
