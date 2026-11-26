"""Declared sharp programs. Every decision, copy and capture consumes a pointer read."""

from .machine import CENTER, SIZE


def repair(bench, left, right):
    result = bench.mean(left, right)
    if result != "accepted":
        raise ArithmeticError(f"Unrepresentable or unreadable repair: {result}")


def accepted(possible, seen):
    candidates = possible.get(seen, set())
    return next(iter(candidates)) if len(candidates) == 1 else None


def meaning(code):
    if code == CENTER:
        return None
    return "low" if code < CENTER else "high"


def read_family(bench, registers):
    return [[bench.read(r, station) for r in registers] for station in (0, 1)]


def adaptive(bench, first, second, mode="native"):
    bench.start([first, second, 32, 96] + [0] * (SIZE - 4))
    repair(bench, 0, 1)
    observation = (lambda v: abs(v - CENTER)) if mode == "fold" else (lambda v: v)
    possible = {}
    for code in (48, 56, 64, 72, 80):
        possible.setdefault(observation(code), set()).add(meaning(code))
    read = bench.read(0, 0)
    decision = accepted(possible, observation(read))
    route = "low" if mode == "force_low" and decision else decision
    if route:
        repair(bench, *((0, 2) if route == "low" else (1, 3)))
    return dict(
        inputs=[first, second],
        mode=mode,
        observed=observation(read),
        chosen=route,
        public=bench.public(route),
        final=read_family(bench, range(4)),
        tape=bench.tape,
    )


STOCK = ((16, 0), (0, 0), (0, 16))


def real_menu(bench, first, second, mode="native"):
    bench.start([*first, *second, *STOCK[2], *STOCK[0]] + [0] * (SIZE - 8))
    for k in (0, 1):
        repair(bench, k, 2 + k)
    u, v = (bench.read(k, 0) for k in (0, 1))
    diff = v * v - (1 if mode == "equal_weights" else 2) * u * u
    route = None if diff == 0 else ("low" if diff < 0 else "high")
    if route:
        source, target = (0, 4) if route == "low" else (2, 6)
        for k in (0, 1):
            repair(bench, source + k, target + k)
    return dict(
        inputs=[first, second],
        mode=mode,
        received=[u, v],
        chosen=route,
        public=bench.public(route),
        final=read_family(bench, range(8)),
        tape=bench.tape,
    )


def balanced(bench, offsets):
    initial = [v for a in offsets for v in (CENTER + a, CENTER - a)]
    bench.start(initial + [CENTER] * (SIZE - len(initial)))


def compiled(bench, instructions):
    for op, i, *others in instructions:
        if op == "clear":
            repair(bench, 2 * i, 2 * i + 1)
        elif op == "copy":
            (j,) = others
            for k in (0, 1):
                repair(bench, 2 * i + k, 2 * j + k)
        else:
            raise ValueError(op)


def bus(bench):
    balanced(bench, (24, -8, 8, 16, -24, 8))
    code = [("clear", 2), ("copy", 2, 3), ("copy", 3, 4), ("copy", 4, 5)]
    compiled(bench, code)
    return dict(
        instructions=code, final=read_family(bench, range(SIZE)), tape=bench.tape
    )


def accumulator(bench):
    balanced(bench, (16, 0, -8))
    # addWord 2 at zero initial scales: shrink accumulator once, copy and clear.
    code = [("copy", 0, 1), ("clear", 1), ("copy", 2, 1), ("copy", 0, 1), ("clear", 1)]
    compiled(bench, code)
    return dict(
        instructions=code,
        scales=[2, 0, 1],
        final=read_family(bench, range(6)),
        tape=bench.tape,
    )


def hop(bench, source, skip_reset=False):
    bench.start([source, 22, 12, 34] + [0] * (SIZE - 4))
    intermediate = []
    bench.put(1, bench.read(0, 0), 0)
    intermediate.append(read_family(bench, range(4)))
    if not skip_reset:
        bench.put(2, 0, 1)
    intermediate.append(read_family(bench, range(4)))
    repair(bench, 1, 2)
    intermediate.append(read_family(bench, range(4)))
    bench.put(3, 2 * bench.read(2, 1), 1)
    intermediate.append(read_family(bench, range(4)))
    bench.put(1, 0, 0)
    intermediate.append(read_family(bench, range(4)))
    bench.put(2, 0, 1)
    intermediate.append(read_family(bench, range(4)))
    return dict(
        source=source, skip_reset=skip_reset, stages=intermediate, tape=bench.tape
    )


def refusal(bench, a, b, control):
    bench.start([a, b] + [0] * (SIZE - 2))
    initial = bench.loads.copy()
    if control == "parked_sensor":
        bench.sensor = "parked"
    result = bench.mean(0, 1, parked=control == "parked_write")
    unchanged = bench.loads == initial
    bench.sensor = "native"
    return dict(
        inputs=[a, b],
        control=control,
        outcome=result,
        preserved=unchanged,
        final=read_family(bench, (0, 1)),
        tape=bench.tape,
    )
