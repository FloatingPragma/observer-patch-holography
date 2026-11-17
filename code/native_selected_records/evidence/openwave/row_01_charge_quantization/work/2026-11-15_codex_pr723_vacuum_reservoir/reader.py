"""Observer: only graph settings, initial vacuum datum and current tapes enter."""
from model import tree


def terminal_relative(n, edges, tape, cap):
    change, relative, cursor = [0]*n, [None]*n, 0
    relative[0] = 0
    for a, b in tree(n, edges):
        correction = change[a] - change[b]
        first_coin = first_current = None
        for _ in range(cap):
            if cursor == len(tape):
                return None
            aa, bb, coin, current = tape[cursor]
            cursor += 1
            assert (aa, bb) == (a, b) and coin in (0, 1)
            assert type(current) is int
            change[a] -= current
            change[b] += current
            if first_coin is None:
                first_coin, first_current = coin, current
            elif coin != first_coin:
                relative[b] = relative[a] - (2*first_current + current - correction)
                break
        else:
            assert cursor == len(tape)
            return None
    assert cursor == len(tape)
    return relative


def observe(n, edges, body, probe, cap, burn, calibration, reference=12):
    final = terminal_relative(n, edges, probe, cap)
    if final is None:
        return {'outcome': 'inconclusive'}
    initial = list(final)
    for a, b, coin, current in reversed(body):
        initial[a] += current
        initial[b] -= current
    # The only absolute datum is INITIAL load at the declared vacuum port = 0.
    offset = -initial[reference]
    initial = [v + offset for v in initial]
    state = list(initial)
    first_units, post_nonunit, expected, observed = None, 0, 0, 0
    calibration_unit, magnitudes = None, set()
    for t in range(len(body)+1):
        low, high = min(state), max(state)
        units = (0 <= low <= high <= 1) or (-1 <= low <= high <= 0)
        if units and first_units is None:
            first_units = t
        if t >= burn:
            post_nonunit += not units
        if t == len(body):
            break
        a, b, coin, current = body[t]
        if burn <= t < burn + calibration and current and calibration_unit is None:
            calibration_unit = abs(current)
        if t >= burn + calibration:
            magnitudes.add(abs(current))
            observed += current != 0
            expected += sum(state[u] != state[v] for u, v in edges)
        state[a] -= current
        state[b] += current
    return {'outcome': 'conclusive', 'initial_from_records': initial,
            'final_from_records': state, 'total_from_records': sum(initial),
            'first_unit_state': first_units, 'postburn_nonunit_states': post_nonunit,
            'calibration_unit': calibration_unit, 'heldout_magnitudes': sorted(magnitudes),
            'heldout_nonzero_events': observed, 'heldout_cut_sum': expected,
            'conditional_probability_denominator': 2*len(edges)}
