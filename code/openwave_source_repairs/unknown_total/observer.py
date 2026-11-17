"""Record-only census and state transport; no total or hidden state argument."""


def transport(state, record):
    a, b, coin, current = record
    out = list(state)
    out[a] -= current
    out[b] += current
    return out


def decode_census(tape, tree, nodes, vacuum_datum, cap=64):
    delta = {v: 0 for v in nodes}
    relative = {tree[0][0]: 0}
    cursor = 0
    for a, b in tree:
        correction = delta[a] - delta[b]
        first = None
        for _ in range(cap):
            if cursor == len(tape):
                return None
            aa, bb, coin, current = tape[cursor]
            cursor += 1
            if (aa, bb) != (a, b):
                raise ValueError('Wrong addressed census seam')
            delta[a] -= current
            delta[b] += current
            if first is None:
                first, j0 = coin, current
            elif coin != first:
                contrast = 2*j0 + current - correction
                relative[b] = relative[a] - contrast
                break
        else:
            return None
    if cursor != len(tape):
        raise ValueError('Unused census records')
    offset = -relative[vacuum_datum]
    initial = [relative[v] + offset for v in nodes]
    terminal = [relative[v] + offset + delta[v] for v in nodes]
    return dict(initial=initial, terminal=terminal, total=sum(initial))
