"""Record-only admission. No seed, case label, filename, initial state or total input."""

def tree(n, edges):
    seen = {0}
    queue = [0]
    out = []
    for a in queue:
        for u, v in edges:
            b = v if a == u else u if a == v else None
            if b is not None and b not in seen:
                seen.add(b)
                queue.append(b)
                out.append((a, b))
    if len(seen) != n:
        raise ValueError('disconnected core')
    return out

def decode(n, edges, records, cap=64):
    changes = [0] * n
    relative = [None] * n
    relative[0] = 0
    cursor = 0
    for a, b in tree(n, edges):
        correction = changes[a] - changes[b]
        first = None
        for _ in range(cap):
            if cursor == len(records):
                return None
            aa, bb, c, j = records[cursor]
            cursor += 1
            if (aa, bb) != (a, b) or c not in (0, 1) or type(j) is not int:
                raise ValueError('bad event')
            changes[a] -= j
            changes[b] += j
            if first is None:
                first = (c, j)
            elif c != first[0]:
                relative[b] = relative[a] - (2 * first[1] + j - correction)
                break
        else:
            return None
    if cursor != len(records):
        raise ValueError('unused records')
    offset = -relative[12]
    initial = [v + offset for v in relative]
    return (initial, [v + d for v, d in zip(initial, changes)])

def path_to(start, predicate, edges):
    queue = [start]
    previous = {start: None}
    end = None
    for a in queue:
        if predicate(a):
            end = a
            break
        for u, v in edges:
            b = v if a == u else u if a == v else None
            if b is not None and b not in previous:
                previous[b] = a
                queue.append(b)
    if end is None:
        raise ValueError('unreachable destination')
    path = [end]
    while previous[path[-1]] is not None:
        path.append(previous[path[-1]])
    return path[::-1]

def admit(request, g, cap=64):
    n = g['core_N']
    probe = []
    for edge in tree(n, g['core']):
        first = None
        for _ in range(cap):
            r = request(edge, 'probe')
            probe.append(r)
            if first is None:
                first = r[2]
            elif r[2] != first:
                break
        else:
            return {'outcome': 'probe_timeout', 'probe': probe}
    initial, state = decode(n, g['core'], probe, cap)
    q = state + [0] * (g['N'] - n)
    for ring in g['rings']:
        q[ring[2 * (len(ring) // 3)]] = 1
    result = {'probe': probe, 'initial_core': initial, 'post_probe_core': state, 'balance': [], 'route': [], 'strict_decreases': []}

    def take(edge, stage):
        r = request(edge, stage)
        a, b, c, j = r
        if (a, b) != tuple(edge) or c not in (0, 1) or type(j) is not int:
            raise ValueError('bad event')
        q[a] -= j
        q[b] += j
        result[stage].append(r)
        return r

    def exchange(a, b, stage):
        wanted = (q[b], q[a])
        for _ in range(cap):
            take((a, b), stage)
            if (q[a], q[b]) == wanted:
                return True
        return False
    while max(q[:n]) - min(q[:n]) > 1:
        maximum = max(q[:n])
        start = q.index(maximum)
        path = path_to(start, lambda v: q[v] <= maximum - 2, g['core'])
        for a, b in zip(path[:-2], path[1:-1]):
            if q[b] == maximum:
                continue
            if q[a] != maximum or q[b] != maximum - 1:
                raise ValueError('bad balancing path')
            if not exchange(a, b, 'balance'):
                return result | {'outcome': 'balance_timeout', 'observed_final': q}
        before = sum((v * v for v in q[:n]))
        take(tuple(path[-2:]), 'balance')
        after = sum((v * v for v in q[:n]))
        result['strict_decreases'].append([before, after])
        if after >= before:
            return result | {'outcome': 'repair_failed', 'observed_final': q}
    result['settled_core'] = q[:n]
    if not any(q[:n]):
        return result | {'outcome': 'no_unit', 'observed_final': q}
    if not (set(q[:n]) <= {0, 1} or set(q[:n]) <= {-1, 0}):
        return result | {'outcome': 'capacity_boundary', 'observed_final': q}
    target = g['target']
    sign = 1 if max(q[:n]) == 1 else -1
    path = path_to(target[1], lambda v: q[v] != 0, g['core'] + [(0, target[1])])[::-1]
    for a, b in zip(path, path[1:]):
        if q[a] != sign or q[b] != 0:
            raise ValueError('occupied extraction corridor')
        if not exchange(a, b, 'route'):
            return result | {'outcome': 'route_timeout', 'observed_final': q}
    if [q[v] for v in target] != [0, sign, 0]:
        raise ValueError('target not isolated')
    return result | {'outcome': 'signed_target_unit', 'sign': sign, 'observed_final': q, 'inferred_initial_total': sum(initial), 'retained_core_total': sum(q[:n])}

def require(condition, message='scientific check failed'):
    if not condition:
        raise AssertionError(message)
