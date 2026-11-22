"""Reconstruct the signed records of a cell from its three addressed boundary bits."""

from collections import Counter

from .closure import certificate
from .equipment import addressed_bits
from .observed import need


def boundary_replay(chart, inner, ports, flags):
    # The receiver reads the internal loads, chart and three boundary outcomes.
    loads = inner.copy()
    for u, hit in zip(ports, flags):
        if hit:
            loads[u] = 1
    marks = []
    for a, b in chart:
        sign = loads[a] - loads[b]
        hit = sign if abs(sign) == 2 else 0
        marks.append(hit)
        if hit:
            loads[a] = loads[b] = 1
    return marks, loads


def cell_inputs(screen, order):
    inner = [[] for _ in range(screen.n // 12)]
    external = [[] for _ in inner]
    for i in order:
        a, b = screen.edges[i]
        if a // 12 == b // 12:
            inner[a // 12].append(i)
        else:
            external[a // 12].append([a, b, i])
            external[b // 12].append([b, a, i])
    return [(group, sorted(bound)) for group, bound in zip(inner, external)]


def run_certificate(screen, service, order, start):
    # This interpreter touches only the addresses in the accepted certificate.
    # Each predecessor completes first, so the patch can use the local ready
    # instrument without asking a nonmember about its state.
    loads = start.copy()
    marks = [0] * len(screen.edges)
    for e in order:
        if service.ready(e, loads, marks):
            service.commit(e, loads, marks)
    return loads, marks


def trial(screen, service, rel, cert, local, chart, boundary, spare, mask):
    ports = [u % 12 for u, _, _ in boundary]
    inner = [2 if j in ports else 0 for j in range(12)]
    offset = boundary[0][0] // 12 * 12
    loads = [1] * screen.n
    cell_ports = slice(offset, offset + 12)
    loads[cell_ports] = inner
    flags = [(mask >> j) & 1 for j in (2, 1, 0)]
    for (_, v, _), w, bit in zip(boundary, spare, flags):
        loads[v] = 1 - bit
        loads[w] = bit
    need(Counter(loads) == {0: 12, 1: screen.n - 15, 2: 3}, "global histogram")
    cut_loads, cut_marks = run_certificate(screen, service, cert["order"], loads)
    if mask == 7:
        # Omitting one accepted boundary write changes an internal readout.
        omitted = boundary[0][2]
        _, altered = run_certificate(screen, service, [i for i in cert["order"] if i != omitted], loads)
        need(any(altered[i] != cut_marks[i] for i in local), "missing boundary did not alter internal readout")
    marks = [0] * len(screen.edges)
    # The full service attempts every address in the order recovered from the hits.
    for e in rel.order:
        if service.enabled(e, loads, marks):
            service.commit(e, loads, marks)
    need(all(cut_marks[i] == marks[i] for i in cert["order"]) and
         cut_loads[cell_ports] == loads[cell_ports], "certified restricted run")
    observed = [addressed_bits()[cut_marks[e]] for _, _, e in boundary]
    need(observed == flags, "boundary environment not realized")
    actual = [marks[i] for i in local]
    replayed, ends = boundary_replay(chart, inner, ports, observed)
    need(actual == replayed and ends == loads[cell_ports], "three-flag local replay")
    word = [j for j, mark in enumerate(actual) if mark]
    covered = {v for j in word for v in chart[j]}
    need([int(u not in covered) for u in ports] == flags, "flags not decoded from internal hits")
    return word, actual, ends


def read_cell(screen, service, relation, board, cell, local, boundary):
    """Run each preparation on the certified patch and on the full service."""
    seeds = set(local)
    cert = certificate(relation, board, seeds)
    need(len(seeds) == 30 and len(boundary) == 3 and
         set(cert["order"]) == seeds | {i for _, _, i in boundary}, "expected a 30+3 closure")
    # Define the readout chart in the certified order consumed by the patch.
    local = [address for address in cert["order"] if address in seeds]
    own = set(range(12 * cell, 12 * cell + 12))
    remote = {v for _, v, _ in boundary}
    ports = [u % 12 for u, _, _ in boundary]
    need(len(remote) == len(set(ports)) == 3 and not own & remote, "overlapping boundary")
    spare = [v for v in range(screen.n) if v not in own | remote][:3]
    chart = [(screen.edges[i][0] % 12, screen.edges[i][1] % 12) for i in local]
    witnesses = []
    for mask in range(8):
        word, signed, ends = trial(screen, service, relation, cert, local, chart, boundary, spare, mask)
        witnesses.append(dict(mask=mask, boundary_hits=[(mask >> j) & 1 for j in (2, 1, 0)],
                              internal_word=word, signed_marks=signed, final_loads=ends))
    words = [row["internal_word"] for row in witnesses]
    need(len(set(map(tuple, words))) == 8, "lower-bound collision")
    return dict(certificate=cert, boundary_addresses=[i for _, _, i in boundary],
                boundary=boundary, compensation=spare, internal_addresses=local,
                ports=ports, chart=chart, classical_bits=3, witnesses=witnesses)
