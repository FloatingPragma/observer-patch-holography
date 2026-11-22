"""Test a local precedence rule on native strict-repair quantum instruments."""
import argparse
from itertools import permutations, product
import json

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import federation
from .__main__ import HERE, local_ops
from .asynchronous import SharpService, induced_path, priority, verify_native_service


def service_column(service, order, input_loads, input_collars):
    loads, collars = list(input_loads), list(input_collars)
    record = [""] * len(collars)
    for i in order:
        if collars[i]:
            record[i] = "occupied"
        elif service.ready(i, loads, collars):
            record[i] = "write"
            service.commit(i, loads, collars)
        else:
            record[i] = "quiet"
    if any(service.ready(i, loads, collars) for i in order):
        raise ArithmeticError("a terminal branch contains a discrepancy")
    return tuple(record), tuple(loads), tuple(collars)


def dependence_depth(level):
    board = federation(level)
    sweep = SharpService(range(board.n), board.seams)
    depth = [0] * len(board.seams)
    for i in sorted(range(len(depth)), key=sweep.rank.__getitem__):
        depth[i] = 1 + max((depth[j] for j in sweep.higher[i]), default=0)
    if max(depth) > 31:
        raise ArithmeticError("a dependency chain exceeds the number of distinct priorities")
    return dict(cells=board.cells, seams=len(depth), visit_layers=max(depth),
                maximum_incident_predecessors=max(map(len, sweep.higher)))


def path_orders(service):
    orders = []
    for order in permutations(range(len(service.edges))):
        when = {edge: step for step, edge in enumerate(order)}
        if all(when[j] < when[i] for i in order for j in service.higher[i]):
            orders.append(order)
    return orders


def phase_pair(board, ports, edges, service, order):
    bridge = next(i for i, e in enumerate(edges) if e in board.external)
    outer = next(i for i, e in enumerate(edges)
                 if set(e).isdisjoint(edges[bridge]) and e[0]//12 == edges[bridge][1]//12)
    initial = []
    terminals = []
    for orientation in (0, 1):
        loads = [0]*len(ports)
        for addr in (bridge, outer):
            loads[service.edges[addr][orientation]] = 2
        initial.append(loads)
        terminals.append(service_column(service, order, loads, [0]*len(edges)))
    if terminals[0][0] != terminals[1][0] or set(i for i, op in enumerate(terminals[0][0])
                                                  if op == "write") != {bridge, outer}:
        raise ArithmeticError("the two coherent histories lack the same public transcript")
    premature = next(i for i in order if service.ready(i, initial[0], [0]*len(edges))
                     != service.ready(i, initial[1], [0]*len(edges)))
    echo, _, hit, _ = local_ops()
    h0, h1 = native_law.local_term_parts()
    half = expm(-.5j*np.pi*(h0+KAPPA*h1))
    launch = [native_law.SRC.index(n) for n in ((2, 0), (0, 2))]
    # For a certified hit, the sign/coin state is hit @ |source>. Unwrite it and
    # mix the restored double rails before interrogating the first port.
    read = (half @ echo.T @ hit[:, launch])[launch, :]
    np.testing.assert_allclose(read.conj().T @ read, np.eye(2), atol=1e-12)
    spare = [next(p for p in board.neighbours[edges[i][0]]
                  if p not in ports and p//12 == edges[i][0]//12) for i in (bridge, outer)]
    if spare[0] == spare[1]:
        raise ArithmeticError("the two readout spares coincide")
    probe = echo.reshape(11,3,3,99)[:,1,1,:][:, [native_law.SRC.index((n,0)) for n in range(3)]]
    np.testing.assert_allclose(probe.T @ probe, np.diag([0,0,1]), atol=1e-12)
    two = np.kron(read, read)
    distributions = [abs((two[:,0] + s*two[:,3])/np.sqrt(2))**2 for s in (1,-1)]
    erasure = (abs(two[:,0])**2 + abs(two[:,3])**2)/2
    even = [float(p[0]+p[3]) for p in (*distributions, erasure)]
    np.testing.assert_allclose(even, [0.,1.,.5], atol=1e-12)
    return dict(written_edges=[list(edges[i]) for i in (bridge,outer)], vacancy_ports=spare,
                premature_seam=list(edges[premature]), common_transcript=list(terminals[0][0]),
                plus_joint=np.round(distributions[0],12).tolist(),
                minus_joint=np.round(distributions[1],12).tolist(),
                premature_probe_joint=np.round(erasure,12).tolist())


def run():
    verify_native_service()
    board = federation(0)
    ports = induced_path(board)
    edges = [e for e in board.seams if set(e) <= set(ports)]
    service = SharpService(ports, edges)
    orders = path_orders(service)
    reference = tuple(sorted(range(len(edges)), key=lambda i: priority(edges[i])))
    if reference not in orders or len(orders) < 2:
        raise ArithmeticError("the patch lacks independently scheduled visits")
    words, inputs, control_count = set(), 0, 0
    reverse = reference[::-1]
    for loads in product(range(3), repeat=len(ports)):
        for markers in product((0, 1), repeat=len(edges)):
            # Marker 1 denotes any preoccupied collar; on that entire block
            # the native instrument is identity. New signs take values ±2.
            target = service_column(service, reference, loads, markers)
            for order in orders:
                if service_column(service, order, loads, markers) != target:
                    raise ArithmeticError("schedule changes a positive instrument column")
            words.add(target[0])
            inputs += 1
            if service_column(service, reverse, loads, markers) != target:
                control_count += 1
    if control_count == 0:
        raise ArithmeticError("removing the precedence guards hid the negative control")
    return dict(board_depths=[dependence_depth(level) for level in range(3)],
                path=list(ports), seams=[list(e) for e in edges], legal_orders=len(orders),
                input_flag_and_load_words=inputs, transcript_classes=len(words),
                changed_by_unguarded_reverse=control_count,
                phase=phase_pair(board, ports, edges, service, reference))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="replace the coherent sweep receipt")
    args = parser.parse_args()
    path = HERE/"coherent_sweep.json"
    result = run()
    if args.write:
        path.write_text(json.dumps(result, indent=2)+"\n")
    elif result != json.loads(path.read_text()):
        raise ValueError("the registered coherent sweep outcome differs")
    print(json.dumps(result["phase"], indent=2))


if __name__ == "__main__":
    main()
