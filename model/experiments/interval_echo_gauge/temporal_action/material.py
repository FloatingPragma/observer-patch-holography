"""Chart incidence, public charge census and independent seam crossings."""

from functools import lru_cache
import numpy as np
import networkx as nx
from .rails import hardware, read_reference
from model.executable_model import federation
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.face_frame import (
    shared_report,
)
from scipy.linalg import expm

from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path
from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import CELL_SEAMS, isomorphism


@lru_cache(maxsize=1)
def geometry():
    path = (
        Path(__file__).resolve().parents[4]
        / "code/maxwell_local_selection/certificate.py"
    )
    spec = spec_from_file_location("maxwell_local_selection_certificate", path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    d, c, _ = module.carrier()
    d, c = np.array(d, dtype=int), np.array(c, dtype=int)
    seams = [
        (int(np.flatnonzero(row == -1)[0]), int(np.flatnonzero(row == 1)[0]))
        for row in d
    ]
    pin = isomorphism(seams, list(CELL_SEAMS))
    if pin is None:
        raise ValueError("Face chart is not a native cell incidence")
    vertices = set(v for e in np.flatnonzero(c[0]) for v in seams[e])
    first = min(vertices)
    next_vertex = next(
        b if a == first else a
        for e, (a, b) in enumerate(seams)
        if (a == first and c[0, e] == 1) or (b == first and c[0, e] == -1)
    )
    last = next(iter(vertices - {first, next_vertex}))
    hole = min(
        v for a, b in seams for v in (a, b) if first in (a, b) and v not in vertices
    )
    return d, c, seams, list(map(int, pin)), [first, next_vertex, last], hole


@lru_cache(maxsize=1)
def reader():
    board = federation(0)
    graph = nx.Graph(board.seams)
    held = {*hardware(), read_reference()}
    local = set(range(12, 24)) - held
    contact = next(
        p for p in sorted(local) if len(set(board.neighbours[p]) & local) >= 2
    )
    standards = sorted(set(board.neighbours[contact]) & local)[:2]
    net = graph.subgraph(set(graph) - held - set(standards))
    paths = [nx.shortest_path(net, i, contact) for i in geometry()[3]]
    for i, path in zip(geometry()[3], paths):
        pairs = list(zip(path, path[1:]))
        labels = np.arange(board.n)
        for a, b in pairs + pairs[-2::-1]:
            labels[a], labels[b] = labels[b], labels[a]
        expected = np.arange(board.n)
        expected[i], expected[contact] = expected[contact], expected[i]
        np.testing.assert_array_equal(labels, expected)
    return dict(
        contact=contact, reference_0=standards[0], reference_2=standards[1], paths=paths
    )


@lru_cache(maxsize=1)
def instrument():
    lookup = shared_report.certify()
    return np.array(
        [[int(lookup[value, ref, False]) for ref in (0, 2)] for value in range(3)]
    )


def readings(stock):
    raw = instrument()[stock]
    # "01" denotes load zero, "00" one, "10" two; "11" refuses.
    if np.any(raw.sum(axis=1) > 1):
        raise ValueError("Refused ternary census")
    return raw, (1 + raw[:, 0] - raw[:, 1])


@lru_cache(maxsize=None)
def source(case):
    d, c, seams, pin, (u, v, w), hole = geometry()
    paths = {
        "loop": [[v, w, u], [], [], []],
        "reverse": [[w, v, u], [], [], []],
        "idle": [[], [], [], []],
        "open_quiet": [[u], [w], [], []],
        "open_swap": [[v], [w], [], []],
    }
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + 8 * contrast
    initial = np.eye(law.DIM)[:, law.SRC.index((2, 1))]
    endpoints = {}
    for kind, t in (("half", np.pi / 2), ("full", np.pi)):
        col = expm(-1j * t * h) @ initial
        expected = np.zeros(len(col))
        expected[law.SRC.index((2, 1))] = 0.5 if kind == "half" else 0
        expected[law.SRC.index((1, 2))] = 0.5 if kind == "half" else 1
        np.testing.assert_allclose(abs(col) ** 2, expected, atol=1e-12)
        endpoints[kind] = [float(expected[law.SRC.index(a)]) for a in ((2, 1), (1, 2))]
    here = u
    stock = np.ones(12, dtype=int)
    stock[hole], stock[u] = 0, 2
    raw, decoded = readings(stock)
    censes = [decoded.tolist()]
    raw_census = [raw.tolist()]
    histories = []
    currents = []
    for n, legs in enumerate(paths[case]):
        current = np.zeros(30, dtype=int)
        history = []
        for destination in legs:
            before = here
            target = v if n == 0 and case.startswith("open_") else destination
            edge = seams.index(tuple(sorted((before, target))))
            event = "half" if n == 0 and case.startswith("open_") else "full"
            observed = 1 if destination != here else 0
            if endpoints[event][observed] <= 0:
                raise ValueError("Impossible sharp return")
            old = stock[[before, target]].tolist()
            if observed:
                stock[before], stock[target] = stock[target], stock[before]
                here = target
            old_raw, old_dec = readings(np.array(old))
            new_raw, new_dec = readings(stock[[before, target]])
            # Current comes from the ordered crossing at the visited contact.
            direction = 1 if seams[edge] == (before, target) else -1
            current[edge] += 2 * direction * int(new_dec[0] - old_dec[0])
            history.append(
                dict(
                    raw_before=old_raw.tolist(),
                    raw_after=new_raw.tolist(),
                    ports=[pin[before], pin[target]],
                    before=old,
                    after=stock[[before, target]].tolist(),
                    pulse=event,
                    probability=endpoints[event][observed],
                )
            )
        raw, decoded = readings(stock)
        censes.append(decoded.tolist())
        raw_census.append(raw.tolist())
        currents.append(current)
        histories.append(history)
    return (
        np.array(censes) - 1,
        np.array(currents),
        dict(
            material_ports=pin,
            census=censes,
            raw_census=raw_census,
            contacts=histories,
            hole=hole,
            half_probabilities=endpoints["half"],
        ),
    )
