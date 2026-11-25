"""Embed the balanced stored-collar ladder beside the actual winding receiver."""

import networkx as nx
import numpy as np

from model.executable_model import federation
from model.experiments.interval_echo_consensus.signed_effects.collar_storage import (
    conserved_columns,
)
from model.experiments.interval_echo_gauge.arrival_ladder import measure as ladder
from ...apparatus import fixture
from ...echo_turn.material import chart, transp
from ..instruments import extract

HUB, TARGET, GUARD, DONOR = 160, 167, 158, 166


def stock(sign):
    f = fixture()
    b = federation(0)
    lab = chart()
    material = nx.Graph(b.seams)
    guide = set(f["ports"])
    outside = set(material) - guide - set(f["donor"])
    graph = nx.Graph(f["blank"])
    graph.add_edges_from(material.subgraph(outside).edges)
    service = lab["service"]
    assert service[1] in outside and service[0] in guide
    special = ((HUB, TARGET), (GUARD, TARGET), (DONOR, TARGET))
    assert all(material.has_edge(*e) for e in special)
    assert all(not graph.has_edge(*e) for e in special)
    assert material.has_edge(HUB, GUARD) and graph.has_edge(HUB, GUARD)
    data_edges = {frozenset(z) for z in f["marked"]}
    assert not data_edges.intersection(frozenset(z) for z in special)
    assert TARGET not in guide and b.degree[TARGET] == 5 and HUB == f["tip"]
    # Two ordinary (2,0) source pairs write the even |d|=2 collars.
    # Two additional balanced pairs provide the live guard and donor stocks.
    slots = lab["storage"] + lab["probes"] + lab["stock"]
    far = [int(x) for x in slots]
    loads = np.ones(b.n, int)
    loads[HUB] = 1 + sign
    loads[f["donor"][0]] = 1 - sign
    labels = [2, 0, 2, 0, 1 - sign, 1 + sign, 1 - sign, 1 + sign]
    loads[far[:8]] = labels
    assert loads.sum() == b.n
    parked = 165
    operations = []

    def trans(u, v):
        route = nx.shortest_path(graph, u, v)
        links = transp(route)
        old = loads.copy()
        for a, z in links:
            loads[a], loads[z] = loads[z], loads[a]
        old[u], old[v] = old[v], old[u]
        np.testing.assert_array_equal(loads, old)
        operations.append(dict(source=u, destination=v, swaps=len(links)))

    trans(HUB, parked)
    for e in (service, *special):
        assert loads[e[0]] == loads[e[1]] == 1
        graph.add_edge(*e)
    conserved_columns()
    for j, (u, v) in enumerate(special[:2]):
        trans(far[2 * j], u)
        trans(far[2 * j + 1], v)
        assert (loads[u], loads[v]) == (2, 0) and u < v
        # The ordinary private reference at 167 realizes the two private holds
        # in the native encoder. The final endpoints are both filled.
        loads[u] = loads[v] = 1
        graph.remove_edge(u, v)
    trans(far[4], GUARD)
    trans(far[6], DONOR)
    trans(parked, HUB)
    assert not data_edges.intersection(frozenset(e) for e in graph.edges)
    assert frozenset((HUB, GUARD)) not in data_edges
    assert tuple(loads[[HUB, TARGET, GUARD, DONOR]]) == (
        1 + sign,
        1,
        1 - sign,
        1 - sign,
    )
    assert loads.sum() == b.n
    mounted = loads.copy()
    storage = extract("tip")[1]
    # Restore the two guide rails for readout on both actual donor branches.
    for converted in (False, True):
        loads[:] = mounted
        if converted:
            loads[TARGET], loads[DONOR] = loads[DONOR], loads[TARGET]
        # Independent sharp load reads via a fresh degree-five private reference.
        for port in (HUB, TARGET, DONOR):
            assert graph.has_node(lab["reference"])
            trans(port, lab["reference"])
            trans(port, lab["reference"])
        trans(GUARD, far[8])
        trans(DONOR, far[9])
        assert loads[GUARD] == loads[DONOR] == 1
        assert loads[HUB] == 1 + sign and loads[TARGET] == (
            1 - sign if converted else 1
        )
        assert loads.sum() == b.n
        assert all(loads[p] == (1 + sign if p == HUB else 1) for p in guide)
        assert TARGET not in storage
    assert ladder.read_loads() == [0, 1, 2]
    return dict(
        hub=HUB,
        target=TARGET,
        guard=GUARD,
        donor=DONOR,
        temporary_carrier_parking=parked,
        private_target=TARGET,
        load_read_reference=lab["reference"],
        service_opening=service,
        stock_ports=far[:8],
        stock_initial_loads=labels,
        readout_storage=storage,
        material_preparation_and_read_transpositions=operations,
    )
