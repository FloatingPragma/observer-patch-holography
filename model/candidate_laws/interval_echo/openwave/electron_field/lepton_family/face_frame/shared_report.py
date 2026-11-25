"""Send one recorded energy/weak event from Bob to Alice over physical Interval Echo seams."""

import argparse
import json
from itertools import pairwise

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_gauge.flag_measure.measure import meter
from . import observers

HERE = observers.HERE
SHOTS = observers.SHOTS
SEED = 178406
TAIL = 0.01 / (2 * 4 * (6 + 3))
REFUSED = 6
INITIAL = (0, 1, 2, 0, 2)
WORDS = [
    (INITIAL[j:3] + INITIAL[:j]) + tuple(INITIAL[3:][:: 1 if w == 0 else -1])
    for j in range(3)
    for w in range(2)
]


def certify():
    law = candidate.native_law
    ex, con = law.local_term_parts()
    full = expm(-1j * np.pi * (ex + candidate.KAPPA * con))
    direct = np.eye(99)[:, [law.SRC.index((b, a)) for a, b in law.SRC]]
    np.testing.assert_allclose(full[:, :9], direct, atol=2e-13)
    assert meter() == [0, 1, 2]
    # Check the two ternary reference reads in both physical endpoint orders.
    echo = law.echo()
    bright = np.tile([w == (1, 1) for w in law.SRC], 11)
    lookup = {}
    for value in range(3):
        for ref in [0, 2]:
            for reversed_edge in [False, True]:
                pair = (ref, value) if reversed_edge else (value, ref)
                word = np.eye(99)[:, law.SRC.index(pair)]
                out = echo @ word
                click = value + ref == 2
                np.testing.assert_allclose(
                    echo.T @ (out * (bright == click)), word, atol=1e-13
                )
                lookup[value, ref, reversed_edge] = click
    return lookup


def path(graph, a, b):
    return min(nx.all_shortest_paths(graph, a, b))


def layout():
    board = federation(0)
    graph = nx.Graph(board.seams)
    owners = nx.Graph((a // 12, b // 12) for a, b in board.external)
    dist = nx.single_source_shortest_path_length(owners, 0)
    alice = 0
    bob = min(b for b, t in dist.items() if t == max(dist.values()))
    readers = {}
    for label, owner in [("alice", alice), ("bob", bob)]:
        readers[label] = sorted(
            (int(a), int(b))
            for edge in board.external
            for a, b in [edge, edge[::-1]]
            if a // 12 == owner
        )[:2]
    protected = {r for pairs in readers.values() for _, r in pairs}
    assert len(protected) == 4
    bob_ports = list(range(12 * bob, 12 * (bob + 1)))
    alice_ports = list(range(12 * alice, 12 * (alice + 1)))
    keep = [p for p in bob_ports if board.degree[p] == 5][:5]
    send = [p for p in bob_ports if p not in keep][:5]
    receive = [p for p in alice_ports if board.degree[p] == 5][:5]
    route_graph = graph.copy()
    route_graph.remove_nodes_from(protected | set(keep))
    transfers = [path(route_graph, a, b) for a, b in zip(send, receive, strict=True)]
    reads = {}
    for label, bank, region in [
        ("bob", keep, bob_ports),
        ("alice", receive, alice_ports),
    ]:
        sub = graph.subgraph(region)
        reads[label] = [[(path(sub, p, g), r) for g, r in readers[label]] for p in bank]
    return graph, dict(
        alice=alice,
        bob=bob,
        retained=keep,
        outbox=send,
        inbox=receive,
        readers=readers,
        transfer_paths=transfers,
        read_paths=reads,
    )


def transpose(word, trail):
    # Move the two endpoints through native full on-bond swaps; every transit
    # load returns to its own port, including registers carrying other digits.
    edges = list(pairwise(trail))
    for u, v in edges + edges[:-1][::-1]:
        word[u], word[v] = word[v], word[u]


def encode(word, targets, initial_source, tag):
    # The classical field pointer controls this local permutation. All six
    # words have the same load histogram, total and input stocks.
    energy, weak = divmod(tag, 2)
    for _ in range(energy):
        for a, b in [(0, 1), (1, 2)]:
            transpose(word, path(initial_source, targets[a], targets[b]))
    if weak:
        transpose(word, path(initial_source, targets[3], targets[4]))
    assert tuple(word[targets]) == WORDS[tag]


def read(word, probes, lookup):
    before = word.copy()
    symbols = []
    for pair in probes:
        clicks = []
        for trail, r in pair:
            transpose(word, trail)
            gate = trail[-1]
            key = (int(word[gate]), int(word[r]), gate > r)
            clicks.append(lookup[key])
            transpose(word, trail)
        if clicks not in ([True, False], [False, False], [False, True]):
            return REFUSED
        symbols.append(2 if clicks[0] else 0 if clicks[1] else 1)
    np.testing.assert_array_equal(word, before)
    record = tuple(symbols)
    return WORDS.index(record) if record in WORDS else REFUSED


def channel():
    lookup = certify()
    graph, geo = layout()
    bob = graph.subgraph(range(12 * geo["bob"], 12 * (geo["bob"] + 1)))
    stocks = np.ones(len(graph), dtype=int)
    for pair in geo["readers"].values():
        stocks[pair[0][1]], stocks[pair[1][1]] = 0, 2
    for bank in (geo["retained"], geo["outbox"]):
        stocks[bank] = INITIAL
    assert stocks.sum() == len(stocks)
    actions = []
    for tag in range(6):
        word = stocks.copy()
        encode(word, geo["retained"], bob, tag)
        # The transmitter reads the retained material word; it is given no
        # access to the field preparation, its probabilities, or a second tag.
        recorded = read(word, geo["read_paths"]["bob"], lookup)
        assert recorded == tag
        encode(word, geo["outbox"], bob, recorded)
        silent = read(word, geo["read_paths"]["alice"], lookup)
        assert silent == REFUSED
        sent = word.copy()
        for trail in geo["transfer_paths"]:
            assert all(graph.has_edge(a, b) for a, b in pairwise(trail))
            transpose(sent, trail)
        alice = read(sent, geo["read_paths"]["alice"], lookup)
        retained = read(sent, geo["read_paths"]["bob"], lookup)
        expected = word.copy()
        expected[geo["outbox"]] = 1
        expected[geo["inbox"]] = WORDS[tag]
        np.testing.assert_array_equal(sent, expected)
        assert sent.sum() == len(sent)
        # Reverse the orientation convention on the three face characters.
        # A copied classical label can be misdecoded without altering a field.
        flipped = 2 * ((-(alice // 2)) % 3) + alice % 2 if alice != REFUSED else REFUSED
        actions.append(
            dict(
                field_pointer=tag,
                bob=retained,
                alice=alice,
                opposite_face=flipped,
                no_send=silent,
            )
        )
    return geo, actions


def run():
    geo, actions = channel()
    (_, ga), (_, gb) = observers.frames()
    field = observers.joint_read(ga, gb, np.random.default_rng(0))
    rng = np.random.default_rng(SEED)
    rows = []
    for setting in field:
        prob = np.array(setting["agreed_probabilities"])
        # A single six-bin draw of the effective pointer for each bank.
        # The material report supplies the conditional *joint* outcomes.
        latent = rng.multinomial(SHOTS, prob)
        bank = {}
        for control in ("alice", "opposite_face", "no_send"):
            counts = np.zeros((7, 7), dtype=np.int64)
            theory = np.zeros((7, 7))
            for n, p, r in zip(latent, prob, actions, strict=True):
                a, b = r["bob"], r[control]
                counts[a, b] += n
                theory[a, b] += p
            valid = counts[:6, :6]
            disagree = int(valid.sum() - np.trace(valid))
            refusal = int(counts[6, :].sum() + counts[:6, 6].sum())
            bank[control] = dict(
                joint=counts.tolist(),
                predicted_joint=theory.tolist(),
                agreements=int(np.trace(valid)),
                disagreements=disagree,
                transmission_refusals=refusal,
                disagreement99=observers.binomial(disagree, SHOTS, TAIL),
            )
        rows.append(
            dict(
                charge=setting["charge"],
                source_port=setting["source"],
                field_event_counts=latent.tolist(),
                field_event_probabilities=prob.tolist(),
                field_event99=[observers.binomial(int(c), SHOTS, TAIL) for c in latent],
                settings=bank,
            )
        )
    return dict(
        starts_per_bank=SHOTS,
        seed=SEED,
        material=geo,
        address_words=WORDS,
        conditional_channel=actions,
        supplied_input="one orthogonal six-outcome pointer of the declared effective field",
        supplied_services="field-pointer-controlled local write; scheduled blank contacts, guarded stocks, and two oriented pointer references per owner",
        digit_order="energy character 0,1,2, each with weak click then weak refusal; final index 6 = invalid report",
        controls_share_input_events=True,
        rows=rows,
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    r = run()
    if ap.parse_args().write:
        (HERE / "shared_report.json").write_text(json.dumps(r, indent=2) + "\n")
    else:
        for row in r["rows"]:
            print(
                row["source_port"],
                row["charge"],
                [
                    (k, v["disagreements"], v["transmission_refusals"])
                    for k, v in row["settings"].items()
                ],
            )
