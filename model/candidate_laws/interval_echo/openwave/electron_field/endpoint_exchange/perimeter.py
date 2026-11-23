"""Exchange two identical charges across opposite sides of a physical six-port loop."""
import argparse
from itertools import product
import json
from pathlib import Path

import networkx as nx

from search import Exchange, Patch, federation, occupied


def injections(patch, gateways, starts, sign):
    graph = nx.DiGraph()
    for p in range(12):
        graph.add_edges_from((p, q) for q, _ in patch.charge_moves(p, patch.seed, sign))
    choices = [[list(path) for end in gateways for path in
                ([[v]] if end == v else nx.all_simple_paths(graph, v, end))]
               for v in starts]
    valid = [(u, v) for u, v in product(*choices) if set(u).isdisjoint(v)]
    selected = min(valid, key=lambda pair: (sum(map(len, pair)), pair))
    return [list(reversed([*path, gateways[path[-1]]])) for path in selected]


def run():
    patch = Patch()
    conjugate = Patch()
    conjugate.seed ^= sum(conjugate.bit.values())
    ring = patch.boundary
    gateways = {p: q for p, q in federation(0).external if p < 12}
    contacts = [occupied(p, ring[(i+1) % 6]) for i, p in enumerate(ring)]
    receipts = []
    for start, direction, sign in product(range(6), (-1, 1), (-1, 1)):
        sites = tuple(ring[(start + direction * i) % 6] for i in range(6))
        origins = (sites[0], sites[3])
        labels = list(origins)
        places, hops = [], []
        for step in range(6):
            who = step % 2
            source = labels[who]
            target = sites[(sites.index(source) + 1) % 6]
            assert target not in labels
            places.append(occupied(*labels))
            hops.append((source, target))
            labels[who] = target
        assert labels == list(reversed(origins))
        protocol = Exchange(patch, sites, hops, places, contacts, sign)
        record = protocol.scan()
        record['injection_paths'] = injections(patch, gateways, origins, sign)
        record['direction'] = direction
        conjugate_result = Exchange(conjugate, sites, hops, places, contacts, -sign).scan()
        assert record['words_per_position'] == conjugate_result['words_per_position']
        assert record['position_only_accepted_words'] == conjugate_result['position_only_accepted_words']
        for group, result in record['resolved'].items():
            control = conjugate_result['resolved'][group]
            for field in ('counts', 'moving_vs_parked', 'initial_face_pairs_with_both_signs',
                          'same_initial_state_examples'):
                assert result.get(field) == control.get(field)
            for outcome, branch in result['max'].items():
                # Complementing every head can rephase a template at a single position.
                # These phases cancel in the product around the closed exchange.
                assert branch['signed_numerator'] == control['max'][outcome]['signed_numerator']
        record['matched_conjugate_control'] = 'opposite charge and complementary background give identical exchange signs and probabilities'
        receipts.append(record)
    return {'physical_perimeter': ring, 'face_ports': patch.faces, 'fixtures': receipts}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = run()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
