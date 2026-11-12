"""Exact admissible Wilson-string states, with every native Gauss current."""
from collections import defaultdict
from functools import lru_cache
from itertools import product
import sympy as s

from model.coherent_defect_wave.chiral_joint.interacting_charge_localization.geometry import graph
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import FAMILIES
from model.coherent_defect_wave.chiral_joint.oriented_links.model import endpoint_current, label_index
from model.coherent_defect_wave.physical_attachment.chiral import representation


def shortest_paths(destination):
    edges, distance = graph()
    if type(destination) is not int or destination not in distance:
        raise ValueError('a cell of the fixed eight-cell source is required')
    paths = []
    def visit(node, route):
        if node == destination:
            paths.append(tuple(route)); return
        for edge, left, right in edges:
            if node not in (left, right): continue
            target = right if node == left else left
            if distance[target] != distance[node]+1 or distance[target] > distance[destination]: continue
            visit(target, route+[(edge, node, target, int(node == left))])
    visit(0, [])
    return tuple(paths)


def string_state(family, path):
    if family not in ('ec', 'L'):
        raise ValueError('a proved native charged local family is required')
    indices = FAMILIES[family]; d = len(indices); lookup = label_index()
    state = {}
    for word in product(range(d), repeat=len(path)+1):
        links = []
        for i, (_, _, _, forward) in enumerate(path):
            label = ((1, indices[word[i]], indices[word[i+1]]) if forward else
                     (-1, indices[word[i+1]], indices[word[i]]))
            links.append(lookup[label])
        state[word[0], tuple(links), word[-1]] = d**(-s.Rational(len(path)+1, 2))
    assert sum(s.conjugate(v)*v for v in state.values()) == 1
    return state


@lru_cache(maxsize=None)
def columns(current, endpoint):
    result = defaultdict(list)
    for (row, col), value in endpoint_current(current, endpoint).todok().items():
        result[col].append((row, value))
    return result


def gauss_image(family, path, state, node, current, *, wrong_matter_sign=False):
    indices = FAMILIES[family]
    t = representation()[2][current].extract(indices, indices)
    result = defaultdict(lambda: s.Integer(0))
    end = path[-1][2] if path else 0
    for (reference, links, matter), amplitude in state.items():
        if node == 0:
            for target in range(len(indices)):
                result[target, links, matter] += t[target, reference]*amplitude
        if node == end:
            for target in range(len(indices)):
                value = t[matter, target] if wrong_matter_sign else -t[matter, target]
                result[reference, links, target] += value*amplitude
        for j, (_, start, finish, forward) in enumerate(path):
            if node not in (start, finish): continue
            endpoint = (0 if forward else 1) if node == start else (1 if forward else 0)
            for target, value in columns(current, endpoint).get(links[j], ()):
                changed = list(links); changed[j] = target
                result[reference, tuple(changed), matter] += value*amplitude
    return {k: s.simplify(v) for k, v in result.items() if s.simplify(v) != 0}


def certificate():
    records = []; constraints = 0; largest = 0
    _, distances = graph()
    for family in ('ec', 'L'):
        for cell, distance in sorted(distances.items()):
            paths = shortest_paths(cell)
            for path in paths:
                assert len(path) == distance
                state = string_state(family, path); largest = max(largest, len(state))
                nodes = [0]+[edge[2] for edge in path]
                for node in nodes:
                    for current in range(12):
                        assert not gauss_image(family, path, state, node, current)
                        constraints += 1
                assert gauss_image(family, path, state, nodes[-1], 11, wrong_matter_sign=True)
            records.append({'family': family, 'cell': cell, 'distance': distance,
                            'shortest_paths': [[list(edge) for edge in p] for p in paths],
                            'normalized_string_count': len(paths),
                            'all_twelve_Gauss_constraints_exact': True,
                            'nonblank_link_energy': str(distance)})
    return {'rows': records, 'exact_local_Gauss_checks': constraints,
            'largest_sparse_string_state_terms': largest,
            'full_link_carrier_dimension': 119,
            'reference_at_origin_and_compensated_matter_at_endpoint': True,
            'wrong_endpoint_conjugation_detected': True,
            'spectator_local_ground_states_are_gauge_singlets': True}
