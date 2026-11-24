"""Joint regional-head and perimeter-interference records on one written cell."""

import argparse
import cmath
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
import json
import math
from pathlib import Path

import networkx as nx

from ..electron_field.flux_string import search


def cut_flux(patch, ports, word):
    """Head-count change relative to the sharp background, using only cut edges."""
    total = 0
    for (a, b), bit in patch.bit.items():
        orientation = int(a in ports) - int(b in ports)
        if orientation:
            total += orientation * (int(bool(word & bit)) - int(bool(patch.seed & bit)))
    return total


def voltage_ports(charges):
    """Two reference phases at each pulse angle, from a material-charge distribution."""
    reads = {}
    for label, theta in [('pi_over_two', math.pi / 2), ('pi', math.pi)]:
        coherence = sum(p * cmath.exp(1j * theta * q) for q, p in charges.items())
        reads[label] = {
            'bright_x': round((1 + coherence.real) / 2, 12),
            'bright_y': round((1 + coherence.imag) / 2, 12),
        }
    return reads


def material_change(pos, region, sign, reference_port):
    return sign * (int(pos in region) - int(reference_port in region))


def boundary(patch, disk):
    counts = Counter(tuple(sorted(edge)) for f in disk for edge in search.around(patch.faces[f]))
    graph = nx.Graph(edge for edge, n in counts.items() if n == 1)
    if not nx.is_connected(graph) or any(degree != 2 for _, degree in graph.degree()):
        raise ValueError('Selected region has no simple perimeter')
    cycle, = nx.cycle_basis(graph)
    return min(tuple(order[i:] + order[:i]) for order in (cycle, cycle[::-1])
               for i in range(len(cycle)))


def contours(patch, near, far):
    disk = {f for f, vertices in enumerate(patch.faces)
            if set(vertices).intersection(patch.faces[near])}
    extra = min(f for f in patch.dual if f not in disk and f != far
                and len(disk.intersection(patch.dual[f])) == 1)
    exterior = min(set(patch.dual) - {near, far})
    return {'enclosing': disk, 'deformed': disk | {extra}, 'exterior': {exterior}}


def antipodal_settings(patch):
    diameter = nx.diameter(patch.dual)
    pairs = [(a, b) for a, b in combinations(range(len(patch.faces)), 2)
             if nx.shortest_path_length(patch.dual, a, b) == diameter]
    return zip(('primary', 'held'), pairs[:2], strict=True)


def trace_tour(patch, words, cycle, sign):
    edges = [tuple(sorted(edge)) for edge in search.around(cycle)]
    sequence = edges + edges[-2::-1]
    final = []
    for word in words:
        pos, head = cycle[0], word
        for edge in sequence:
            pos, head = patch.move(pos, head, edge, sign)
        final.append((pos, head))
    if len(set(final)) != len(words):
        raise ValueError('The controlled tour must preserve orthogonal inputs')
    return final


def joint_ports(words, final, start, cuts, amplitudes):
    """Integer bright/dark weights, with one common denominator 4 times the input size."""
    active, idle = {}, {}
    for word, (pos, head) in zip(words, final, strict=True):
        tag = cuts[word]
        idle[tag, start, word] = amplitudes[word]
        active[tag, pos, head] = amplitudes[word]
    bins = defaultdict(lambda: [0, 0])
    for key in idle.keys() | active.keys():
        tag, pos, _ = key
        a, b = idle.get(key, 0), active.get(key, 0)
        weights = bins[tag, pos]
        weights[0] += (a + b) ** 2
        weights[1] += (a - b) ** 2
    denominator = 4 * len(words)
    if sum(sum(weights) for weights in bins.values()) != denominator:
        raise ValueError('Incomplete terminal read')
    return {
        'denominator': denominator,
        'by_cuts_and_final_port': [
            {'cut': tag, 'port': pos, 'bright': bright, 'dark': dark}
            for (tag, pos), (bright, dark) in sorted(bins.items())
        ],
        'inclusive_bright': sum(weights[0] for weights in bins.values()) / denominator,
    }


def endpoint_pair(patch, words, near, far):
    regions = [patch.faces[near], patch.faces[far]]
    cuts = {w: tuple(cut_flux(patch, r, w) for r in regions) for w in words}
    mask = patch.path(near, far)
    bare = dict.fromkeys(words, 1)
    string = {w: search.phase(w, mask) for w in words}
    energy = [0] * len(patch.faces)
    for w in words:
        for face, other in patch.ring_moves(w):
            if w < other:
                energy[face] += (string[w] - string[other]) ** 2
    panels = []
    for name, disk in contours(patch, near, far).items():
        cycle = boundary(patch, disk)
        for sign in (-1, 1):
            final = trace_tour(patch, words, cycle, sign)
            voltage = {}
            for source, weights in [('bare', bare), ('string', string)]:
                # Evaluate the addressed voltage on each prepared material/head
                # basis amplitude before flight, keeping the preloaded probe.
                prepared = {(cycle[0], w): amp for w, amp in weights.items()}
                voltage[source] = []
                for region in regions:
                    material = Counter()
                    for (pos, _), amp in prepared.items():
                        material[material_change(pos, region, sign, cycle[0])] += amp ** 2
                    voltage[source].append(voltage_ports({q: n / len(words) for q, n in material.items()}))
            counts = Counter()
            winding_overlap = 0
            for w, (pos, head) in zip(words, final, strict=True):
                if pos != cycle[0]:
                    counts['away'] += 1
                elif head == w:
                    counts['idle'] += 1
                else:
                    counts['winding'] += 1
                    winding_overlap += string[w] * string[head]
            # Complete dephasing before the arm split makes the return overlap
            # the sum of diagonal entries of the addressed permutation.
            dephased = (1 + counts['idle'] / len(words)) / 2
            panels.append({
                'contour': name, 'disk': sorted(disk), 'cycle': cycle, 'sign': sign,
                'compensator_path': patch.compensator(cycle[0], sign),
                'pre_tour_material_voltage_by_region': voltage,
                'path_counts': {k: counts[k] for k in ('idle', 'winding', 'away')},
                'conditional_winding_overlap': (
                    winding_overlap / counts['winding'] if counts['winding'] else None),
                'bare': joint_ports(words, final, cycle[0], cuts, bare),
                'string': joint_ports(words, final, cycle[0], cuts, string),
                'fully_dephased_bright_both_inputs': dephased,
            })
    return {
        'faces': [near, far], 'regions': regions,
        'cut_pair_multiplicities': [
            {'cut': tag, 'words': n} for tag, n in sorted(Counter(cuts.values()).items())],
        'ring_energy_by_face': [n / (2 * len(words)) for n in energy],
        'panels': panels,
    }


def first_fork(patch, words, region, sign, inside_count):
    for word in words:
        for anchor in range(patch.cell.n):
            if anchor in region:
                continue
            tips = sorted(pos for pos, _ in patch.charge_moves(anchor, word, sign))
            for a, b in combinations(tips, 2):
                if int(a in region) + int(b in region) == inside_count:
                    return word, anchor, a, b
    raise ValueError('No eligible fork for calibration')


def partial_hop(patch, state, edge, sign, probability, write_head):
    """Exponentiate the Laplacian on each addressed two-state block."""
    angle = math.asin(math.sqrt(probability))
    stay = cmath.exp(-1j * angle) * math.cos(angle)
    swap = 1 - stay
    result = defaultdict(complex)
    for (pos, word), amplitude in state.items():
        if write_head:
            target = patch.move(pos, word, edge, sign)
        else:
            a, b = edge
            target = ({a: b, b: a}.get(pos, pos), word)
        result[pos, word] += stay * amplitude
        result[target] += swap * amplitude
    return result


def control_read(patch, state, region, sign, anchor):
    ports, flux, joint = Counter(), Counter(), Counter()
    for (pos, word), amplitude in state.items():
        probability = abs(amplitude) ** 2
        if probability < 1e-25:  # cos(pi/2) roundoff at a complete swap
            continue
        material = material_change(pos, region, sign, anchor)
        heads = cut_flux(patch, region, word)
        ports[material] += probability
        flux[heads] += probability
        joint[material, heads] += probability
    if not math.isclose(sum(joint.values()), 1, abs_tol=1e-12):
        raise ValueError('Lost probability in control')
    mean = sum(q * p for q, p in ports.items())
    return {
        'joint': [{'material': q, 'cut': h, 'probability': round(p, 12)}
                  for (q, h), p in sorted(joint.items())],
        'mean_cut': round(sum(q * p for q, p in flux.items()), 12),
        'mean_material': round(mean, 12),
        'material_voltage': voltage_ports(ports),
        'hypothetical_sharp_charge_at_same_mean': voltage_ports({mean: 1}),
    }


def calibration(patch, words, region, sign, inside_count):
    word, anchor, a, b = first_fork(patch, words, region, sign, inside_count)
    sharp = {(anchor, word): 1}
    inward = a if a in region else b
    reads = {}
    for name, write_head in [('written', True), ('unwritten', False)]:
        unit = partial_hop(patch, sharp, tuple(sorted((anchor, inward))), sign, 1, write_head)
        state = partial_hop(patch, sharp, tuple(sorted((anchor, a))), sign, 1 / 3, write_head)
        state = partial_hop(patch, state, tuple(sorted((anchor, b))), sign, 1 / 2, write_head)
        reads[name] = {'sharp': control_read(patch, unit, region, sign, anchor),
                       'spread': control_read(patch, state, region, sign, anchor)}
    return {'sign': sign, 'tips_inside': inside_count, 'word': word,
            'ports': [anchor, a, b], 'compensator_path': patch.compensator(anchor, sign),
            'reads': reads}


def run():
    patch = search.Patch()
    words = patch.field_component()
    settings = []
    for label, (near, far) in antipodal_settings(patch):
        result = endpoint_pair(patch, words, near, far)
        result['setting'] = label
        result['calibrations_by_region'] = {
            end: [calibration(patch, words, patch.faces[face], sign, count)
                  for sign in (-1, 1) for count in (1, 2)]
            for end, face in [('near', near), ('far', far)]}
        settings.append(result)
    here = Path(__file__).resolve()
    inputs = [here, here.with_name('DESIGN.md'), Path(search.__file__),
              search.REPO / 'model/executable_model/carrier.py',
              search.REPO / 'model/executable_model/carriers/FEDERATION_LEVEL0.npz']
    return {'law': 'positive written-field Laplacian, fixed one-cell component',
            'inputs_sha256': {str(p.relative_to(search.REPO)): sha256(p.read_bytes()).hexdigest()
                              for p in inputs},
            'field_words': len(words), 'seed_head_word': patch.seed, 'settings': settings}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    target = parser.parse_args().output
    result = run()
    with target.open('x', encoding='utf-8') as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write('\n')
