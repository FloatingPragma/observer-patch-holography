"""An orientation-sensitive read of an imposed one-third-turn field phase."""

import argparse
from collections import defaultdict
from hashlib import sha256
import json
import math
from pathlib import Path

import networkx as nx

from . import measure


def directed_edges(face):
    return {tuple(sorted((a, b))): 1 if a < b else -1 for a, b in measure.search.around(face)}


def oriented_cell(patch):
    faces = {0: patch.faces[0]}
    for a, b in nx.bfs_edges(patch.dual, 0):
        face = patch.faces[b]
        previous = directed_edges(faces[a])
        candidate = directed_edges(face)
        edge, = previous.keys() & candidate.keys()
        faces[b] = face[::-1] if previous[edge] == candidate[edge] else face
    for a, b in patch.dual.edges:
        left, right = directed_edges(faces[a]), directed_edges(faces[b])
        edge, = left.keys() & right.keys()
        assert left[edge] == -right[edge]
    return faces


def path_phases(patch, words, oriented, near, far):
    path = nx.shortest_path(patch.dual, near, far)
    crossings = []
    for a, b in zip(path, path[1:]):
        seam, = directed_edges(oriented[a]).keys() & directed_edges(oriented[b]).keys()
        crossings.append((patch.bit[seam], -directed_edges(oriented[a])[seam]))
    exponent = {w: sum(s for bit, s in crossings if w & bit) for w in words}
    roots = (1, complex(-0.5, math.sqrt(3) / 2), complex(-0.5, -math.sqrt(3) / 2))
    phases = {w: roots[k % 3] for w, k in exponent.items()}
    energy = [0.] * len(patch.faces)
    for w in words:
        for face, other in patch.ring_moves(w):
            mask, clockwise = patch.orientation(oriented[face])
            if w & mask == clockwise:
                assert exponent[other] - exponent[w] == int(face == near) - int(face == far)
                energy[face] += 0.5 * abs(phases[other] - phases[w]) ** 2 / len(words)
    return phases, energy


def quadratures(patch, words, final, cycle, phase):
    """Outcome probabilities after recombination, tracing out unobserved word bits."""
    start = cycle[0]
    mask, forward = patch.orientation(cycle)
    idle = {(start, w): phase[w] for w in words}
    active = {end: phase[w] for w, end in zip(words, final, strict=True)}
    moments = defaultdict(lambda: [0., 0j])
    for state in sorted(idle.keys() | active.keys()):
        pos, word = state
        orientation = {forward: 'forward', mask ^ forward: 'reverse'}.get(word & mask, 'mixed')
        a, b = idle.get(state, 0), active.get(state, 0)
        row = moments[pos, orientation]
        row[0] += abs(a) ** 2 + abs(b) ** 2
        row[1] += a.conjugate() * b
    panels = []
    for (pos, direction), (weight, overlap) in sorted(moments.items()):
        panel = {'port': pos, 'heads': direction}
        for axis, part in [('x', overlap.real), ('y', overlap.imag)]:
            panel['bright_' + axis] = round((weight + 2 * part) / (4 * len(words)), 12)
            panel['dark_' + axis] = round((weight - 2 * part) / (4 * len(words)), 12)
        panels.append(panel)
    return panels


def run():
    patch = measure.search.Patch()
    words = patch.field_component()
    # Each singleton is a physical signed cut. Summing the singleton reads
    # cancels each internal edge and gives the cut of any region.
    local_cuts = [sorted({measure.cut_flux(patch, (port,), w) for w in words})
                  for port in range(patch.cell.n)]
    orientations = oriented_cell(patch)
    settings = []
    for label, (near, far) in measure.antipodal_settings(patch):
        phases, energy = path_phases(patch, words, orientations, near, far)
        panels = []
        for name, disk in measure.contours(patch, near, far).items():
            cycle = measure.boundary(patch, disk)
            for sign in (-1, 1):
                final = measure.trace_tour(patch, words, cycle, sign)
                reads = {source: quadratures(patch, words, final, cycle, weights)
                         for source, weights in [('uniform', dict.fromkeys(words, 1)),
                                                 ('positive', phases),
                                                 ('negative', {w: z.conjugate() for w, z in phases.items()})]}
                tv_y = sum(abs(a[k] - b[k]) for a, b in zip(reads['positive'], reads['negative'], strict=True)
                           for k in ('bright_y', 'dark_y')) / 2
                panels.append({'contour': name, 'cycle': cycle, 'sign': sign, 'outcomes': reads,
                               'conjugate_total_variation_y': round(tv_y, 12)})
        settings.append({'setting': label, 'faces': [near, far], 'energy_by_face': [round(e, 12) for e in energy],
                         'panels': panels})
    here = Path(__file__).resolve()
    inputs = [here, here.with_name('PHASE_CONTROL.md'), Path(measure.__file__),
              Path(measure.search.__file__), measure.search.REPO / 'model/executable_model/carrier.py']
    return {'source_phase_turns': '1/3 and -1/3 supplied', 'field_words': len(words),
            'inputs_sha256': {str(p.relative_to(measure.search.REPO)): sha256(p.read_bytes()).hexdigest()
                              for p in inputs}, 'singleton_cut_values': local_cuts, 'settings': settings}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    output = parser.parse_args().output
    result = run()
    with output.open('x', encoding='utf-8') as file:
        json.dump(result, file, indent=2, allow_nan=False)
        file.write('\n')
