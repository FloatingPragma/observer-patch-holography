"""Signed contact and collective colour reads on a released three-load face."""

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import block_diag, expm

from model.candidate_laws.contrast_memory_spherical_v1 import law, sectors
from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier
from model.experiments.interval_echo_gauge import __main__ as gauge
from model.experiments.interval_echo_gauge.singlet_flight.measure import colour_meter


def contact_calibration(kappa):
    echo = law.echo()
    exchange, cost = law.local_term_parts()
    h = exchange + kappa * cost
    pulse = expm(-1j * np.pi / 32 * h)
    swap = expm(-1j * np.pi * h)[:, :9]
    permutation, _ = law._load_ops()
    np.testing.assert_allclose(swap[:9], permutation, atol=8e-13)
    np.testing.assert_allclose(swap[9:], 0, atol=8e-13)
    table = {}
    for column, (n, reference) in enumerate(law.SRC):
        words = abs(echo[9:, column].reshape(len(law.LABELS), 9)) ** 2
        signs = Counter()
        for (_, difference), p in zip(law.LABELS, words.sum(axis=1), strict=True):
            signs[str(difference)] += p
        # The reference-arm overlap includes the full unitary on the active arm.
        overlap = pulse[column, column]
        table[f'{n},{reference}'] = dict(
            signed_echo=dict(signs), bright_x=(1 + overlap.real) / 2,
            bright_y=(1 + overlap.imag) / 2, changed_word=1 - abs(overlap) ** 2)
    return table


def signed_ports(weights, sector, triangle, calibration, reference_ports):
    reads = []
    for p in triangle:
        occupied = sum(weights[i] for i, sites in enumerate(sector.states) if p in sites)
        controls = {}
        for reference in (0, 1, 2):
            empty, loaded = (calibration[f'{n},{reference}'] for n in (0, 1))
            controls[str(reference)] = {
                k: (1 - occupied) * empty[k] + occupied * loaded[k]
                for k in ('bright_x', 'bright_y', 'changed_word')}
            controls[str(reference)]['signed_echo'] = {
                key: (1 - occupied) * empty['signed_echo'][key] + occupied * loaded['signed_echo'][key]
                for key in empty['signed_echo']}
        reads.append(dict(port=p, reference_port=reference_ports[p], occupied=occupied, probes=controls))
    return reads


def routing(board, triangle):
    owners = sorted({b // 12 for a, b in board.external if a < 12})
    assert len(owners) == 3
    slot = {p: i for i, p in enumerate(triangle)}
    remaining = sorted(set(range(12)) - set(slot))
    slot.update({p: i % 3 for i, p in enumerate(remaining)})
    destination = [12 * owners[slot[p]] + p for p in range(12)]
    network = nx.Graph(board.seams)
    labels = np.arange(board.n)
    expected = labels.copy()
    paths = []
    reference = []
    for p, dest in enumerate(destination):
        path = nx.shortest_path(network, p, dest)
        edges = list(zip(path, path[1:]))
        for a, b in edges + edges[-2::-1]:
            labels[a], labels[b] = labels[b], labels[a]
        expected[p], expected[dest] = expected[dest], expected[p]
        paths.append(path)
        reference.append(min(q for q in board.neighbours[dest]
                             if q // 12 == dest // 12 and q not in destination))
    np.testing.assert_array_equal(labels, expected)
    return slot, reference, dict(receivers=owners, destinations=destination, transposition_paths=paths)


def colour_words():
    P, centre, antipodes, info = colour_meter()
    _, sources, _, _, bands, _, _, _, addresses, _ = currents._service(carrier.cell(), False, 1.)
    embedding = np.column_stack(bands)
    axis = bands[0].T @ sources[:, 0]
    axis /= np.linalg.norm(axis)
    generator = block_diag(np.eye(3) / 3 - np.outer(axis, axis), np.zeros((3, 3)))
    coeff = np.linalg.lstsq(gauge.coordinates(addresses).T, gauge.coordinates([generator])[0], rcond=1e-12)[0]
    half = gauge.ordered_word(addresses, coeff, np.pi, info['pulse_subdivisions'])
    undo = gauge.ordered_word(addresses, coeff, -2 * np.pi, info['pulse_subdivisions'])
    ident = np.eye(12)
    half = ident + embedding @ (half - np.eye(6)) @ embedding.T
    inverse = ident + embedding @ (undo - np.eye(6)) @ embedding.T
    ideal_half = ident + embedding @ (expm(1j * np.pi * generator) - np.eye(6)) @ embedding.T
    info['half_word_error'] = np.linalg.norm(ideal_half - half, 2)
    info['inverse_error'] = np.linalg.norm(inverse @ centre - ident, 2)
    assert max(info['half_word_error'], info['compiled_centre_error']) < 1e-4
    singlet = np.zeros((12, 12, 12))
    for cols in permutations(range(3)):
        sign = (-1) ** sum(a > b for a, b in combinations(cols, 2))
        singlet += sign * np.einsum('i,j,k->ijk', *(bands[0][:, j] for j in cols)) / np.sqrt(6)
    axis_port = bands[0] @ axis
    orbital = P[:, 1] - axis_port * (axis_port @ P[:, 1])
    orbital /= np.linalg.norm(orbital)
    control = np.einsum('i,j,k->ijk', orbital, orbital, orbital)
    words = {'sham': ident, 'centre': centre, 'centre_and_inverse': inverse @ centre, 'half_turn': half}
    return P, singlet, words, antipodes, info, control


def transform(ops, state):
    return np.einsum('ia,jb,kc,abc->ijk', *ops, state, optimize=True)


def colour_read(state, P, singlet, words):
    mass = np.vdot(state, state).real
    classes = np.zeros(4)
    for bits in product((0, 1), repeat=3):
        projected = transform([P if b else np.eye(12) - P for b in bits], state)
        classes[sum(bits)] += np.vdot(projected, projected).real
    ports = {}
    for word, u in words.items():
        out = transform([u] * 3, state)
        ports[word] = {'bright': np.linalg.norm(out + state) ** 2 / 4,
                       'dark': np.linalg.norm(out - state) ** 2 / 4}
        np.testing.assert_allclose(sum(ports[word].values()), mass, atol=2e-11)
    epsilon = abs(np.vdot(singlet, state)) ** 2
    return dict(accepted=mass, triplet_count_weights=classes.tolist(),
                three_triplet_singlet=epsilon, full_singlet=classes[0] + epsilon, ports=ports)


def sorted_slots(sites, slot):
    return tuple(sorted(sites, key=slot.__getitem__))


def route_state(psi, sector, slot):
    tensor = np.zeros((12, 12, 12), complex)
    occupancy = Counter()
    for sites, amp in zip(sector.states, psi, strict=True):
        counts = tuple(sum(slot[p] == j for p in sites) for j in range(3))
        occupancy[','.join(map(str, counts))] += abs(amp) ** 2
        if counts == (1, 1, 1):
            tensor[sorted_slots(sites, slot)] = amp
    return tensor, dict(sorted(occupancy.items()))


def face_bound(faces, slot, P, antipodes, word_error):
    addresses = [sorted_slots(face, slot) for face in faces if len({slot[p] for p in face}) == 3]
    ideal = np.eye(12) + (np.exp(2j * np.pi / 3) - 1) * P
    columns = np.asarray(addresses).T
    sub = np.prod([ideal[np.ix_(c, c)] for c in columns], axis=0)
    effect = (2 * np.eye(len(addresses)) - sub - sub.conj().T) / 4

    # The ideal dark effect is 3/4 (sum P_i - sum_{i<j} P_i P_j).
    # Every diagonal is 27/64. Bound the absolute off-diagonals by rationals.
    radii = []
    for x in addresses:
        radius = Fraction(0)
        for y in addresses:
            off = [Fraction(1, 4) if antipodes[a] == b else Fraction(1, 8)
                   for a, b in zip(x, y, strict=True) if a != b]
            if len(off) == 1:
                radius += Fraction(3, 8) * off[0]
            elif len(off) == 2:
                radius += Fraction(3, 4) * off[0] * off[1]
        radii.append(radius)
    lower = Fraction(27, 64) - max(radii)
    return dict(accepted_face_addresses=addresses, ideal_min_eigenvalue=np.linalg.eigvalsh(effect)[0],
                ideal_rational_lower_bound=str(lower),
                compiled_lower_estimate=float(lower) - 1.5 * word_error - 1e-9)


def panel(psi, source_accept, sector, faces, triangle, slot, refs, calibration, apparatus):
    P, singlet, words = apparatus
    weights = abs(psi) ** 2
    face_weights = [weights[sector.index[f]] for f in faces]
    tensor, occupancy = route_state(psi * np.sqrt(source_accept), sector, slot)
    colour = colour_read(tensor, P, singlet, words)
    np.testing.assert_allclose(sum(occupancy.values()) + 1 - source_accept, 1, atol=1e-12)
    birth = sector.index[triangle]
    selected = np.zeros_like(weights)
    selected[birth] = weights[birth]
    conditioned = (signed_ports(selected / weights[birth], sector, triangle, calibration, refs)
                   if weights[birth] else None)
    return dict(source_failure=1 - source_accept, occupancy_per_launch=occupancy,
                face_probabilities_given_source=face_weights, any_face_given_source=sum(face_weights),
                target_ports_given_source=signed_ports(weights, sector, triangle, calibration, refs),
                birth_face_herald_probability=source_accept * weights[birth],
                target_ports_given_birth_face=conditioned,
                colour_per_launch=colour)


def run():
    cell, board = carrier.cell(), carrier.federation(0)
    faces = cell.faces()
    sector = sectors.BinarySector(cell, 3)
    P, singlet, words, antipodes, info, symmetric = colour_words()
    apparatus = (P, singlet, words)
    calibrations = {k: contact_calibration(k) for k in (candidate.KAPPA, 0)}
    cases = []
    for triangle in faces[:2]:
        slot, refs, hardware = routing(board, triangle)
        lower = face_bound(faces, slot, P, antipodes, info['compiled_centre_error'])
        runs = []
        for kappa, calibration in calibrations.items():
            h = sector.hamiltonian(kappa)
            energy, vectors = np.linalg.eigh(h)
            sharp = np.eye(len(h))[:, sector.index[triangle]]
            flight = vectors @ (np.exp(-1j * energy) * (vectors.T @ sharp))
            ground = vectors[:, 0]
            assert energy[1] - energy[0] > 1e-9
            for source, psi, prob in [('sharp', sharp, 1), ('flight_at_one', flight, 1),
                                      ('ground_filter', ground, abs(ground @ sharp) ** 2)]:
                data = panel(psi, prob, sector, faces, triangle, slot, refs, calibration, apparatus)
                runs.append(dict(kappa=kappa, source=source, **data))
        cases.append(dict(birth_face=triangle, routing=hardware, all_accepted_faces=lower, runs=runs))
    return currents._portable(dict(
        protocol_sha256=sha256(Path(__file__).with_name('DESIGN.md').read_bytes()).hexdigest(),
        board_sha256=sha256((carrier.HERE / 'carriers/FEDERATION_LEVEL0.npz').read_bytes()).hexdigest(),
        face_order=faces, colour_calibration=info, signed_calibration=calibrations,
        receiver_controls={name: colour_read(state, *apparatus)
                           for name, state in [('antisymmetric', singlet), ('three_identical', symmetric)]},
        cases=cases))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = run()
    with args.output.open('x') as file:
        json.dump(result, file, indent=2, allow_nan=False)
        file.write('\n')
