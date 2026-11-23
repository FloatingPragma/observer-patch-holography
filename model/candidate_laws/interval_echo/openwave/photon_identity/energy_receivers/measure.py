"""Read local field energy and coincidences after a signed source on an intact sphere."""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[5]))
from model.candidate_laws.interval_echo.openwave.photon_identity import audit
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.decay.measure import ring_projection

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

SHOTS = 2**20
TIMES = (0, 2, 8, 32, 128)
TOL = 2e-11


def energy_readers(original, sphere):
    """Joint eigen-projectors on rings with disjoint physical owners."""
    rings = []
    owners = [set(port // 12 for port in ring) for ring in sphere.faces]
    for face in range(sphere.n):
        pairs = (original.face_tag == face) & (original.src < original.dst)
        rings.append((original.src[pairs], original.dst[pairs]))
    squares = {}
    lookup = {word: i for i, word in enumerate(original.words)}
    for a, b in combinations(range(sphere.n), 2):
        if owners[a] & owners[b]:
            continue
        quads = []
        low_bits = (1 << a) | (1 << b)
        for i, j in zip(*rings[a], strict=True):
            w = original.words[i]
            v = original.words[j]
            if w & (1 << a):
                w, v, i, j = v, w, j, i
            if w & (1 << b):
                continue
            wb = w | (1 << b)
            wab = w | low_bits
            if wb in lookup and wab in lookup:
                quads.append((i, j, lookup[wb], lookup[wab]))
        squares[a, b] = np.array(quads, dtype=int).reshape(-1, 4)
    return owners, rings, squares


def square_probability(wave, quads):
    if len(quads) == 0:
        return 0.
    amplitudes = wave[quads]
    differences = amplitudes[:, 0] - amplitudes[:, 1] - amplitudes[:, 2] + amplitudes[:, 3]
    return float(np.vdot(differences, differences).real / 4)


def single_probabilities(wave, rings):
    return np.array([float(np.sum(abs(wave[i] - wave[j])**2) / 2) for i, j in rings])


def joint_record(p_stay, click_a, click_b, both):
    stay = [p_stay - click_a - click_b + both, click_b - both, click_a - both, both]
    table = np.array(stay + [1 - p_stay, 0., 0., 0.])
    assert np.min(table) >= -TOL and abs(sum(table) - 1) < TOL
    table = np.maximum(table, 0.)
    return table / table.sum()


def link_interface(instrument):
    """Can one material hop transport both independent head inputs coherently?"""
    from measure import repaired_instrument

    rows = []
    positions = np.eye(9)
    for label, h in [('contact', instrument.h8 / 8), ('positive_extension', repaired_instrument(instrument).h8 / 8)]:
        pulse = expm(-1j * np.pi / 8 * h)
        for sign in (1, -1):
            for reverse in (False, True):
                start = (1, 1+sign) if reverse else (1+sign, 1)
                incoming = np.kron(instrument.basis, positions[:, [instrument.pairs.index(start)]])
                outgoing = np.kron(instrument.basis, positions[:, [instrument.pairs.index(start[::-1])]])
                blocked = int(not ((sign == 1) != reverse))
                idle = incoming[:, blocked]
                eigenvalue = np.vdot(idle, h @ idle)
                # One incoming written orientation stays at the original matter endpoint for all pulse lengths.
                np.testing.assert_allclose(h @ idle, eigenvalue * idle, atol=TOL, rtol=0)
                crossing = outgoing.conj().T @ pulse @ incoming
                singular = np.linalg.svd(crossing, compute_uv=False)
                np.testing.assert_allclose(singular, [1., 0.], atol=TOL, rtol=0)
                rows.append(dict(law=label, charge=sign, reverse=reverse, blocked_head=blocked,
                                 blocked_energy=float(eigenvalue.real), crossing_singular_values=singular.tolist(),
                                 operator_norm_distance_to_any_lossless_two_head_link_at_least=1))
    return rows


def board(level, instrument, rng):
    from flight import LaplacianSector, Sector, Sphere
    from waves import federation

    sphere = Sphere(federation(level), "depth-first within shell")
    original = Sector(sphere)
    transitions = original.verify_native_moves(sphere)
    sector = LaplacianSector(original)
    sources, _, _, distances = audit.station_matrix(sphere, instrument)
    audit.verify_source_cuts(sphere, original, sources)
    owners, rings, squares = energy_readers(original, sphere)
    ground = np.ones(len(original.words)) / np.sqrt(len(original.words))
    # On each commuting four-word square the joint projector has rank one.
    # Verify this basis-level geometry against successive ring operations.
    verify = rng.normal(size=len(ground)) + 1j * rng.normal(size=len(ground))
    verify /= np.linalg.norm(verify)
    first = [ring_projection(verify, ring) for ring in rings]
    for (a, b), quads in squares.items():
        forward = ring_projection(first[a], rings[b])
        backward = ring_projection(first[b], rings[a])
        np.testing.assert_allclose(forward, backward, atol=TOL, rtol=0)
        np.testing.assert_allclose(np.vdot(forward, forward).real, square_probability(verify, quads),
                                   atol=TOL, rtol=0)
    assert max(single_probabilities(ground, rings)) == 0
    np.testing.assert_allclose(sector.h @ ground, 0, atol=TOL)
    starts = []
    eligible = {}
    for face, seam in sources.items():
        x = sector.flipped[:, face] ^ bool(sphere.head[seam])
        # For opposite signs the two unnormalized stay branches partition the ground.
        # Moved branches are these selectors on the cut graph with the source seam reversed.
        starts.extend([ground * (~x), ground * x])
        for launch in starts[-2:]:
            np.testing.assert_allclose(sector.hamiltonian(cut=(face,)) @ launch, 0., atol=TOL)
            clicks = single_probabilities(launch, rings)
            np.testing.assert_allclose(np.delete(clicks, face), 0., atol=TOL)
        eligible[face] = [(a, b) for a, b in squares if not owners[face] & (owners[a] | owners[b])]
    initial = np.column_stack(starts)
    energies = np.sum(initial * (sector.h @ initial), axis=0)
    best = None
    records = []
    for t in TIMES:
        evolved = expm_multiply((-1j * t) * sector.h, initial) if t else initial
        for s, face in enumerate(sources):
            plus, minus = evolved[:, 2*s], evolved[:, 2*s+1]
            np.testing.assert_allclose(plus + minus, ground, atol=TOL, rtol=0)
            clicks = single_probabilities(plus, rings)
            np.testing.assert_allclose(clicks, single_probabilities(minus, rings), atol=TOL, rtol=0)
            np.testing.assert_allclose(clicks.sum(), energies[2*s:2*s+2], atol=TOL, rtol=0)
            pairs = eligible[face]
            pair_probs = np.array([square_probability(plus, squares[pair]) for pair in pairs])
            winner = int(np.argmax(pair_probs))
            a, b = pairs[winner]
            joint = pair_probs[winner]
            remote = [f for f in range(sphere.n) if not owners[face] & owners[f]]
            receiver = max(remote, key=lambda f: clicks[f])
            records.append(dict(source=face, time=t, energy=float(clicks.sum()),
                                all_face_click_probabilities=clicks.tolist(),
                                strongest_disjoint_receiver=receiver,
                                remote_energy=float(clicks[remote].sum()),
                                tested_pairs=len(pairs), pairs_with_two_click_support=sum(bool(len(squares[p])) for p in pairs),
                                pair_probabilities=pair_probs.tolist(),
                                strongest_pair=[a, b], strongest_double=float(joint)))
            if best is None or joint > best[0]:
                heralds = np.linalg.norm(initial[:, 2*s:2*s+2], axis=0)**2
                best = joint, a, b, face, t, clicks[a], clicks[b], heralds
    joint, a, b, face, t, ca, cb, heralds = best
    calibration = np.zeros_like(ground)
    assert len(squares[a, b]) > 0
    calibration[squares[a, b][0]] = np.array([1, -1, -1, 1]) / 2
    np.testing.assert_allclose(single_probabilities(calibration, rings)[[a, b]], [1, 1], atol=TOL)
    np.testing.assert_allclose(square_probability(calibration, squares[a, b]), 1, atol=TOL)
    counts = {}
    for name, table in [('positive', joint_record(heralds[0], ca, cb, joint)),
                        ('negative', joint_record(heralds[1], ca, cb, joint)),
                        ('sham', joint_record(1., 0., 0., 0.)),
                        ('source_face_disabled_positive', joint_record(heralds[0], 0., 0., 0.)),
                        ('calibrated_two_rings', joint_record(1., 1., 1., 1.))]:
        counts[name] = dict(probabilities=table.tolist(), counts=rng.multinomial(SHOTS, table).tolist())
    return dict(cells=sphere.board.cells, reachable_words=len(ground), field_faces=sphere.n,
                checked_directed_ring_transitions=transitions, sources=list(sources),
                eligible_receiver_pairs={str(s): ps for s, ps in eligible.items()},
                outcome_order=['stay00', 'stay01', 'stay10', 'stay11', 'move00', 'move01', 'move10', 'move11'],
                survey=records, held_counts=dict(source=face, receiver_faces=[a, b], time=t,
                    face_distances=[distances[face][a], distances[face][b]], source_owner_cells=sorted(owners[face]),
                    receiver_owner_cells=[sorted(owners[a]), sorted(owners[b])],
                    square_rank=len(squares[a, b]), unconditional_singles=[float(ca), float(cb)],
                    unconditional_double=float(joint), ideal_calibration_word_indices=squares[a, b][0].tolist(),
                    samples=counts))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    with TemporaryDirectory() as work:
        manifest = audit.checked_archive(Path(work))
        from waves import Instrument, pulse_oracle
        instrument = Instrument()
        pulse_oracle(instrument)
        rng = np.random.default_rng(837142)
        receipt = dict(archive_files=len(manifest), shots_per_arm=SHOTS, seed=837142,
                       design_sha256=sha256((HERE/'DESIGN.md').read_bytes()).hexdigest(),
                       code_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                       single_pulse_link_comparison=link_interface(instrument),
                       boards=[board(level, instrument, rng) for level in (1, 2)])
    with args.output.open('x') as output:
        json.dump(receipt, output, separators=(',', ':'), allow_nan=False)
        output.write('\n')


if __name__ == '__main__':
    main()
