"""Measure hard downward transitions on the existing phase-wall daughter source."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix, eye
from scipy.sparse.linalg import expm_multiply

from .. import experiment as ex, walls, wall_flight

HERE = Path(__file__).resolve().parent
SHOTS, SEED, BINS, STEP = 2**26, 81148, 2**15, 1 / 64
TIMES = (8, 32)
ARMS = ((8, False), (8, True), (0, False))


def check(condition, message):
    if not condition:
        raise ArithmeticError(message)


def sample(raw, rng):
    probabilities = np.asarray(raw, dtype=float)
    check(np.all(np.isfinite(probabilities)) and min(probabilities) >= 0,
          "negative or nonfinite detector probabilities")
    error = float(abs(probabilities.sum() - 1))
    check(error < 3e-10, "detector normalization error exceeds budget")
    last = int(np.argmax(probabilities))
    positive = probabilities.copy()
    positive[last] = 1 - sum(probabilities[i] for i in range(len(positive)) if i != last)
    check(positive[last] >= 0, "negative explicit complementary probability")
    order = [i for i in range(len(positive)) if i != last] + [last]
    counts = np.empty(len(positive), dtype=np.int64)
    counts[order] = rng.multinomial(SHOTS, positive[order])
    return dict(probabilities=probabilities.tolist(), counts=counts.tolist(),
                sampling_l1_complement_bound=float(np.sum(abs(positive - probabilities))))


def compare(actual, expected, path="daughter_spectrum"):
    """Read-only receipt check that also fails with Python optimization enabled."""
    if isinstance(expected, dict):
        check(isinstance(actual, dict) and actual.keys() == expected.keys(), path)
        for key in expected:
            compare(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        check(isinstance(actual, list) and len(actual) == len(expected), path)
        for i, (a, b) in enumerate(zip(actual, expected, strict=True)):
            compare(a, b, f"{path}[{i}]")
    elif isinstance(expected, float):
        check(isinstance(actual, (int, float)) and np.isfinite(actual)
              and np.isclose(actual, expected, rtol=1e-10, atol=5e-9), path)
    else:
        check(type(actual) is type(expected) and actual == expected, path)


def local_probe(h, target, bins=BINS, step=STEP):
    """Finite Fourier average of a far-dimer swap; no source state enters."""
    energies, vectors = eigh(h)
    local_x = np.kron(np.eye(4), [[0, 1], [1, 0]])
    phases = step * (energies[:, None] - energies[None, :] + target)
    weight = (np.exp(0.5j * (bins - 1) * phases)
              * np.sinc(bins * phases / (2 * np.pi)) / np.sinc(phases / (2 * np.pi)))
    kraus = weight * (vectors.T @ local_x @ vectors)
    effect = vectors @ (kraus.conj().T @ kraus) @ vectors.T
    check(np.linalg.eigvalsh(effect)[0] > -2e-13
          and np.linalg.eigvalsh(effect)[-1] <= 1 + 2e-13, "probe is not an effect")
    return effect, vectors @ kraus @ vectors.T


def group_spectrum(local, joint, phase):
    """Group equal energy eigenvalues; no eigenvector inside a degenerate level is read."""
    grouped = [np.unique(np.round(levels - levels[0], 9), return_inverse=True)
               for levels, _ in local]
    weights, labels = [], []
    for il, ir in product(range(len(grouped[0][0])), range(len(grouped[1][0]))):
        mask = ((grouped[0][1] == il)[:, None] & (grouped[1][1] == ir)[None, :])
        weights.append(float(joint[mask].sum()))
        labels.append(dict(ring_phase=phase, left_energy=float(grouped[0][0][il]),
                           right_energy=float(grouped[1][0][ir])))
    return weights, labels


def reference():
    _, anchors, _ = wall_flight.extended_geometry()
    ref = walls.reference(0)
    boundary = np.asarray(walls.loads(ref))[anchors]
    ground = np.outer(*(walls.wing_ref(8, int(b), False, 3)[2] for b in boundary))
    basis = list(product((0, 1), repeat=4))
    states = np.zeros((16, 8, 8, 2))
    for case in (0, 1):
        source = tuple(bit ^ int(case and j >= 2) for j, bit in enumerate(ref))
        states[basis.index(source), :, :, case] = ground
    order = [basis.index(ref), basis.index(tuple(1 - x for x in ref))]
    return boundary, states.reshape(1024, 2), order


def run():
    rng = np.random.default_rng(SEED)
    boundary, initial, order = reference()
    native = ex.wing(8, 3, int(boundary[0]))
    e, _ = eigh(native)
    targets = (e[1:4] - e[0]).tolist()
    check(all(16 < x < 48 for x in targets), "native targets leave the hard band")
    calibrations = []
    for b in (0, 2):
        h = ex.wing(8, 3, b)
        energies, vectors = eigh(h)
        readings = []
        for target in targets:
            effect, _ = local_probe(h, target)
            response = np.diag(vectors.T @ effect @ vectors).real
            check(max(response[1:4]) > 0.20 and response[0] < 1e-6,
                  "hard-state probe sensitivity or ground leakage fails")
            readings.append(dict(target=float(target), responses_by_level=response.tolist(),
                                 calibration_hit_counts=[int(rng.binomial(SHOTS, x)) for x in response]))
        calibrations.append(dict(boundary=b, levels=(energies - energies[0]).tolist(), settings=readings))
    arms = []
    comparisons = 3 * len(TIMES) * len(targets)**2
    halfwidth = float(np.sqrt(np.log(2 * comparisons / 0.01) / SHOTS))
    for k, cut in ARMS:
        h = csr_matrix(walls.operators(k, cut, 3))
        centered = h - eye(h.shape[0]) * np.mean(h.diagonal())
        flights = [expm_multiply(-1j * t * centered, initial).reshape(16, 8, 8, 2) for t in TIMES]
        panels = []
        for t, flight in zip(TIMES, flights, strict=True):
            spectra, detectors = [], []
            for source in (0, 1):
                found = flight[order, :, :, source]
                success = float(np.vdot(found, found).real)
                joint_weights, spectrum_labels, phase_energies, phase_success = [], [], [], []
                detector_probabilities = np.zeros((2, 3, 3, 4))
                conditional_energy = np.zeros(2)
                for phase, state in enumerate(found):
                    hs = [walls.wing_ref(k, int(b if phase == 0 else 2 - b), cut, 3)[0]
                          for b in boundary]
                    local = [eigh(wing) for wing in hs]
                    spectra_by_phase = abs(local[0][1].T @ state @ local[1][1])**2
                    weights, labels = group_spectrum(local, spectra_by_phase, phase)
                    joint_weights.extend(weights)
                    spectrum_labels.extend(labels)
                    norm = float(spectra_by_phase.sum())
                    local_energy = [float((levels - levels[0]) @
                                          spectra_by_phase.sum(axis=1 - side))
                                    for side, (levels, _) in enumerate(local)]
                    conditional_energy += local_energy
                    phase_success.append(norm)
                    phase_energies.append((np.asarray(local_energy) / norm).tolist())
                    effects = [[local_probe(wing, target)[0] for target in targets] for wing in hs]
                    for left, right in product(range(3), repeat=2):
                        p = ex.read_joint(state, effects[0][left], effects[1][right])
                        check(min(p) >= 0, "joint detector returned a negative branch")
                        detector_probabilities[phase, left, right] = p
                spectra.append(dict(wall_source=bool(source), ring_ordered_probability=success,
                                    ring_phase_probabilities=phase_success,
                                    spectrum_outcome_labels=spectrum_labels,
                                    joint_spectrum=sample(np.r_[joint_weights, 1 - success], rng),
                                    conditional_mean_excitation_by_phase=phase_energies,
                                    conditional_mean_excitation=(conditional_energy / success).tolist()))
                detectors.append([dict(left=i, right=j,
                                       **sample(np.r_[detector_probabilities[:, i, j].ravel(), 1 - success], rng))
                                  for i, j in product(range(3), repeat=2)])
            rates = []
            for sham, wall in zip(detectors[0], detectors[1], strict=True):
                observed = (sum(wall['counts'][i] for i in (3, 7))
                            - sum(sham['counts'][i] for i in (3, 7))) / SHOTS
                exact = (sum(wall['probabilities'][i] for i in (3, 7))
                         - sum(sham['probabilities'][i] for i in (3, 7)))
                rates.append(dict(left=wall['left'], right=wall['right'],
                                  measured_wall_minus_sham_double_click=observed,
                                  exact_wall_minus_sham_double_click=exact,
                                  simultaneous_99_halfwidth=halfwidth,
                                  resolved_positive=bool(observed > halfwidth)))
            cut_factorization_residual = None
            if cut:
                differences = []
                for sham, wall in zip(detectors[0], detectors[1], strict=True):
                    conditional = [
                        np.asarray(record['probabilities'][:-1]).reshape(2, 4).sum(axis=0)
                        / spectra[source]['ring_ordered_probability']
                        for source, record in enumerate((sham, wall))
                    ]
                    differences.append(float(max(abs(conditional[0] - conditional[1]))))
                cut_factorization_residual = max(differences)
                check(cut_factorization_residual < 1e-10, "cut guide state depends on wall source")
            panels.append(dict(exposure=t, spectra=spectra, detector_records=detectors,
                               cut_conditional_detector_factorization_residual=cut_factorization_residual,
                               source_excess=rates))
        arms.append(dict(contrast=k, attachments_cut=cut, readings=panels))
    return dict(source="Existing guarded phase-wall ring and two three-dimer wings",
                shots=SHOTS, seed=SEED, clock_bins=BINS, clock_step=STEP,
                target_frequencies=targets, calibrations=calibrations, arms=arms,
                spectrum_outcomes="Listed grouped energies and ring phase; source miss last",
                detector_outcomes="phase0(neither,left only,right only,both),phase1(same),source miss",
                scope="Local downward material-transition response; radiation identity unresolved")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--receipt', type=Path, default=HERE / 'results.json')
    args = parser.parse_args()
    report = run()
    if args.output:
        with args.output.open('x') as stream:
            json.dump(report, stream, indent=2, allow_nan=False)
            stream.write('\n')
    else:
        compare(report, json.loads(args.receipt.read_text()), 'daughter_spectrum')
        print('Replayed far-dimer absorption on the native phase-wall source and all controls.')
