"""Prepare a written-core orbital by a finite energy tag and test native erasure."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from model.experiments import count_intervals
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo import candidate
from model.executable_model import EMITTER_SEAM, federation
from ..written_transitions import measure as transition

HERE = Path(__file__).resolve().parent
SHOTS, SEED = 2**22, 317031
TARGET = -6.136
TIMES = (0, 1, 4, 16, 64)
OPERATIONS = ('retain', 'one_native_echo', 'two_native_echoes', 'supplied_memory_erase')


def check(value, message):
    if not value:
        raise ArithmeticError(message)


def sample(raw, rng):
    p = np.asarray(raw, float)
    check(np.all(np.isfinite(p)) and np.min(p) >= 0, 'Invalid terminal probability')
    check(abs(p.sum() - 1) < 2e-10, 'Incomplete terminal distribution')
    order = np.argsort(p, kind='stable')
    positive = p[order].copy()
    positive[-1] = 1 - positive[:-1].sum()
    check(positive[-1] >= 0, 'Negative explicit complement')
    counts = np.empty(len(p), dtype=np.int64)
    counts[order] = rng.multinomial(SHOTS, positive)
    return dict(probabilities=p.tolist(), counts=counts.tolist(),
                sampling_complement_l1=float(abs(positive[-1] - p[order[-1]])))


def native_certificate():
    law = candidate.native_law
    echo = law.echo()
    vacuum = law.SRC.index((0, 0))
    blank = np.eye(law.DIM)[:, vacuum]
    written = echo @ blank
    record = written.reshape(-1, law.N_SRC)[:, vacuum]
    h0, hc = law.local_term_parts()
    h = h0 + candidate.KAPPA * hc
    errors = []
    erasure = []
    for word in ((0, 0), (1, 0), (0, 1)):
        index = law.SRC.index(word)
        state = np.kron(record, np.eye(law.N_SRC)[:, index])
        expected = np.eye(law.DIM)[:, index] if word == (0, 0) else state
        errors.extend([np.linalg.norm(echo @ state - expected), np.linalg.norm(h @ state)])
        erasure.append(float(np.sum(abs((echo @ state).reshape(-1, law.N_SRC)[0])**2)))
    check(max(errors) < 1e-12, 'Native memory or zero-action reduction fails')
    check(np.linalg.norm(echo @ echo - np.eye(law.DIM)) < 1e-12, 'Echo is not involutive')
    return dict(formation_input='blank memory and 00 endpoints',
                formation_probabilities_by_memory=abs(record)**2,
                formation_work=float(np.vdot(written, h @ written).real - np.vdot(blank, h @ blank).real),
                blank_probability_after_echo_for_00_10_01=erasure,
                maximum_local_identity_residual=float(max(errors)))


def prepare(h, port):
    e, v = eigh(h)
    sharp = v[port].copy()
    kernel = transition.kernel(e - TARGET).conj()
    transmitted = kernel * sharp
    reflected = np.sqrt(1 - abs(kernel)**2) * sharp
    psi = v @ transmitted
    success = float(np.vdot(psi, psi).real)
    failure = float(np.vdot(reflected, reflected).real)
    initial_mean = float(h[port, port])
    admitted_energy = float(np.sum(e * abs(transmitted)**2))
    refused_energy = float(np.sum(e * abs(reflected)**2))
    check(abs(success + failure - 1) < 2e-12, 'Source tag does not conserve probability')
    check(abs(admitted_energy + refused_energy - initial_mean) < 2e-11, 'Source tag energy ledger fails')
    return psi, failure, dict(target=TARGET, bins=transition.BINS, step=transition.STEP,
        acceptance=success, refusal=failure, conditional_ground_fidelity=float(abs(transmitted[0])**2 / success),
        initial_relative_energy=initial_mean, admitted_relative_energy=admitted_energy,
        refused_relative_energy=refused_energy,
        conditional_relative_energy=admitted_energy / success,
        largest_bound_line_weight=float(abs(transmitted[0])**2))


def branch_states(psi, seam, operation):
    blank, written = np.zeros_like(psi), psi.copy()
    if operation == 'one_native_echo':
        blank = psi.copy()
        blank[list(seam)] = 0
        written = psi - blank
    elif operation == 'supplied_memory_erase':
        blank, written = psi.copy(), np.zeros_like(psi)
    elif operation not in ('retain', 'two_native_echoes'):
        raise ValueError(operation)
    return blank, written


def evolve(eigen, state, t):
    e, v = eigen
    return v @ (np.exp(-1j * t * e) * (v.T @ state))


def panel(level, rng):
    board = federation(level)
    seam = tuple(EMITTER_SEAM[level])
    h = transition.generator(board, seam, candidate.KAPPA)
    psi, failure, preparation = prepare(h, seam[0])
    success = preparation['acceptance']
    distances = shortest_path(board.adjacency(), directed=False, unweighted=True,
                              indices=list(seam)).min(axis=0).astype(int)
    r = native_certificate()['formation_probabilities_by_memory']
    arms = []
    for contrast in (candidate.KAPPA, 0):
        blank_h = board.laplacian().toarray() / 2
        written_h = transition.generator(board, seam, contrast)
        blank_spectrum, written_spectrum = eigh(blank_h), eigh(written_h)
        law_quench = float(np.vdot(psi, (written_h - h) @ psi).real)
        cases = []
        for operation in OPERATIONS:
            blank, written = branch_states(psi, seam, operation)
            after_energy = float((np.vdot(blank, blank_h @ blank) + np.vdot(written, written_h @ written)).real)
            before_energy = float(np.vdot(psi, written_h @ psi).real)
            blank_mass = float(np.vdot(blank, blank).real)
            readings = []
            for wait in TIMES:
                a, b = evolve(blank_spectrum, blank, wait), evolve(written_spectrum, written, wait)
                joint = r[:, None] * abs(b)[None, :]**2
                joint[0] += abs(a)**2
                reading = sample(np.r_[joint.ravel(), failure], rng)
                coarse = abs(a)**2 + abs(b)**2
                port_counts = np.array(reading['counts'][:-1]).reshape(11, board.n).sum(axis=0)
                admitted = int(port_counts.sum())
                core_count = int(port_counts[list(seam)].sum())
                tail = .01 / (2 * 2 * 2 * len(OPERATIONS) * len(TIMES))
                lower, upper = count_intervals.clopper_pearson(core_count, admitted, tail)
                interval = [float(lower), float(upper)]
                public_probabilities = np.r_[joint.sum(axis=0), failure]
                public_counts = np.r_[port_counts, reading['counts'][-1]]
                final_energy = float((np.vdot(a, blank_h @ a) + np.vdot(b, written_h @ b)).real)
                check(abs(final_energy - after_energy) < 2e-10, 'Released flight energy drift')
                reading.update(wait=wait, core_probability_per_start=float(coarse[list(seam)].sum()),
                    public_probabilities_by_port_then_source_refusal=public_probabilities.tolist(),
                    public_counts_by_port_then_source_refusal=public_counts.tolist(),
                    measured_core_given_admission=core_count / admitted,
                    simultaneous_99_core_given_admission_interval=interval,
                    measured_conditional_rms_radius=float(np.sqrt(np.dot(distances**2, port_counts) / admitted)),
                    blank_memory_probability_per_start=float(np.sum(abs(a)**2)),
                    far_probability_per_start=float(coarse[distances >= 3].sum()),
                    conditional_rms_radius=float(np.sqrt(np.dot(distances**2, coarse) / success)),
                    released_relative_energy=final_energy)
                readings.append(reading)
            cases.append(dict(operation=operation,
                              blank_memory_probability_before_flight=blank_mass,
                              intervention_work_per_start=after_energy - before_energy,
                              readings=readings))
        arms.append(dict(contrast=contrast,
                         contrast_quench_relative_work_per_start=law_quench,
                         contrast_quench_work_per_start=law_quench + 6 * (contrast - candidate.KAPPA) * float(np.vdot(psi, psi).real),
                         cases=cases))
    e, v = eigh(h)
    transitions = []
    for label, ports, target in (('endpoint_native_line', [seam[0]], .47993),
                                  ('common_phase_native_line', seam, .47993),
                                  ('hydrogen_target', [seam[0]], .75 * abs(TARGET))):
        filtered = transition.filtered(e, v, psi, ports, target)
        weight = float(np.vdot(filtered, filtered).real)
        record = sample(np.r_[abs(filtered)**2, success - weight, failure], rng)
        record.update(setting=label, frequency=target, acceptance_per_start=weight)
        transitions.append(record)
    return dict(cells=board.cells, ports=board.n, seam=list(seam),
                radius_by_port=distances.tolist(), preparation=preparation,
                arms=arms, transition_records=transitions)


def run():
    rng = np.random.default_rng(SEED)
    certificate = native_certificate()
    certificate['formation_probabilities_by_memory'] = certificate['formation_probabilities_by_memory'].tolist()
    formation_calibration = sample(certificate['formation_probabilities_by_memory'], rng)
    panels = [panel(level, rng) for level in (0, 1)]
    paths = [Path(__file__), Path(transition.__file__), Path(candidate.native_law.__file__)]
    return dict(shots=SHOTS, seed=SEED, native_certificate=certificate,
                independent_formation_memory_diagnostics=formation_calibration,
                native_uniform_single_load_reference=6 * candidate.KAPPA,
                memory_labels=['blank'] + [list(x) for x in candidate.native_law.LABELS], panels=panels,
                diagnostic_outcome_order='Private collar diagnostics: 11 native memory labels, each all ports; source-tag refusal last',
                public_outcome_order='All terminal load ports and source-tag refusal; collar probabilities require a separate public readback',
                transition_outcome_order='all ports, transition-clock refusal, source-tag refusal',
                code_sha256={str(p.relative_to(HERE.parents[3])): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                scope='Stationary native memory core with a clock-tagged load; no mobile nucleus or hydrogenic shells')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--receipt', type=Path, default=HERE / 'results.json.gz')
    args = parser.parse_args()
    data = run()
    text = json.dumps(data, indent=1, allow_nan=False) + '\n'
    if args.output:
        encoded = gzip.compress(text.encode(), mtime=0) if args.output.suffix == '.gz' else text.encode()
        with args.output.open('xb') as stream:
            stream.write(encoded)
    else:
        archived = gzip.decompress(args.receipt.read_bytes()).decode() if args.receipt.suffix == '.gz' else args.receipt.read_text()
        transition.compare(data, json.loads(archived), 'native_core_release')
        print('Verified native source formation, finite orbital preparation, erasure and released records.')


if __name__ == '__main__':
    main()
