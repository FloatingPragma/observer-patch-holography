"""Account for passive energy uptake by supplied local two-level receivers."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix, eye
from scipy.sparse.linalg import expm_multiply

from .. import experiment as ex, walls
from ..daughter_spectrum import measure as parent

HERE = Path(__file__).resolve().parent
DELAYS = (0, 16, 32, 64)
COUPLING = 1 / 16
SEED = 81168
SHOTS = parent.SHOTS
X = np.array([[0., 1.], [1., 0.]])
N = np.diag([0., 1.])
FAR = np.kron(np.eye(4), X)


def check(condition, message):
    if not condition:
        raise ArithmeticError(message)


def receiver(h, omega, delay, coupling=COUPLING):
    """Wing-major tensor order; only the far dimer couples to the detector."""
    levels, vectors = eigh(h)
    free_wing = np.kron(h - levels[0] * np.eye(8), np.eye(2))
    free_detector = np.kron(np.eye(8), omega * N)
    interaction = coupling * np.kron(FAR, X)
    generator = free_wing + free_detector + interaction
    energy, modes = eigh(generator)
    unitary = (modes * np.exp(-1j * delay * energy)) @ modes.T
    kraus = [unitary.reshape(8, 2, 8, 2)[:, d, :, 0] for d in (0, 1)]
    return dict(levels=levels, vectors=vectors, unitary=unitary, kraus=kraus,
                operators=(free_wing, free_detector, interaction), generator=generator)


def ledger(a, apparatus):
    """Unnormalized branch energy expectations; no extra source postselection."""
    total = np.zeros(5)
    for side, meter in enumerate(apparatus):
        rho = a @ a.conj().T if side == 0 else a.T @ a.conj()
        initial = np.kron(rho, np.diag([1., 0.]))
        final = meter['unitary'] @ initial @ meter['unitary'].conj().T
        wing, detector, interaction = meter['operators']
        wing_loss = np.trace(wing @ (initial - final)).real
        detector_gain = np.trace(detector @ final).real
        residual = np.trace(interaction @ final).real
        drift = np.trace(meter['generator'] @ (final - initial)).real
        total += [wing_loss, detector_gain, residual, drift,
                  np.trace(wing @ initial).real]
    check(abs(total[3]) < 2e-9, 'receiver total energy drift')
    check(abs(total[0] - total[1] - total[2]) < 2e-9, 'receiver energy ledger')
    return total


def run():
    rng = np.random.default_rng(SEED)
    boundary, initial, ordered = parent.reference()
    levels = eigh(ex.wing(8, 3, int(boundary[0])), eigvals_only=True)
    omega = float(levels[2] - levels[0])
    calibrations = []
    for k, cut in parent.ARMS:
        for b in (0, 2):
            h = walls.wing_ref(k, b, cut, 3)[0]
            settings = []
            for delay in DELAYS:
                meter = receiver(h, omega, delay)
                responses = []
                for vector in meter['vectors'].T:
                    final = meter['unitary'] @ np.kron(vector, [1., 0.])
                    responses.append(float(np.sum(abs(final.reshape(8, 2)[:, 1])**2)))
                settings.append(dict(delay=delay, click_probabilities_by_level=responses))
            calibrations.append(dict(contrast=k, cut=cut, boundary=b,
                                     excitation_levels=(eigh(h, eigvals_only=True) - eigh(h, eigvals_only=True)[0]).tolist(),
                                     settings=settings))
    arms = []
    for k, cut in parent.ARMS:
        h = csr_matrix(walls.operators(k, cut, 3))
        centered = h - eye(1024) * np.mean(h.diagonal())
        # Pinch the ring's complete sharp read, including every missed source word.
        coo = h.tocoo()
        pinched = csr_matrix((coo.data * (coo.row // 64 == coo.col // 64),
                              (coo.row, coo.col)), shape=h.shape)
        panels = []
        for exposure in parent.TIMES:
            flight = expm_multiply(-1j * exposure * centered, initial).reshape(16, 8, 8, 2)
            cohorts = []
            for source in (0, 1):
                full = flight[:, :, :, source].ravel()
                energy_initial = float(np.vdot(initial[:, source], h @ initial[:, source]).real)
                energy_before = float(np.vdot(full, h @ full).real)
                energy_after = float(np.vdot(full, pinched @ full).real)
                check(abs(energy_before - energy_initial) < 2e-7, 'free source energy drift')
                found = flight[ordered, :, :, source]
                success = float(np.vdot(found, found).real)
                records = []
                for delay in DELAYS:
                    weights, energies = [], np.zeros(5)
                    for phase, a in enumerate(found):
                        apparatus = [receiver(walls.wing_ref(k, int(b if phase == 0 else 2 - b), cut, 3)[0],
                                              omega, delay) for b in boundary]
                        energies += ledger(a, apparatus)
                        for left, right in ((0, 0), (1, 0), (0, 1), (1, 1)):
                            state = apparatus[0]['kraus'][left] @ a @ apparatus[1]['kraus'][right].T
                            weights.append(float(np.vdot(state, state).real))
                    weights.append(float(np.sum(abs(flight[[j for j in range(16) if j not in ordered], :, :, source])**2)))
                    record = parent.sample(weights, rng)
                    record.update(delay=delay,
                                  energy_per_original_launch=dict(zip(
                                      ('wing_loss', 'detector_gain', 'interaction_final', 'total_drift', 'initial_wing_excitation'),
                                      [float(x) for x in energies], strict=True)),
                                  detector_gain_per_ordered_source=float(energies[1] / success),
                                  switch_off_work_per_original_launch=float(-energies[2]))
                    records.append(record)
                cohorts.append(dict(source=('ordered_sham', 'wall_pair')[source],
                                    source_order_probability=success,
                                    before_ring_read_energy=energy_before,
                                    initial_total_energy=energy_initial,
                                    ring_read_and_hold_work=energy_after - energy_before,
                                    energy_after_complete_ring_read=energy_after,
                                    records=records))
            panels.append(dict(exposure=exposure, cohorts=cohorts))
        arms.append(dict(contrast=k, attachments_cut=cut, readings=panels))
    native = arms[0]['readings'][1]['cohorts']
    at32 = [c['records'][2] for c in native]
    both = lambda r: r['probabilities'][3] + r['probabilities'][7]
    difference = both(at32[1]) - both(at32[0])
    counted_both = lambda r: r['counts'][3] + r['counts'][7]
    measured_difference = (counted_both(at32[1]) - counted_both(at32[0])) / SHOTS
    comparisons = len(parent.ARMS) * len(parent.TIMES) * len(DELAYS)
    width = float(np.sqrt(np.log(2 * comparisons / .01) / SHOTS))
    check(measured_difference > width, 'native wall-minus-sham coincidence unresolved')
    inputs = [Path(__file__), Path(parent.__file__), Path(ex.__file__), Path(walls.__file__)]
    return dict(shots=SHOTS, seed=SEED, target_frequency=omega, coupling=COUPLING,
                outcome_order='phase0(neither,left,right,both),phase1(same),source miss',
                calibration=calibrations, arms=arms,
                native_late_double_click_excess=difference,
                measured_native_late_double_click_excess=measured_difference,
                simultaneous_99_difference_halfwidth=width,
                input_sha256={p.name if p == Path(__file__) else str(p.relative_to(HERE.parent)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                scope='Supplied passive detector and protected source; material daughter energy transfer, radiation identity open')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--receipt', type=Path, default=HERE / 'results.json')
    args = parser.parse_args()
    report = run()
    serialized = json.dumps(report, indent=2, allow_nan=False) + '\n'
    if args.output:
        with args.output.open('x') as stream:
            stream.write(serialized)
    else:
        parent.compare(report, json.loads(args.receipt.read_text()), 'passive_receiver')
        print('Verified static local energy receivers, complete source cohorts and work ledgers.')


if __name__ == '__main__':
    main()
