"""Clock records and mass-response checks for the unchanged mobile-memory source.

The cage clock is an addressed, prepared source-time device. After a stored
current is fixed, the charge and load-cage operators act on disjoint factors.
A coherent mass-offset tag is a supplied Ramsey diagnostic. With one charge
the offset is a multiple of the identity, including along every branch of
any mass-independent schedule and instrument. This argument excludes carrier
conversion sectors in which the charged occupation changes.
"""

import numpy as np
from scipy import sparse

from model.executable_model import carrier
from model.experiments.cage.construction import ENERGY
from model.experiments.radiative_repair.source import evolve
from ..cage_transfer import PORTS_PER_CELL, transfer
from ..law import MASS, charged_operator
from ..mobility import home_seam, lowest_band_preparation


SOURCE_TICK = 0.5
TRAINING = 8
QUERIES = np.array([*range(TRAINING), 47, 96])
WAITS = QUERIES * SOURCE_TICK
CAGE_LINK = 1.0
MASS_OFFSET = 0.25
TOLERANCE = 2e-10


def _clock_reader(full_fuel):
    """Recover a line from training clicks, then predict separate waits."""
    training = full_fuel[:TRAINING]
    differences = np.diff(training)
    middle = differences[1:-1]
    norm = middle @ middle
    if norm < 1e-12:
        raise ValueError("the training counts have no resolved clock line")
    cosine = middle @ (differences[:-2] + differences[2:]) / (2 * norm)
    if not -1 < cosine < 1:
        raise ValueError("the line is outside the enrolled unaliased interval")
    frequency = np.arccos(cosine)
    phases = frequency * QUERIES
    design = np.column_stack((np.ones(len(QUERIES)), np.cos(phases), np.sin(phases)))
    coefficients = np.linalg.lstsq(design[:TRAINING], training, rcond=None)[0]
    error = float(np.max(abs(design @ coefficients - full_fuel)))
    if error > TOLERANCE:
        raise ValueError("the clock line does not predict the held-out waits")
    return float(frequency), error


def _cell_records(states, cells):
    return (abs(states) ** 2).reshape(-1, cells, PORTS_PER_CELL).sum(axis=-1)


def _native_clock(car, seam, initial, active, parked):
    """Read the existing finite fuel pulse and compare the full product records."""
    written = [transfer(seam=seam, memory=True, coupling=CAGE_LINK, wait=t) for t in WAITS]
    blank = [transfer(seam=seam, memory=False, coupling=CAGE_LINK, wait=t) for t in WAITS]
    fuel = np.array([row["receiver_record"] for row in written])
    controls = {
        "no_sender": np.array([row["reference_record"] for row in written]),
        "blank_memory": np.array([row["receiver_record"] for row in blank]),
        "no_sender_blank_memory": np.array([row["reference_record"] for row in blank]),
    }
    for stationary in controls.values():
        np.testing.assert_allclose(stationary, np.broadcast_to(stationary[0], stationary.shape),
                                   atol=TOLERANCE, rtol=0)
    ideal = np.array([row["exciton_record"][3] for row in written])
    np.testing.assert_allclose(ideal, np.sin(CAGE_LINK * WAITS) ** 2, atol=TOLERANCE, rtol=0)
    observed, error = _clock_reader(fuel[:, 3])
    np.testing.assert_allclose(observed, 2 * CAGE_LINK * SOURCE_TICK, atol=TOLERANCE, rtol=0)

    receiver = seam[1] // PORTS_PER_CELL
    acceptance = float(np.vdot(initial, initial).real)
    cage_start = [1.0, 0.0]
    cage_source = sparse.csr_matrix(2 * ENERGY * np.eye(2) + CAGE_LINK * np.array([[1, -1], [-1, 1]]))
    arms = {}
    for label, source in (("traveling", active), ("parked_hops", parked)):
        charge_cells = _cell_records(np.array(list(evolve(source, initial, WAITS))), car.cells)
        joint = charge_cells[:, :, None] * fuel[:, None, :]
        np.testing.assert_allclose(joint.sum(axis=(1, 2)), acceptance, atol=TOLERANCE, rtol=0)
        receiver_fuel = joint[:, receiver] / charge_cells[:, receiver, None]
        rate, heldout = _clock_reader(receiver_fuel[:, 3])

        # Independent sparse propagation on actual carrier charge ports times
        # the two exact cage positions, followed by joint occupation reads.
        tensor_source = sparse.kronsum(cage_source, source, format="csr")
        tensor_start = np.kron(initial, cage_start)
        propagated = np.array(list(evolve(tensor_source, tensor_start, WAITS)))
        tensor_records = (abs(propagated) ** 2).reshape(-1, car.cells, PORTS_PER_CELL, 2).sum(axis=2)
        ideal_joint = charge_cells[:, :, None] * np.column_stack((1 - ideal, ideal))[:, None, :]
        tensor_residual = float(np.max(abs(tensor_records - ideal_joint)))
        if tensor_residual > TOLERANCE:
            raise ValueError("direct tensor source did not give the joint occupation records")
        arms[label] = {
            "charge_cell_probabilities": charge_cells.tolist(),
            "joint_charge_cell_and_fuel_probabilities": joint.tolist(),
            "receiver_conditioned_fuel_probabilities": receiver_fuel.tolist(),
            "read_radians_per_tick": rate,
            "clock_fit_max_error": heldout,
            "ideal_joint_tensor_residual": tensor_residual,
        }

    gain = arms["traveling"]["charge_cell_probabilities"][-1][receiver] - arms["parked_hops"]["charge_cell_probabilities"][-1][receiver]
    ratio = arms["traveling"]["read_radians_per_tick"] / arms["parked_hops"]["read_radians_per_tick"]
    np.testing.assert_allclose(ratio, 1, atol=TOLERANCE, rtol=0)
    if gain < 0.1:
        raise ValueError("the charge control did not resolve transport")
    return {
        "receiver_cell": receiver,
        "filter_rejection_probability": 1 - acceptance,
        "fuel_bits": ["00", "01", "10", "11"],
        "fuel_probabilities": fuel.tolist(),
        "control_probabilities": {name: rows.tolist() for name, rows in controls.items()},
        "ideal_cage_receiver_probability": ideal.tolist(),
        "read_radians_per_tick": observed,
        "clock_fit_max_error": error,
        "arms": arms,
        "traveling_minus_parked_charge_at_final_receiver": gain,
        "traveling_over_parked_clock_rate": ratio,
    }


def _expected_analyzer(phase):
    contrast = np.column_stack((np.cos(phase), np.sin(phase)))
    return np.stack(((1 + contrast) / 2, (1 - contrast) / 2), axis=-1)


def _tagged_clock(source, initial, cells):
    """Propagate both tagged mass arms and retain both analyzer outcomes."""
    identity = sparse.eye(source.shape[0], format="csr")
    tag_source = sparse.block_diag((source, source + MASS_OFFSET * identity), format="csr")
    stock = np.concatenate((initial, initial)) / np.sqrt(2)
    states = np.array(list(evolve(tag_source, stock, WAITS)))
    left, right = np.split(states, 2, axis=-1)
    analyzers = ((left + right, left - right), (left + 1j * right, left - 1j * right))
    records = np.stack([
        np.stack([_cell_records(vector / np.sqrt(2), cells) for vector in pair], axis=-1)
        for pair in analyzers
    ], axis=2)
    cells_at_read = records.sum(axis=-1)[:, :, 0]
    expected = cells_at_read[:, :, None, None] * _expected_analyzer(MASS_OFFSET * WAITS)[:, None, :, :]
    np.testing.assert_allclose(records, expected, atol=TOLERANCE, rtol=0)

    acceptance = float(np.vdot(initial, initial).real)
    np.testing.assert_allclose(records.sum(axis=(1, 3)), acceptance, atol=TOLERANCE, rtol=0)
    contrast = (records[..., 0] - records[..., 1]).sum(axis=1) / acceptance
    phase = np.unwrap(np.angle(contrast[:TRAINING, 0] + 1j * contrast[:TRAINING, 1]))
    rate = float((phase[-1] - phase[0]) / (TRAINING - 1))
    return {
        "joint_cell_analyzer_outcome_probabilities": records.tolist(),
        "received_clock_contrast": contrast.tolist(),
        "read_radians_per_tick": rate,
        "filter_rejection_probability": 1 - acceptance,
    }


def _mass_clock(car, source, traveling):
    rest = np.ones(car.n) / np.sqrt(car.n)
    a = _tagged_clock(source, rest, car.cells)
    b = _tagged_clock(source, traveling, car.cells)
    ratio = b["read_radians_per_tick"] / a["read_radians_per_tick"]
    np.testing.assert_allclose(ratio, 1, atol=TOLERANCE, rtol=0)
    np.testing.assert_allclose(a["read_radians_per_tick"], MASS_OFFSET * SOURCE_TICK,
                               atol=TOLERANCE, rtol=0)

    # Same positive kinetic eigenmode in both mass arms of a spectral
    # square-root comparator. Its kinetic term 2*MASS*K is held fixed
    # under the offset, matching the original slope at rest mass MASS.
    energy = np.linalg.eigvalsh(source.toarray())[car.cells - 1]
    kinetic = float(energy - MASS)
    if kinetic <= 0:
        raise ValueError("the enrolled low-band spectral mode has no positive kinetic energy")
    momentum_term = 2 * MASS * kinetic
    benchmark_ratio = (2 * MASS + MASS_OFFSET) / (
        np.sqrt((MASS + MASS_OFFSET) ** 2 + momentum_term) + np.sqrt(MASS ** 2 + momentum_term)
    )
    if benchmark_ratio >= ratio - 1e-3:
        raise ValueError("the matched positive mode does not separate the two clock laws")
    return {
        "analyzer_settings": ["X", "negative_Y"],
        "analyzer_outcomes": [1, -1],
        "onsite_mass_offset": MASS_OFFSET,
        "rest": a,
        "traveling": b,
        "ratio": ratio,
        "benchmark_kinetic_eigenvalue": kinetic,
        "benchmark_fixed_squared_energy_term": momentum_term,
        "matched_spectral_mode_square_root_ratio": float(benchmark_ratio),
        "matched_spectral_mode_benchmark_probabilities": _expected_analyzer(
            MASS_OFFSET * benchmark_ratio * WAITS
        ).tolist(),
    }


def run():
    """Replay the addressed source-time clock and the mass-offset audit."""
    car = carrier.federation(0)
    seam = home_seam(car)
    initial = lowest_band_preparation(car)
    active = charged_operator(car)
    parked = charged_operator(car, parked=True)
    return {
        "oriented_seam": list(seam),
        "source_tick": SOURCE_TICK,
        "query_ticks": QUERIES.tolist(),
        "training_ticks": QUERIES[:TRAINING].tolist(),
        "native_cage_clock": _native_clock(car, seam, initial, active, parked),
        "mass_offset_ramsey": _mass_clock(car, active, initial),
    }
