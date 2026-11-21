"""Co-active charge, cage and native load records on a fixed incidence metric."""

from dataclasses import dataclass
import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import law as native_law
from ..contrast_binding import _native_controls
from ..mobility import PORTS_PER_CELL, home_seam
from ..relativity.incidence import (
    CAGE_MASS, CHARGE_MASS, cell_incidence, dirac_operator, positive_embedding,
)


KAPPA = 128
TOLERANCE = 2e-10


@dataclass(frozen=True)
class Source:
    hamiltonian: sparse.csr_matrix
    incidence: sparse.csr_matrix
    load_operator: np.ndarray
    charge_operator: sparse.csr_matrix
    cage_operator: sparse.csr_matrix
    initial: np.ndarray
    charge_rest: np.ndarray
    cage_rest: np.ndarray
    charge_grading: np.ndarray
    cage_grading: np.ndarray
    shape: tuple
    cell_ids: tuple
    seam_ids: tuple
    addressed_seam: tuple
    pairs: tuple
    swap_seam: tuple
    contact_index: int
    lengths: tuple
    length: float
    speed: float
    coupling: float
    charge_mass: float
    cage_mass: float


def _native_pair_chart(car, addressed, parked_swap):
    """Hold every load move except one native swap inside an actual face."""
    sender = addressed[0]
    face = next(face for face in car.faces() if sender in face)
    anchor, alternate = [port for port in face if port != sender]
    pairs = (tuple(sorted((anchor, sender))), tuple(sorted((anchor, alternate))))
    swap = tuple(sorted((sender, alternate)))

    exchange, _ = native_law.local_term_parts()
    rows = [native_law.SRC.index(tuple(int(port in pair) for port in swap)) for pair in pairs]
    load = exchange[np.ix_(rows, rows)]
    np.testing.assert_allclose(load, [[0.5, -0.5], [-0.5, 0.5]], atol=TOLERANCE, rtol=0)

    blank_weights, _ = _native_controls(KAPPA)
    weights = np.array([blank_weights[native_law.SRC.index(tuple(int(port in pair) for port in addressed))]
                        for pair in pairs])
    np.testing.assert_array_equal(weights, [1, 0])
    # The common internal-cell contrast cost is removed. The external square
    # appears as the contact activation, with no external blank-energy bias.
    if parked_swap:
        load = np.diag(np.diag(load))
    return load, weights, pairs, swap


def build(car, lengths, speed, *, selected_seam_only=False, coupling=1.0, parked_swap=False,
          charge_mass=CHARGE_MASS, cage_mass=CAGE_MASS):
    """Build the same calibrated generator on a full or addressed chart.

    The source's only load transition is the supplied retained native swap.
    Its contrast record multiplies the charge–cage contact on the addressed
    seam register. The calibrated incidence itself is independent of that
    record, including in the coupling-zero and parked-swap controls.
    """
    if not math.isfinite(coupling) or coupling < 0:
        raise ValueError("contact strength must be finite and nonnegative")

    addressed = home_seam(car)
    load, activations, pairs, swap = _native_pair_chart(car, addressed, parked_swap)
    selected = (addressed,) if selected_seam_only else tuple(car.external)
    if any(seam not in lengths or not math.isfinite(lengths[seam]) or lengths[seam] <= 0 for seam in selected):
        raise ValueError("every enabled seam needs a finite positive calibrated length")
    incidence, seams = cell_incidence(car, {seam: lengths[seam] ** -2 for seam in selected})
    cell_ids = tuple(port // PORTS_PER_CELL for port in addressed) if selected_seam_only else tuple(range(car.cells))
    if selected_seam_only:
        incidence = incidence[list(cell_ids), :][:, [seams.index(addressed)]]
        seams = selected

    cells, edges = incidence.shape
    mode_count = cells + edges
    edge_mode = cells + seams.index(addressed)
    D_charge = dirac_operator(incidence, charge_mass, speed)
    D_cage = dirac_operator(incidence, cage_mass, speed)
    pair_source = sparse.kronsum(D_cage, D_charge, format="csr")
    H = sparse.kron(sparse.csr_matrix(load), sparse.eye(mode_count ** 2), format="csr")
    H += sparse.kron(sparse.eye(2), pair_source, format="csr")
    contact_index = edge_mode * mode_count + edge_mode
    diagonal = np.zeros(2 * mode_count ** 2)
    diagonal[contact_index::mode_count ** 2] = coupling * activations
    H += sparse.diags(diagonal, format="csr")

    sender_cell = cell_ids.index(addressed[0] // PORTS_PER_CELL)
    seed_charge = positive_embedding(incidence, charge_mass, speed)[:, sender_cell]
    seed_cage = positive_embedding(incidence, cage_mass, speed)[:, sender_cell]
    rest = np.r_[np.ones(cells) / math.sqrt(cells), np.zeros(edges)]
    grading = np.r_[np.ones(cells), -np.ones(edges)]
    initial_load = [1, 0]
    shape = (2, mode_count, mode_count)
    return Source(
        hamiltonian=H, incidence=incidence, load_operator=load, charge_operator=D_charge,
        cage_operator=D_cage, initial=np.kron(initial_load, np.kron(seed_charge, seed_cage)),
        charge_rest=np.kron(initial_load, np.kron(rest, seed_cage)),
        cage_rest=np.kron(initial_load, np.kron(seed_charge, rest)),
        charge_grading=np.broadcast_to(grading[None, :, None], shape).ravel(),
        cage_grading=np.broadcast_to(grading[None, None, :], shape).ravel(),
        shape=shape, cell_ids=cell_ids, seam_ids=seams, addressed_seam=addressed, pairs=pairs,
        swap_seam=swap, contact_index=contact_index, lengths=tuple(lengths[seam] for seam in seams),
        length=lengths[addressed], speed=speed,
        coupling=coupling, charge_mass=charge_mass, cage_mass=cage_mass,
    )


def leakage_witness(model):
    """Certify the leading charge-negative probability on the 18-state chart."""
    if model.shape != (2, 3, 3):
        raise ValueError("the exact leakage witness needs the single-seam chart")

    kinetic = math.sqrt(2) * model.speed / model.length
    E_charge = math.hypot(model.charge_mass, kinetic)
    E_cage = math.hypot(model.cage_mass, kinetic)
    # The sender decomposes equally into the symmetric rest mode and the
    # antisymmetric mode. Its seam probability is half the latter's.
    edge_charge = (1 - model.charge_mass / E_charge) / 4
    edge_cage = (1 - model.cage_mass / E_cage) / 4
    leading = model.coupling ** 2 * (1 - 2 * edge_charge) * edge_charge * edge_cage

    frame = positive_embedding(model.incidence, model.charge_mass, model.speed)
    charge_negative = np.eye(3) - frame @ frame.conj().T
    preparation = model.initial.reshape(model.shape)
    edge_reads = [np.linalg.norm(preparation[:, 2, :]) ** 2, np.linalg.norm(preparation[:, :, 2]) ** 2]
    np.testing.assert_allclose(edge_reads, [edge_charge, edge_cage], atol=TOLERANCE, rtol=0)
    projected = np.kron(np.eye(2), np.kron(charge_negative, np.eye(3)))
    np.testing.assert_allclose(projected @ model.initial, 0, atol=TOLERANCE, rtol=0)
    first_order = np.linalg.norm(projected @ (model.hamiltonian @ model.initial)) ** 2
    np.testing.assert_allclose(first_order, leading, atol=TOLERANCE, rtol=0)

    spectral_bound = 1 + E_charge + E_cage + model.coupling
    # Taking imaginary parts of the real source gives an odd Taylor remainder.
    wait = min(1 / 128, math.sqrt(3 * math.sqrt(leading) / spectral_bound ** 3)) if leading else 1 / 128
    lower = max(0, wait * math.sqrt(leading) - (wait * spectral_bound) ** 3 / 6) ** 2
    step = -1j * wait * model.hamiltonian
    evolved = expm_multiply(step, model.initial, traceA=step.diagonal().sum())
    observed = float(np.linalg.norm(projected @ evolved) ** 2)
    if observed < lower - TOLERANCE or (leading > 0 and (lower <= 0 or observed <= 0)):
        raise ValueError("the charge-negative detector failed its certified lower bound")
    return {
        "wait": wait,
        "charge_and_cage_initial_seam_probabilities": [edge_charge, edge_cage],
        "negative_charge_probability_leading_coefficient": leading,
        "negative_charge_probability_lower_bound": lower,
        "negative_charge_probability": observed,
    }
