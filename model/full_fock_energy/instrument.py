"""Energy readout of complete native evolution and a number-pointer instrument.

No Hamiltonian is assigned to a terminated proper-code export. Its probability
is retained as an explicit missing output-energy domain, never normalized out.
"""
from __future__ import annotations

from fractions import Fraction
import math
import numpy as np
from scipy import sparse

from ..field_stress.energy import GaussianInteger, commutator
from ..full_fock_field import field
from ..full_fock_field.bounds import complete_history_bound, stability_bounds
from ..observer.source_scalar import compiler
from ..reference_recovery.recovery import execute_complete
from .energy import local_energies


def expectation(columns, operator):
    return float(np.vdot(columns, operator@columns).real)


def pointer_work(geometry, mode):
    """Certify E*(H)-H=-[n,[n,H]] for the complete binary number read.

    The result is an exact integer numerator divided by energy.scale, not
    a zero-work assumption for an instrument that preserves number.
    """
    energy = local_energies(geometry)
    space = field.FullFock(geometry.modes)
    n = GaussianInteger.from_exact(space.number(mode))
    empty = GaussianInteger.from_exact(sparse.eye(space.dimension))-n
    direct = n@energy.total@n+empty@energy.total@empty-energy.total
    negative_double_commutator = GaussianInteger.zero(space.dimension)-commutator(
        n, commutator(n, energy.total))
    if not (direct-negative_double_commutator).is_zero():
        raise ArithmeticError('number-pointer energy identity failed')
    return direct, energy.scale


def execute(geometry, initial, duration=Fraction(1, 8), steps=4, *, mode=0):
    """Execute all source branches on a normalized state/reference matrix.

    Initial columns describe the prepared field and any untouched reference;
    preparation is supplied to this protocol. The source program is executed,
    not replaced by a matrix exponential. The latter is an independent target.
    """
    if not isinstance(geometry, field.FieldGeometry) or geometry.volume > 2:
        raise ValueError('one or two complete sites in the explicit source domain required')
    if type(duration) not in (int, Fraction):
        raise ValueError('exact rational Hamiltonian model duration required')
    if type(steps) is not int or not 1 <= steps <= 128:
        raise ValueError('one through 128 split steps required')
    space = field.FullFock(geometry.modes)
    columns = space.state(initial)
    if columns.ndim == 1:
        columns = columns[:, None]
    if columns.shape[1] > 16 or abs(np.vdot(columns, columns).real-1) > 2e-10:
        raise ValueError('normalized input with at most sixteen untouched-reference columns required')
    work, scale = pointer_work(geometry, mode)
    h = space.hamiltonian(geometry)
    target = space.exact_evolve(geometry, columns, float(duration))
    word = field.source_word(geometry, float(duration), steps)
    evolved = execute_complete(word, geometry.modes, columns)
    pointer_initial = np.zeros((2*space.dimension, columns.shape[1]), complex)
    pointer_initial[:space.dimension] = evolved.success
    pointer_word = compiler.cnot(mode, geometry.modes)
    read = execute_complete(pointer_word, geometry.modes+1, pointer_initial)
    outputs = tuple(read.success[j*space.dimension:(j+1)*space.dimension] for j in (0, 1))
    n = space.number(mode)
    ideal_outputs = ((sparse.eye(space.dimension)-n)@evolved.success,
                     n@evolved.success)
    branch_probabilities = [float(np.vdot(a, a).real) for a in outputs]
    abort_probability = evolved.failure_probability+read.failure_probability
    total_probability = math.fsum(branch_probabilities)+abort_probability
    e0 = expectation(columns, h)
    et = expectation(target, h)
    e1 = expectation(evolved.success, h)
    branches = [expectation(a, h) for a in outputs]
    predicted_read_change = expectation(evolved.success, work.matrix(scale))
    measured_read_change = math.fsum(branches)-e1
    evolution_error = math.sqrt(float(np.linalg.norm(evolved.success-target))**2+
                                evolved.failure_probability)
    pointer_error = math.sqrt(math.fsum(float(np.linalg.norm(a-b))**2
                                       for a, b in zip(outputs, ideal_outputs, strict=True))+
                              read.failure_probability)
    c, u = Fraction(geometry.speed), Fraction(geometry.interaction)
    split = complete_history_bound(geometry.sides, ((duration, steps),), c, u)
    b = stability_bounds(geometry.sides, c, u)['generator_norm_sum_upper']
    # Half-diamond distance is at most one, but a coherent dilation-vector
    # difference need not be. Keep its operator bound separate.
    vector_bound = min(Fraction(2), b**3*abs(duration)**3/(3*steps**2))
    diagnostics = [total_probability, e0, et, e1, *branches, evolution_error,
                   pointer_error, measured_read_change, predicted_read_change]
    if (not all(math.isfinite(v) for v in diagnostics) or
            abs(total_probability-1) > 2e-8 or abs(et-e0) > 2e-9 or
            evolution_error > float(vector_bound)+2e-8 or pointer_error > 2e-8 or
            abs(measured_read_change-predicted_read_change) > 2e-8):
        raise ArithmeticError('complete source energy-accounting diagnostic failed')
    return {'geometry': list(geometry.sides), 'full_Fock_dimension': space.dimension,
            'Hamiltonian_model_duration': str(Fraction(duration)), 'split_steps': steps,
            'initial_field_energy': e0, 'exact_target_energy': et,
            'evolution_success_field_energy': e1,
            'pointer_branch_probabilities': branch_probabilities,
            'pointer_branch_unnormalized_energies': branches,
            'predicted_pointer_energy_change': predicted_read_change,
            'executed_pointer_success_energy_change': measured_read_change,
            'exact_target_energy_drift': abs(et-e0),
            'source_success_energy_drift': abs(e1-e0),
            'ideal_split_half_diamond_upper': str(split),
            'ideal_split_operator_norm_upper': str(vector_bound),
            'ideal_unitary_energy_drift_upper': str(2*b*split),
            'source_complete_evolution_error_diagnostic': evolution_error,
            'source_complete_pointer_error_diagnostic': pointer_error,
            'retained_export_abort_probability': abort_probability,
            'complete_probability': total_probability,
            'evolution_resources': compiler.resources(word),
            'pointer_resources': compiler.resources(pointer_word),
            'initial_state_and_reference_supplied': True,
            'success_normalization_applied': False,
            'abort_output_energy_assigned': False,
            'physical_total_energy_balance_closed': False,
            'outward_floating_or_hardware_error_certified': False}, evolved, read
