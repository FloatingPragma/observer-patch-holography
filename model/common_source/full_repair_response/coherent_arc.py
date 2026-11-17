"""A specified coherent continuation of the ideal mixed-cap native law.

This generalizes observer/coherent_continuation.py with the existing ideal
mixed-cap seam weights. It is a different operational contract from the
finite selector with committed terminal timeouts. Its coherent event register
is never measured or reset during the walk. No timeout is restarted.

Native branch reversal is evaluated as an exact source-defined permutation;
the coin reflection is also executed on all64 columns with the existing pulse
compiler. This is not a compilation of the whole reversal into pulse gates.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np
import sympy as s

from ...native_modular_reference import common_caps
from ...observer import native
from ...observer.source_scalar import compiler
from ..mass_response import history
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))
    from tools.pin_scope import bound_pins


ROOT = Path(__file__).resolve().parents[3]


@lru_cache(maxsize=1)
def construction():
    core, weights, certificate = common_caps.cap_schedule()
    states = native.states(binary=True)
    indices = {state: i for i, state in enumerate(states)}
    events = tuple((seam, coin) for seam in sorted(weights) for coin in (0, 1))
    event_index = {event: j for j, event in enumerate(events)}
    probabilities = tuple(weights[seam]/2 for seam, _ in events)
    if len(events) != 40 or sum(probabilities, F()) != 1:
        raise ArithmeticError("the existing ideal mixed-cap event law changed")
    width = 64  # all padding amplitudes retained in the source coin execution
    reversal = list(range(len(states)*width))
    transition = Counter()
    for i, state in enumerate(states):
        for j, ((seam, coin), probability) in enumerate(zip(events, probabilities, strict=True)):
            output, (_, current) = native.native(state, seam, coin)
            other = indices[output]
            reverse_coin = coin if current == 0 else 1-coin
            back = event_index[seam, reverse_coin]
            reversal[i*width+j] = other*width+back
            if probabilities[back] != probability:
                raise ArithmeticError("branch reversal changed the coherent event weight")
            transition[other, i] += s.Rational(probability)
    if any(reversal[reversal[i]] != i for i in range(len(reversal))):
        raise ArithmeticError("native record reversal is not an involution")
    t = s.ImmutableSparseMatrix(66, 66, dict(transition))
    if t != t.T or t*s.ones(66, 1) != s.ones(66, 1):
        raise ArithmeticError("native coherent discriminant is not symmetric stochastic")
    gap = F(certificate["gap_lower_per_attempt"])
    if gap != F(1, 600):
        raise ArithmeticError("the inherited ideal completion gap bound changed")
    # Each local population map is (I+Swap_e)/2, an orthogonal projector.
    # Therefore 0<=T<=I. The source graph theorem supplies its nonzero gap.
    return {"core": core, "events": events, "probabilities": probabilities,
            "reversal": tuple(reversal), "transition": t, "width": width,
            "gap_certificate": certificate}


def ideal_coin():
    data = construction()
    vector = np.zeros(data["width"], complex)
    vector[:len(data["events"])] = np.sqrt(list(map(float, data["probabilities"])))
    return vector


@lru_cache(maxsize=1)
def source_coin():
    vector = ideal_coin()
    preparation = history.preparation_word(vector)
    # U (2|0><0|-I) U†; words execute from left to right on the state.
    reflection = (compiler.inverse(preparation)
                  + compiler.compile_diagonal([0.]+[math.pi]*63)
                  + preparation)
    actual_preparation = compiler.execute(preparation, 6, np.eye(64, dtype=complex))
    actual_reflection = compiler.execute(reflection, 6, np.eye(64, dtype=complex))
    for matrix in (actual_preparation, actual_reflection):
        matrix.flags.writeable = False
    return actual_preparation, actual_reflection, {
        "all_coin_input_columns_executed": 64,
        "preparation_pulses": len(preparation), "reflection_pulses": len(reflection),
        "complete_native_record_reversal_pulse_compilation": False,
        "native_record_reversal_exact_permutation_executed": True,
    }


def embed(vectors, *, source=False):
    vectors = np.asarray(vectors, complex)
    if vectors.ndim == 1:
        vectors = vectors[:, None]
    if vectors.ndim != 2 or vectors.shape[0] != 66 or not np.isfinite(vectors).all():
        raise ValueError("finite66-row input vectors, with optional reference columns, required")
    coin = source_coin()[0][:, 0] if source else ideal_coin()
    return (vectors[:, None, :]*coin[None, :, None]).reshape(66*64, -1)


def compress(vectors):
    vectors = np.asarray(vectors, complex)
    if vectors.ndim != 2 or vectors.shape[0] != 66*64:
        raise ValueError("complete4224-row arc input required")
    return np.einsum("e,ier->ir", ideal_coin().conj(), vectors.reshape(66, 64, -1))


def advance(vectors, *, phase=-1, source=False):
    """One W=S[Pi+phase(I-Pi)] step, retaining event and padding workspace."""
    if type(phase) is not int or phase not in (-1, 1):
        raise ValueError("one of the existing continuation phases±1 required")
    vectors = np.asarray(vectors, complex)
    if vectors.ndim != 2 or vectors.shape[0] != 66*64 or not np.isfinite(vectors).all():
        raise ValueError("complete finite4224-row arc input required")
    if phase == 1:
        reflected = vectors
    elif source:
        matrix = source_coin()[1]
        reflected = np.einsum("ef,ifr->ier", matrix, vectors.reshape(66, 64, -1)).reshape(vectors.shape)
    else:
        reflected = 2*embed(compress(vectors))-vectors
    return reflected[np.array(construction()["reversal"])]


def spectrum():
    """All arc eigenphase multiplicities, with exact algebraic dimension counts."""
    data = construction()
    t = data["transition"]
    # The Dirichlet form is sum_ij T_ij |f_i-f_j|^2/2: its kernel is
    # exactly the constants on connected positive-weight components.
    remaining, components = set(range(66)), []
    while remaining:
        reached, frontier = set(), {min(remaining)}
        while frontier:
            reached |= frontier
            frontier = {j for i in frontier for j in range(66) if t[i, j] > 0}-reached
        remaining -= reached
        components.append(sorted(reached))
    fixed = len(components)
    if fixed != 4:
        raise ArithmeticError("the four protected-core population components changed")
    # Exact strict diagonal dominance is a direct Gershgorin lower bound;
    # avoids unnecessary rational66x66 elimination. No spectrum is clipped.
    lower = min(2*t[i, i]-1 for i in range(66))
    if lower <= 0:
        raise ArithmeticError("strict positive transfer certificate changed")
    values = np.linalg.eigvalsh(np.array(t, float))
    nonfixed = values[:66-fixed]
    if not (np.all(nonfixed > 0) and np.all(nonfixed < 1)
            and np.max(abs(values[-fixed:]-1)) < 1e-12):
        raise ArithmeticError("numerical spectrum disagrees with exact fixed-space count")
    phases = np.arccos(nonfixed)[::-1]
    dimension = len(data["reversal"])
    trace_shift = sum(i == j for i, j in enumerate(data["reversal"]))
    shift_plus, shift_minus = (dimension+trace_shift)//2, (dimension-trace_shift)//2
    nontrivial = 66-fixed
    flat_zero, flat_pi = shift_minus-nontrivial, shift_plus-66
    if fixed+2*nontrivial+flat_zero+flat_pi != dimension or min(flat_zero, flat_pi) < 0:
        raise ArithmeticError("complementary arc eigenmodes were lost")
    # Since gap>=1/600, arccos(1-gap)>=sqrt(2gap)>1/18.
    # All principal absolute phases are <=pi<22/7, including flat pi modes.
    if not F(1, 300) > F(1, 18**2):
        raise ArithmeticError("elementary angle-gap certificate failed")
    return {
        "discriminant_dimension": 66,
        "exact_positive_transfer_lower_bound": str(lower),
        "fixed_component_sizes": sorted(map(len, components)),
        "physical_native_arcs": 66*len(data["events"]),
        "coin_padding_arcs_retained": 66*(64-len(data["events"])),
        "total_executed_arc_dimension": dimension,
        "discriminant_eigenvalues_diagnostic": values.tolist(),
        "cyclic_nonzero_positive_principal_phases_diagnostic": phases.tolist(),
        "cyclic_zero_phase_multiplicity": fixed,
        "cyclic_positive_phase_count_including_multiplicity": nontrivial,
        "cyclic_negative_phase_count_including_multiplicity": nontrivial,
        "complement_zero_phase_multiplicity": flat_zero,
        "complement_pi_phase_multiplicity": flat_pi,
        "all_arc_modes_accounted_for": True,
        "cyclic_principal_phase_spread_diagnostic": float(phases[-1]/phases[0]),
        "full_nonzero_principal_phase_spread_diagnostic": float(math.pi/phases[0]),
        "rigorous_principal_phase_spread_upper": "396/7 < 57",
        "angle_gap_argument": "T is a convex sum of (I+Swap)/2; gap>=1/600; acos(1-gap)>=sqrt(2gap)>1/18; |phase|<=pi<22/7",
        "bound_scope": "nonzero elementary principal absolute eigenphases referenced to phase zero; not arbitrary differences, composite fields or other logarithm branches",
        "physical_rest_mass_identification": False,
    }


def report(*, with_arrays=False):
    data = construction()
    # Two untouched reference columns retain the complete initial entanglement.
    initial = np.zeros((66, 2), complex)
    initial[0, 0] = 1/np.sqrt(2)
    initial[-1, 1] = 1j/np.sqrt(2)
    actual = embed(initial, source=True)
    expected = embed(initial)
    anchor = actual.copy()
    t = np.array(data["transition"], float)
    previous, recurrence = None, initial.copy()
    rows, trajectory = [], []
    for step in range(9):
        trajectory.append(actual.copy())
        plus, minus = (anchor+actual)/2, (anchor-actual)/2
        # Both Ramsey work outputs, including all padding and reference
        # columns, determine the probabilities; no branch is normalized.
        pplus, pminus = float(np.vdot(plus, plus).real), float(np.vdot(minus, minus).real)
        rows.append({"step": step, "source_norm": float(np.vdot(actual, actual).real),
                     "source_vs_ideal_state_residual": float(np.linalg.norm(actual-expected)),
                     "ideal_compression_recurrence_residual": float(np.max(abs(compress(expected)-recurrence))),
                     "Ramsey_complete_outcome_probabilities": [pplus, pminus],
                     "retained_padding_probability": float(np.sum(abs(actual.reshape(66, 64, 2)[:, 40:])**2))})
        if step == 8:
            break
        next_recurrence = t@recurrence if previous is None else 2*t@recurrence-previous
        previous, recurrence = recurrence, next_recurrence
        actual, expected = advance(actual, source=True), advance(expected)
    prepared = embed(initial)
    first = [advance(prepared, phase=phase) for phase in (-1, 1)]
    second = [advance(first[i], phase=phase) for i, phase in enumerate((-1, 1))]
    source_prep, source_reflection, resources = source_coin()
    target_reflection = 2*np.outer(ideal_coin(), ideal_coin().conj())-np.eye(64)
    paths = [ROOT/name for name in common_caps.SOURCE_PATHS]
    paths += [Path(__file__), ROOT/"model/observer/coherent_continuation.py",
              ROOT/"model/observer/native.py", ROOT/"model/common_source/mass_response/history.py",
              ROOT/"model/observer/source_scalar/compiler.py"]
    result = {"schema": "oph.full_repair_response.coherent_arc.v1",
            "construction": "weighted native branch reversal S with W=S(2AA†−I), extending the existing coherent continuation",
            "event_weights_exact": [str(x) for x in data["probabilities"]],
            "selector_contract": "ideal coherent preparation of the existing mixed-cap law; no committed rejection/timeout protocol is continued",
            "time_parameter": "coherent walk steps; no physical or modular-clock calibration derived",
            "source_execution": resources,
            "floating_residuals_are_physical_error_certificates": False,
            "source_coin_preparation_residual": float(np.max(abs(source_prep[:, 0]-ideal_coin()))),
            "source_coin_reflection_residual": float(np.max(abs(source_reflection-target_reflection))),
            "source_reflection_unitarity_residual": float(np.max(abs(source_reflection.conj().T@source_reflection-np.eye(64)))),
            "continuing_arc_diagnostics": rows,
            "readout_contract": "one terminal declared Ramsey instrument at a chosen horizon; earlier rows are alternative endpoints, not repeated committed measurements",
            "controlled_walk_Ramsey_pulse_word_compiled": False,
            "reference_columns_retained": 2,
            "source_workspace_reset_between_steps": False,
            "same_first_event_different_continuation": {
                "first_step_state_residual": float(np.linalg.norm(first[0]-first[1])),
                "second_step_state_distance": float(np.linalg.norm(second[0]-second[1])),
                "minus_one_reflection_phase_is_declared": True},
            "spectrum": spectrum(),
            "interpretation": {"coherent_oscillatory_continuation_exists": True,
                               "continuation_selected_by_committed_native_records": False,
                               "same_C_E_committed_history_as_repair": False,
                               "Standard_Model_charged_field_identification": False,
                               "Lorentzian_3plus1_dispersion_derived": False,
                               "observed_lepton_hierarchy_as_elementary_principal_phases": False,
                               "empirical_inputs_used_to_choose_weights_or_phases": False},
            "source_bindings": bound_pins({path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                                for path in sorted(set(paths))})}
    if with_arrays:
        final = trajectory[-1]
        return result, {"source_arc_states": np.stack(trajectory),
                        "terminal_plus": (anchor+final)/2,
                        "terminal_minus": (anchor-final)/2}
    return result
