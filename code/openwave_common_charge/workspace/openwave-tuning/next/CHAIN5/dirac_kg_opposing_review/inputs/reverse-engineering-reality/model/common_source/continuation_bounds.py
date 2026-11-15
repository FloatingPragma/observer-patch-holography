"""Directed occupation and temporal bounds for finite smooth continuations.

The certificate follows one retained history through every Hamiltonian segment,
smooth Ramsey read, gauge read and committed-record gauge feedback.  It uses
the declared vacuum preparation and coupling 1/4.  Sharp reads, postselection,
resets, arbitrary inputs, floating execution and physical services are outside
this theorem.  See COMPOSITION.md for the stagewise proof and owned premises.
"""
from fractions import Fraction as F

from ..observer.source_scalar.bounds import root_upper
from ..observer.source_scalar.cutoff import exponential_upper
from . import refinement as r


KINDS = frozenset(("evolve", "gauge_read", "smooth_read", "feedback"))


def _rational(value, description):
    if type(value) not in (int, F) or value < 0:
        raise ValueError(f"nonnegative exact rational {description} required")
    return F(value)


def _coupling(coupling):
    if type(coupling) not in (int, F) or coupling != r.COUPLING:
        raise ValueError("this continuation certificate requires coupling exactly 1/4")


def _round_upper(value, scale):
    """Round a nonnegative rational upward to the shared certificate grid."""
    return F(-(-value.numerator * scale // value.denominator), scale)


def validated_stages(operations):
    """Normalize a finite operation list without depending on the runtime.

    Runtime Operation objects use kind/time/record/enabled attributes.  The
    standalone pair format is ("evolve", rational_time), ("smooth_read", 1),
    ("gauge_read", 0), ("feedback", earlier_record_index).  Disabled feedback
    still requires an existing record.  No branch is conditioned or removed.
    """
    if not isinstance(operations, (tuple, list)):
        raise ValueError("a finite tuple or list of protocol operations required")
    stages, records = [], 0
    for operation in operations:
        if isinstance(operation, (tuple, list)):
            if len(operation) != 2:
                raise ValueError("each certificate stage must be a kind/value pair")
            kind, value = operation
        else:
            kind = getattr(operation, "kind", None)
            time = getattr(operation, "time", None)
            record = getattr(operation, "record", None)
            enabled = getattr(operation, "enabled", None)
            if type(enabled) is not bool or (kind != "feedback" and not enabled):
                raise ValueError("only feedback admits a boolean enable switch")
            time = _rational(time, "operation time")
            if kind != "evolve" and time != 0:
                raise ValueError("only evolve operations carry Hamiltonian model time")
            if kind != "feedback" and record is not None:
                raise ValueError("only feedback operations reference an earlier record")
            value = (time if kind == "evolve" else record if kind == "feedback"
                     else 1 if kind == "smooth_read" else 0)
        if not isinstance(kind, str) or kind not in KINDS:
            raise ValueError("only evolve, gauge_read, smooth_read and feedback are certified; sharp reads are excluded")
        if kind == "evolve":
            value = _rational(value, "Hamiltonian model time")
        elif kind == "feedback":
            if type(value) is not int or not 0 <= value < records:
                raise ValueError("feedback requires the index of an already committed record")
        else:
            expected = 1 if kind == "smooth_read" else 0
            if type(value) not in (int, F) or value != expected:
                raise ValueError(f"{kind} requires certificate parameter {expected}")
            value = F(value)
            records += 1
        stages.append((kind, value))
    return tuple(stages)


def occupation_bound(d, operations, *, coupling=r.COUPLING):
    """Bound the complete ideal recorded output, uniformly in gauge feedback.

    Every even d >= 2 defines an analytic Galerkin comparison.  This does not
    imply source compilation at those dimensions.  The recurrence stores a
    projected vector error; it adds the outside tail only at the final stage.
    """
    d = r._dimension(d)
    _coupling(coupling)
    stages = validated_stages(operations)
    scale = 10**30 * (d + 2)**6
    amplitude = root_upper(r.initial_fourth_moment(), scale=scale)
    omitted_exact = r.preparation_tail_upper(d)
    # Keep the receipt polynomial in size even when the exact exponential
    # preparation tail has thousands of denominator digits.  Its upward
    # rounding contributes O(d^-6) to the vector bound, with no fixed floor.
    tail_scale = scale**2
    omitted = _round_upper(omitted_exact, tail_scale)
    preparation = root_upper(2 * omitted, scale=scale)
    projected = preparation
    inverse_root = root_upper(F(1, d + 2), scale=scale)
    h_leakage = (r.KAPPA + abs(r.COUPLING) * inverse_root) / (d + 2)
    pointer_leakage = inverse_root / (2 * (d + 2))
    rows, total_time, pointer_count, records = [], F(0), 0, 0
    for index, (kind, value) in enumerate(stages):
        incoming = amplitude
        growth, integral, residual = F(1), F(0), F(0)
        if kind == "evolve":
            total_time += value
            growth = exponential_upper(r.MOMENT_RATE * value / 2)
            integral = incoming * (growth - 1) / (r.MOMENT_RATE / 2)
            residual = integral * h_leakage
        elif kind == "smooth_read":
            pointer_count += 1
            records += 1
            growth = exponential_upper(F(1))
            integral = incoming * (growth - 1)
            residual = integral * pointer_leakage
        elif kind == "gauge_read":
            records += 1
        # Bound denominator growth without reducing either upper bound.
        # Each rounding adds less than 1/scale.  For a fixed finite program,
        # propagated amplitude rounding and accumulated residual rounding
        # remain O(d^-6).  Rounding the increment preserves the exact
        # preparation + sum(reported increments) projected-error ledger.
        amplitude = _round_upper(amplitude * growth, scale)
        residual = _round_upper(residual, scale)
        projected += residual
        rows.append({
            "index": index, "kind": kind, "parameter_exact": str(value),
            "incoming_sqrt_fourth_moment_upper": str(incoming),
            "sqrt_moment_growth_factor_upper": str(growth),
            "integrated_sqrt_fourth_moment_upper": str(integral),
            "projected_duhamel_increment_upper": str(residual),
            "outgoing_sqrt_fourth_moment_upper": str(amplitude),
            "projected_history_vector_error_upper": str(projected),
            "committed_binary_records": records,
        })
    outside = amplitude / (d + 2)**2
    vector = root_upper(projected**2 + outside**2, scale=scale)
    distance = min(F(1), vector)
    return {
        "schema": "oph.common_source.smooth_continuation_occupation.v1",
        "dimension": d, "coupling_exact": str(r.COUPLING),
        "total_hamiltonian_model_time": str(total_time),
        "smooth_pointer_count": pointer_count,
        "total_pointer_auxiliary_parameter": str(pointer_count),
        "committed_binary_records": records, "stages": rows,
        "initial_fourth_number_moment_exact": str(r.initial_fourth_moment()),
        "preparation_omitted_probability_upper": str(omitted),
        "preparation_tail_rounded_outward": True,
        "preparation_vector_error_upper": str(preparation),
        "final_fourth_number_moment_upper": str(amplitude**2),
        "final_projected_history_vector_error_upper": str(projected),
        "final_outside_vector_norm_upper": str(outside),
        "complete_stinespring_vector_error_upper": str(vector),
        "complete_instrument_trace_distance_upper": str(distance),
        "complete_classical_history_tv_error_upper": str(distance),
        "each_bounded_effect_probability_error_upper": str(distance),
        "trace_distance_convention": "one half of the trace norm",
        "all_unnormalized_branches_and_quantum_outputs_retained": True,
        "same_history_evolution_after_smooth_reads_included": True,
        "bound_is_an_arbitrary_input_diamond_norm_claim": False,
        "sharp_detector_or_state_dependent_modular_flow_included": False,
        "source_splitting_floating_or_physical_service_error_included": False,
        "spatial_continuum_or_physical_clock_calibration_claimed": False,
        "scope": "declared vacuum family, fixed two-site Hamiltonian and fixed finite smooth protocol",
    }


def temporal_bound(d, time, steps):
    """Ideal palindrome versus continuous finite H_d, on every input.

    The same existing generator norm estimate extends analytically to every
    even cutoff.  Source execution itself is available only at d=2,4.
    """
    d, time = r._dimension(d), _rational(time, "Hamiltonian model time")
    if type(steps) is not int or steps < 1:
        raise ValueError("a positive integer split count required")
    if d in (2, 4):
        from . import scalar_coupling
        return {"steps": steps, **scalar_coupling.time_error_bound(
            d, steps, time=time, coupling=r.COUPLING)}
    norm_sum = (3 + (1 + 2 * r.KAPPA) * (2 * d - 1)
                + abs(r.COUPLING) * root_upper(F(4 * d - 2), scale=10**30))
    error = norm_sum**3 * time**3 / (3 * steps**2)
    return {
        "model_time_exact": str(time), "steps": steps,
        "sum_generator_norms_upper": str(norm_sum),
        "unclipped_operator_error_upper": str(error),
        "ideal_channel_half_diamond_error_upper": str(min(F(1), error)),
        "method": "palindromic degree-two agreement and unitary telescoping",
    }


def complete_error_budget(d, operations, steps=64, *, coupling=r.COUPLING):
    """Add every same-history temporal error to the full occupation error.

    This is an ideal product-formula budget.  Numerical defects and physical
    service errors need their own certificates and are not assigned zero here.
    """
    stages = validated_stages(operations)
    occupation = occupation_bound(d, stages, coupling=coupling)
    # Validate even if a program contains no Hamiltonian segment.
    temporal_bound(d, F(0), steps)
    segments = [{"operation_index": index, **temporal_bound(d, value, steps)}
                for index, (kind, value) in enumerate(stages) if kind == "evolve"]
    temporal = sum((F(item["ideal_channel_half_diamond_error_upper"])
                    for item in segments), F(0))
    cutoff = F(occupation["complete_instrument_trace_distance_upper"])
    total = min(F(1), cutoff + temporal)
    return {
        "schema": "oph.common_source.smooth_continuation_budget.v1",
        "dimension": d, "steps_per_segment": steps,
        "occupation_certificate": occupation,
        "ideal_temporal_segment_certificates": segments,
        "ideal_temporal_instrument_error_upper": str(min(F(1), temporal)),
        "complete_ideal_product_to_infinite_instrument_error_upper": str(total),
        "complete_classical_history_tv_error_upper": str(total),
        "nonvacuous_combined_bound": total < 1,
        "source_compiler_available_at_this_dimension": d in (2, 4),
        "floating_execution_or_physical_service_error_included": False,
        "scope": "same complete smooth continuation, including all later Hamiltonian segments and gauge reads",
    }


def program_bound(program):
    """Validate the runtime Program's declared preparation and return its budget."""
    if type(program.native_current) is not int or program.native_current not in (-1, 0, 1):
        raise ValueError("the declared native-current preparations are -1, 0 and 1")
    return complete_error_budget(program.dimension, program.operations, program.steps,
                                 coupling=program.coupling)


def required_service_keys(operations):
    """Obligations for preparation, every complete operation and record storage."""
    stages = validated_stages(operations)
    return ("preparation", *(f"operation_{index}:{kind}" for index, (kind, _) in enumerate(stages)),
            "storage")


def service_error_budget(ideal_error_upper, service_errors=None, *, required_stage_keys):
    """Conditional telescoping arithmetic for explicitly supplied services.

    Values are exact bounds on preparation trace distance or complete-stage
    half-diamond distance, uniform over retained records and references.  A
    label denotes a whole stage (or an already summed repeated service).
    Supply required_service_keys(operations), extended with any surrounding
    native stages.  Completeness follows from exact key coverage, never a
    caller boolean.  This function does not establish the numerical premises
    or infer missing values.  Incomplete coverage has no complete bound.
    """
    ideal = _rational(ideal_error_upper, "ideal error upper bound")
    if (not isinstance(required_stage_keys, (tuple, list))
            or not all(isinstance(key, str) and key.strip() for key in required_stage_keys)
            or len(set(required_stage_keys)) != len(required_stage_keys)
            or not {"preparation", "storage"}.issubset(required_stage_keys)):
        raise ValueError("distinct explicit required stage keys including preparation and storage required")
    required = tuple(required_stage_keys)
    if service_errors is None:
        supplied = {}
    elif not isinstance(service_errors, dict):
        raise ValueError("a mapping of labelled exact service errors required")
    else:
        supplied = {}
        for label, value in service_errors.items():
            if not isinstance(label, str) or not label.strip():
                raise ValueError("nonempty service labels required")
            if label not in required:
                raise ValueError(f"service bound for unknown stage: {label}")
            supplied[label] = _rational(value, "service error upper bound")
    total = sum(supplied.values(), F(0))
    missing = [key for key in required if key not in supplied]
    complete = not missing
    return {
        "ideal_instrument_error_upper": str(ideal),
        "supplied_service_error_bounds": {label: str(value) for label, value in supplied.items()},
        "sum_of_supplied_service_error_bounds": str(total),
        "service_bounds_supplied": service_errors is not None,
        "required_service_stage_keys": list(required),
        "missing_service_stage_keys": missing,
        "complete_service_key_coverage": complete,
        "conditional_complete_error_upper": str(min(F(1), ideal + total)) if complete else None,
        "service_bounds_independently_established_by_this_function": False,
        "scope": "conditional arithmetic; preparation and all implemented stages must be covered to assert a complete error bound",
    }
