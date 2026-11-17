"""Outward rational comparison of the fixed two-read scalar history.

This is a new consequence of the parent's full one-particle residual and
compact-tail certificates. Its single-read probability bound alone would
not control the two-bit law, including the zero-current vacuum branches.
No floating field execution or sampled outcome enters this calculation.
"""
from __future__ import annotations

from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))
    from tools.pin_scope import bound_pins


ROOT = Path(__file__).resolve().parents[3]
PARENT = "code/source_scalar_packet/source_common_scalar_receipt.json"
PARENT_SHA256 = "c050889ea4b8e46b4635ed6ab7b850f2b8af24f3cfcac51e4a2b349d865ba7de"
NATIVE = "model/observer/native.py"
NATIVE_SHA256 = "dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d"


def _sources():
    raw = (ROOT / PARENT).read_bytes()
    if sha256(raw).hexdigest() != PARENT_SHA256:
        raise ValueError("two-read bound parent pin mismatch")
    parent = json.loads(raw)
    for name, pin in bound_pins(parent["source_pins"]).items():
        if sha256((ROOT / name).read_bytes()).hexdigest() != pin["sha256"]:
            raise ValueError(f"two-read bound parent source pin mismatch: {name}")
    if sha256((ROOT / NATIVE).read_bytes()).hexdigest() != NATIVE_SHA256:
        raise ValueError("two-read bound native prior pin mismatch")
    # Load only after checking the immutable parent's arithmetic source.
    name = "_oph_scalar_continuation_intervals"
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "code/source_scalar_packet/source_scalar_intervals.py")
    arithmetic = importlib.util.module_from_spec(spec)
    sys.modules[name] = arithmetic
    spec.loader.exec_module(arithmetic)
    return parent, arithmetic


def walsh_tv_upper(a, b, c):
    """TV box bound for p(y,z) differences (y*a+z*b+y*z*c)/4.

    For the actual signed differences the same expression in their absolute
    values is an equality. No probability normalization or positivity is
    inferred from arbitrary supplied Walsh coefficients.
    """
    values = (a, b, c)
    if any(type(x) not in (int, Q) or x < 0 for x in values):
        raise ValueError("nonnegative exact rational coefficient bounds required")
    return max(*values, sum(values) / Q(2)) / 2


def run():
    """Return the certificate for t1=19/20, t2=1 and gamma=4 only.

    Accept no unvalidated receipt or alternative setting. The complete
    history includes initial state and coins; subsequent native events
    form the same retained-record stochastic continuation in both models.
    """
    parent, arithmetic = _sources()
    I, sin, exp, PI = arithmetic.I, arithmetic.sin, arithmetic.exp, arithmetic.PI
    upper = lambda x: I.of(Q(x[1]))
    exact_upper = lambda x: Q(x.hi, arithmetic.SCALE)
    export_upper = lambda x: x.encode()[1]
    export_lower = lambda x: x.encode()[0]
    t1, t2, gamma = Q(19, 20), Q(1), Q(4)
    delay = t2 - t1
    continuum = parent["continuum"]
    b = upper(continuum["coherent_mode_norm"])
    tail = (upper(continuum["smooth_comparison_tail_mass"]) / 2).sqrt().upper()
    # Reconstruct the declared smooth comparison, using exact coefficients
    # and the same certified pi enclosure, without importing the producer.
    raw = list(map(Q, parent["preparation"]["sine_coefficients"]))
    norm = I.of(Q(parent["preparation"]["coefficient_squared_norm"])).sqrt()
    coefficients = [I.of(x) / norm for x in raw]
    frequencies = [(1 + PI * PI * (k * k + 2)).sqrt() for k in range(1, 10)]
    variance = sum((a * a / (2 * w) for a, w in zip(coefficients, frequencies)), I.of(0))
    covariance = sum((a * a * sin(PI / 2 + w * delay) / (2 * w)
                      for a, w in zip(coefficients, frequencies)), I.of(0))
    commutator = sum((a * a * sin(w * delay) / w
                      for a, w in zip(coefficients, frequencies)), I.of(0))
    means = [gamma * sum((a * a * (-1) ** k * sin(w * t) / w
                         for k, (a, w) in enumerate(zip(coefficients, frequencies))), I.of(0))
             for t in (t1, t2)]
    continuum_kernel_error = tail * (2 * b + tail)

    def enclose(center, error):
        return I(center.lo - error.hi, center.hi + error.hi)

    def witnesses(var_error, kernel_error, response):
        v = enclose(variance, var_error)
        v = I(max(0, v.lo), v.hi)  # variance is nonnegative independently
        cross = enclose(covariance, kernel_error)
        comm = enclose(commutator, 2 * kernel_error)
        # The joint vacuum correlation is exp(-V)*sinh(C). Use sinh C >= C
        # and the directed exponential; no monotonicity outside C>0 is used.
        # Both witnesses are signed: an outcome-product expectation and a
        # difference of probabilities. Without the sign premises below,
        # only their universal lower bound -1 is available from this proof.
        corr = I.of(-1)
        if cross.lo > 0:
            corr = exp(-I(v.hi, v.hi)) * I(cross.lo, cross.lo)
        disturbance = I.of(-1)
        # For 0<c<1, 1-cos(c/2) >= c^2/8-c^4/384. The polynomial is
        # increasing there. Multiply by the positive lower bound on -r2.
        if 0 < comm.lo <= comm.hi < arithmetic.SCALE and response.hi < 0:
            low = I(comm.lo, comm.lo)
            loss = low ** 2 / 8 - low ** 4 / 384
            disturbance = loss * I(-response.hi, -response.hi)
        return {
            "variance_interval": v.encode(),
            "real_cross_covariance_interval": cross.encode(),
            "commutator_over_i_interval": comm.encode(),
            "vacuum_bit_correlation_lower": export_lower(corr),
            "second_plus_disturbance_given_positive_current_lower": export_lower(disturbance),
            "disturbance_reference": "same prepared field at t2 with first read omitted",
        }

    cont_response = I(I.of(Q(continuum["compact_induced_response"][0])).lo,
                      I.of(Q(continuum["compact_induced_response"][1])).hi)
    continuum_witnesses = witnesses(continuum_kernel_error, continuum_kernel_error, cont_response)
    rows = []
    for level in parent["levels"]:
        db = (upper(level["full_preparation_frequency_residual_upper"]) / 2
              + upper(level["sampling_isometry_error_upper"])) / I.of(2).sqrt()
        eta = upper(level["full_coherent_frequency_residual_upper"]) / I.of(2).sqrt()
        eg = db + (upper(level["sampled_detector_tail_mass"]) / 2).sqrt()
        ef = db + (upper(level["sampled_preparation_tail_mass"]) / 2).sqrt()
        graph_variance_error = eg * (2 * b + eg)
        graph_kernel_error = b * (2 * eg + delay * eta) + eg * (eg + delay * eta)
        dv = graph_variance_error + continuum_kernel_error
        dg = graph_kernel_error + continuum_kernel_error
        da, da_cont = gamma * ef, gamma * tail
        dm = []
        for t in (t1, t2):
            d = eg + t * eta
            dm.append(2 * (b * da + gamma * b * d + da * d
                           + b * da_cont + gamma * b * tail + da_cont * tail))
        # |sin mean * exp(-V/2)| sensitivity, optionally sharpened by the
        # existing probability certificate, which covers both chosen times.
        old = 2 * exact_upper(upper(level["graph_vs_compact_continuum_error_upper"]))
        ds = [min(Q(2), old, exact_upper(dv / 2 + x)) for x in dm]
        cmax = commutator.abs() + 2 * I.of(max(exact_upper(graph_kernel_error),
                                              exact_upper(continuum_kernel_error)))
        # |cos(c_h/2)-cos(c_inf/2)| <= max(|c_h|,|c_inf|)*|dc|/4.
        dd = min(Q(2), exact_upper(cmax * dg / 2))
        second = min(Q(2), ds[1] + dd)
        dr0 = min(Q(2), exact_upper(dv + dg))
        dr1 = min(Q(2), exact_upper(dv + dg + dm[0] + dm[1]))
        tv0 = min(Q(1), dr0 / 2)
        tv1 = min(Q(1), walsh_tv_upper(ds[0], second, dr1))
        tv = Q(28, 33) * tv0 + Q(5, 33) * tv1
        response = I(I.of(Q(level["graph_compact_induced_response"][0])).lo,
                     I.of(Q(level["graph_compact_induced_response"][1])).hi)
        row = {
            "q": level["q"],
            "detector_vector_to_band_error_upper": export_upper(eg),
            "preparation_vector_to_band_error_upper": export_upper(ef),
            "detector_evolution_residual_upper": export_upper(eta),
            "graph_to_band_variance_error_upper": export_upper(graph_variance_error),
            "graph_to_band_cross_kernel_error_upper": export_upper(graph_kernel_error),
            "variance_difference_upper": export_upper(dv),
            "cross_kernel_difference_upper": export_upper(dg),
            "positive_current_mean_difference_upper": [export_upper(x) for x in dm],
            "first_sine_moment_difference_upper": export_upper(I.of(ds[0])),
            "second_sequential_sine_moment_difference_upper": export_upper(I.of(second)),
            "attenuation_difference_upper": export_upper(I.of(dd)),
            "vacuum_bit_correlation_difference_upper": export_upper(I.of(dr0)),
            "nonzero_current_bit_correlation_difference_upper": export_upper(I.of(dr1)),
            "vacuum_conditional_history_tv_upper": export_upper(I.of(tv0)),
            "nonzero_current_conditional_history_tv_upper": export_upper(I.of(tv1)),
            "complete_history_tv_upper": export_upper(I.of(tv)),
            "finite_field_witnesses": witnesses(graph_variance_error, graph_kernel_error, response),
        }
        rows.append(row)
    return {
        "schema": "oph.observer.scalar_continuation_bounds.v1",
        "parent": {"path": PARENT, "sha256": PARENT_SHA256},
        "native_prior_source": {"path": NATIVE, "sha256": NATIVE_SHA256},
        "settings": {"times": [str(t1), str(t2)], "gamma": str(gamma),
                     "detector_couplings": ["1", "1"], "field_reset_between_reads": False},
        "scope": {"directed_rational_ideal_model_bound": True,
                  "floating_execution_roundoff_certified": False,
                  "parent_single_read_error_alone_suffices": False,
                  "fixed_band_tail_bound_proves_vanishing_rate": False,
                  "comparison": "same two centered field instruments and retained native continuation"},
        "native_first_current_probabilities": {"-1": "5/66", "0": "28/33", "1": "5/66"},
        "smooth_band_reference": {"variance": variance.encode(), "real_cross_covariance": covariance.encode(),
                                  "commutator_over_i": commutator.encode(),
                                  "positive_current_means": [x.encode() for x in means]},
        "continuum_to_band_kernel_error_upper": export_upper(continuum_kernel_error),
        "continuum_witnesses": continuum_witnesses,
        "levels": rows,
    }


def verify(candidate):
    """Replay all mathematics and source/settings bindings, not just hashes."""
    def exact(left, right):
        if type(left) is not type(right):
            return False
        if isinstance(left, dict):
            return left.keys() == right.keys() and all(exact(left[key], right[key]) for key in left)
        if isinstance(left, list):
            return len(left) == len(right) and all(exact(a, b) for a, b in zip(left, right))
        return left == right

    if not exact(candidate, run()):
        raise ValueError("two-read bound replay mismatch")
    return True
