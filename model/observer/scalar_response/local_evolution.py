"""Execute local scalar updates and compose their detector with native feedback.

Position recurrences execute the declared canonical kick--drift--kick (KDK)
map using only a site and its six neighbors. Three columns carry the prepared
mean and two detector probes. The latter recover the evolved original-vacuum
variance without allocating a dense covariance matrix. Spectral transforms
evaluate the supplied initial vacuum only; they do not advance the columns.

These are floating diagnostics of a supplied free field. The independently
verified parent bounds concern the ideal mathematical map, not roundoff,
physical hardware, or a native quantum-control compilation. Checkpoints are
summaries, not a complete authenticated per-site event history.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from numbers import Integral
from pathlib import Path

import numpy as np

from . import experiment
from .field import ScalarField
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))
    from tools.pin_scope import bound_pins

TEMPORAL_PARENT = "code/source_scalar_time_refinement/source_scalar_time_refinement_receipt.json"
TEMPORAL_SHA = "a7b8c0a8581968494bf67cde2fd854f624829b9fd58e074ed276d3bf5d211467"
LOCAL_GRIDS = (13, 34)
REFERENCE_TIME = Q(39, 40)


def temporal_source():
    experiment.sources()
    path = experiment.ROOT / TEMPORAL_PARENT
    if hashlib.sha256(path.read_bytes()).hexdigest() != TEMPORAL_SHA:
        raise ValueError("Temporal parent receipt changed")
    parent = json.loads(path.read_text())
    for name, digest in bound_pins(parent["source_pins"]).items():
        if hashlib.sha256((experiment.ROOT / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Temporal parent source changed: {name}")
    return parent


def transfer_bound(level, reference_time=REFERENCE_TIME):
    """Transfer the parent's bound to this exact inside-window step readout."""
    if not isinstance(reference_time, Q) or not Q(19, 20) <= reference_time <= 1:
        raise ValueError("reference time must be an exact fraction in [19/20,1]")
    step = Q(level["time_step"])
    nearest = (reference_time / step + Q(1, 2)).__floor__()
    steps = min(level["last_window_step"], max(level["first_window_step"], nearest))
    model_time = steps * step
    temporal = Q(level["temporal_probability_error_upper"])
    same_time = Q(level["same_grid_total_error_upper"])
    # The continuum probability has Lipschitz constant at most 2 on this
    # apparatus. Use the actual exact rounding distance, not a fitted clock.
    reference = min(Q(1), same_time + 2 * abs(model_time-reference_time))
    return {"q": level["q"], "step": str(step), "steps": steps,
            "model_time": str(model_time), "reference_time": str(reference_time),
            "reference_time_distance": str(abs(model_time-reference_time)),
            "temporal_probability_error_upper": str(temporal),
            "same_time_probability_error_upper": str(same_time),
            "reference_probability_error_upper": str(reference),
            "same_time_history_tv_upper": str(5*same_time/33),
            "same_time_current_error_upper": str(5*same_time/132),
            "reference_history_tv_upper": str(5*reference/33),
            "reference_current_error_upper": str(5*reference/132),
            "split_native_current_magnitude_lower": str(5*Q(level["split_response_lower"])/132)}


class LocalEvolution:
    """Three finite grid columns advanced by nearest-neighbor recurrence.

The positive scalar action and mass-whitened phase space are supplied.
Independent compatibility is equality of shared field coordinates; fixed
Dirichlet boundary values are zero. The schedule reads each complete old
layer and commits the new layer simultaneously. This is conservative field
propagation, not the dissipative native integer averaging rule.
"""

    def __init__(self, q):
        if isinstance(q, (bool, np.bool_)) or not isinstance(q, Integral) or not 2 <= q <= 34:
            raise ValueError("local execution requires an integer grid from 2 through 34")
        self.field = ScalarField(q)
        self.shape = (int(q)-1,)*3
        self.diagonal = self.field.stiffness_diagonal / self.field.weights
        self.off = self.field.stiffness_off_diagonal / np.sqrt(
            self.field.weights[:-1]*self.field.weights[1:])
        mass_root = np.sqrt(self.field.weights)
        transverse = mass_root*self.field.transverse_profile
        product = np.outer(transverse, transverse)
        self.preparation = (mass_root*self.field.preparation_profile)[:, None, None]*product
        self.detector = (mass_root*self.field.detector_profile)[:, None, None]*product
        self.operator_bound = 1 + 12/min(self.field.gaps)**2

    def apply_action(self, values):
        """A=I+Lx+Ly+Lz, assembled as a seven-point local stencil."""
        if np.ma.isMaskedArray(values):
            raise ValueError("masked field values are not supported")
        values = np.asarray(values)
        if (values.shape[-3:] != self.shape or values.dtype.kind not in "fiu"
                or not np.all(np.isfinite(values))):
            raise ValueError("real field columns on this three-dimensional grid required")
        result = values.astype(float, copy=True)  # mass-one onsite potential
        for axis in range(values.ndim-3, values.ndim):
            shape = [1]*values.ndim
            shape[axis] = len(self.diagonal)
            result += self.diagonal.reshape(shape)*values
            left, right = [slice(None)]*values.ndim, [slice(None)]*values.ndim
            left[axis], right[axis] = slice(None, -1), slice(1, None)
            shape[axis] -= 1
            coupling = self.off.reshape(shape)
            result[tuple(left)] += coupling*values[tuple(right)]
            result[tuple(right)] += coupling*values[tuple(left)]
        return result

    def _modes(self, values):
        result = values
        for axis in range(3):
            result = np.moveaxis(np.tensordot(self.field.eigenvectors.T, result,
                                            axes=(1, axis)), 0, axis)
        return result

    def evolve(self, step, steps):
        if not isinstance(step, Q) or step <= 0:
            raise ValueError("positive exact rational model step required")
        if (isinstance(steps, (bool, np.bool_)) or not isinstance(steps, Integral)
                or not 0 <= steps <= 32768
                or steps*self.field.mode_count > 400_000_000):
            raise ValueError("nonnegative integer steps within the local execution budget required")
        tau = float(step)
        if not math.isfinite(tau) or tau*tau*self.operator_bound >= 4:
            raise ValueError("local canonical step must obey the strict stability bound")
        # Columns: physical prepared mean, T_qq*g, T_qp*g. Symmetry of
        # these polynomials in A makes the last two the detector pullback.
        position = np.stack((np.zeros(self.shape), self.detector, np.zeros(self.shape)))
        momentum = np.stack((4*self.preparation, np.zeros(self.shape), self.detector))
        previous = position - tau*momentum - .5*tau*tau*self.apply_action(position)
        checkpoints = []
        selected = {1, max(1, steps//2), int(steps)}
        for tick in range(1, int(steps)+1):
            following = 2*position - previous - tau*tau*self.apply_action(position)
            previous, position = position, following
            if tick in selected:
                checkpoints.append({"step": tick, "model_time": str(tick*step),
                                    "mean": float(np.sum(self.detector*position[0])),
                                    "position_column_norms": [float(np.linalg.norm(x)) for x in position]})
        mean = float(np.sum(self.detector*position[0]))
        reciprocal_mean = float(4*np.sum(self.preparation*position[2]))
        a, b = self._modes(position[1]), self._modes(position[2])
        lambdas = self.field.laplacian_eigenvalues
        omega = np.sqrt(1 + lambdas[:, None, None] + lambdas[None, :, None]
                        + lambdas[None, None, :])
        variance = float(np.sum(a*a/(2*omega) + b*b*omega/2))
        plus = (1+math.exp(-variance/2)*math.sin(mean))/2
        if not all(math.isfinite(x) for x in (mean, variance, reciprocal_mean, plus)):
            raise ArithmeticError("unresolved local field calculation")
        return {"q": self.field.q, "step": str(step), "steps": int(steps),
                "model_time": str(step*steps), "oscillator_count": self.field.mode_count,
                "mean": mean, "reciprocal_mean": reciprocal_mean,
                "variance": variance, "plus_probability": plus,
                "response": plus-.5, "checkpoints": checkpoints,
                "executed_position_column_updates": 3*int(steps)*self.field.mode_count,
                "propagated_columns": ["prepared_mean", "detector_position", "detector_momentum"],
                "full_per_site_event_log_retained": False}


def split_oracle(field, step, steps):
    """Independent analytic modal power of the SAME local KDK gates."""
    eigenvalues = field.laplacian_eigenvalues
    transverse = np.outer(field.transverse_modes**2, field.transverse_modes**2)
    base = 1 + eigenvalues[:, None] + eigenvalues[None, :]
    means, variances = [], []
    for eigenvalue, f, g in zip(eigenvalues, field.preparation_modes, field.detector_modes, strict=True):
        omega = np.sqrt(base+eigenvalue)
        z = float(step)*omega/2
        if np.any(z >= 1):
            raise ValueError("unstable modal split step")
        s = np.sqrt(1-z*z)
        phase = 2*steps*np.arcsin(z)
        sine, cosine = np.sin(phase), np.cos(phase)
        means.append(float(np.sum(4*f*g*transverse*sine/(omega*s))))
        variances.append(float(np.sum(g*g*transverse*(cosine*cosine+sine*sine/(s*s))/(2*omega))))
    mean, variance = math.fsum(means), math.fsum(variances)
    return {"mean": mean, "variance": variance,
            "plus_probability": (1+math.exp(-variance/2)*math.sin(mean))/2}


def terminal_branch(evolution, current, outcome):
    """Two Weyl-shifted copies of the split-evolved squeezed Gaussian state."""
    if type(current) is not int or current not in (-1, 0, 1) or type(outcome) is not int or outcome not in (-1, 1):
        raise ValueError("native preparation current and binary outcome required")
    mean, variance = current*evolution["mean"], evolution["variance"]
    overlap = math.exp(-variance/2)*complex(math.cos(mean), math.sin(mean))
    probability = (1+outcome*overlap.imag)/2
    scale = 1/(2*math.sqrt(probability))
    return {"current": current, "outcome": outcome, "probability": probability,
            "q": evolution["q"], "step": evolution["step"], "steps": evolution["steps"],
            "state_program": "U_KDK(step)^steps D(momentum=4*current*f)|original_vacuum>",
            "basis": "exp(i*s*Phi(g))|split_gaussian>; s=-1/2,+1/2",
            "coefficients": [[scale, 0.], [0., -outcome*scale]],
            "gram": [[[1., 0.], [overlap.real, overlap.imag]],
                     [[overlap.real, -overlap.imag], [1., 0.]]]}


def run(trials=experiment.TRIALS, seed=experiment.SEED):
    parent = temporal_source()
    bounds = [transfer_bound(row) for row in parent["levels"]]
    executions, checks = [], {}
    for level, bound in zip(parent["levels"], bounds, strict=True):
        if level["q"] not in LOCAL_GRIDS:
            continue
        local = LocalEvolution(level["q"])
        step, steps = Q(bound["step"]), bound["steps"]
        value = local.evolve(step, steps)
        oracle = split_oracle(local.field, step, steps)
        continuous = local.field.response(float(Q(value["model_time"])))
        for key in ("mean", "variance", "plus_probability"):
            checks[f"q{level['q']}_local_{key}_matches_modal_power"] = abs(value[key]-oracle[key]) < 2e-10
        checks[f"q{level['q']}_reciprocity"] = abs(value["mean"]-value["reciprocal_mean"]) < 2e-10
        checks[f"q{level['q']}_original_vacuum_squeezing"] = value["variance"] > continuous["vacuum_variance"]
        probability_error = abs(value["plus_probability"]-continuous["plus_probability"])
        checks[f"q{level['q']}_temporal_bound"] = probability_error <= float(Q(bound["temporal_probability_error_upper"]))
        sampled = experiment.sample(value["response"], trials=trials, seed=seed)
        for name, statistic in sampled["statistics"].items():
            checks[f"q{level['q']}_{name}_sample"] = (
                abs(statistic["mean"]-statistic["forecast"]) <= 6*statistic["standard_error"]+1/trials)
        branches = [terminal_branch(value, j, y) for j in (-1, 0, 1) for y in (-1, 1)]
        executions.append({"local": value, "modal_oracle": oracle,
                           "continuous_same_time": continuous,
                           "absolute_probability_error": probability_error,
                           "terminal_field_branches": branches, "sample": sampled})
        if level["q"] == 13:
            finer = local.evolve(step/2, 2*steps)
            fine_error = abs(finer["plus_probability"]-continuous["plus_probability"])
            # A measured refinement check, not the source of the proof.
            checks["q13_half_step_improves_same_time_response"] = fine_error < probability_error/3
            executions[-1]["half_step_refinement"] = {
                "local": finer, "absolute_probability_error": fine_error,
                "temporal_probability_error_upper": str(Q(bound["temporal_probability_error_upper"])/4)}
    result = {"schema": "oph.observer.scalar_local_evolution.v1",
              "source_pins": bound_pins({TEMPORAL_PARENT: TEMPORAL_SHA, **experiment.PINS}),
              "reference_time": str(REFERENCE_TIME), "locally_executed_grids": list(LOCAL_GRIDS),
              "scope": {
                  "evolution": "three position columns, local seven-point KDK recurrence; all sites and steps executed",
                  "vacuum": "supplied original action vacuum, including split squeezing; spectral initial covariance evaluation",
                  "record": "native sampled histories and field checkpoints; no full authenticated per-site field log",
                  "higher_grids": "q89 and q233 have transferred analytic bounds only; no local trajectory execution here",
                  "time": "supplied dimensionless scalar model coordinate c*t/L; update count is an algorithmic parameter",
                  "numerics": "binary64 diagnostics checked independently; exact bounds exclude implementation roundoff",
                  "physical_status": "conditional free scalar apparatus; no empirical OPH validation or common microscopic law"},
              "refinement_bounds": bounds, "executions": executions,
              "checks": checks, "all_pass": all(checks.values())}
    if not result["all_pass"]:
        raise ValueError("Local evolution check failed: " + ", ".join(k for k, v in checks.items() if not v))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=experiment.TRIALS)
    parser.add_argument("--seed", type=int, default=experiment.SEED)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path, help="create a result without overwriting")
    group.add_argument("--check", type=Path, help="replay exact fields and floating diagnostics within 2e-10")
    args = parser.parse_args()
    result = run(args.trials, args.seed)
    if args.write:
        with args.write.open("x") as stream:
            stream.write(experiment.canonical(result))
    if args.check:
        experiment.compare_receipt(json.loads(args.check.read_text()), result)
    print(json.dumps({"all_pass": result["all_pass"],
                      "locally_executed_grids": result["locally_executed_grids"],
                      "runs": [{"q": x["local"]["q"], "steps": x["local"]["steps"],
                                "model_time": x["local"]["model_time"],
                                "predicted_current": x["sample"]["statistics"]["current"]["forecast"],
                                "probability_error": x["absolute_probability_error"]}
                               for x in result["executions"]]}, indent=2))


if __name__ == "__main__":
    main()
