"""Separate literal-rule and direct-matrix reference for the C10 clock adapter.

This numerical reference imports the supplied C10 source but no clock class,
Kraus builder, sparse Liouville generator, or counting implementation.
"""
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
C10 = ROOT / "openwave-tuning/next/C10"
sys.path.insert(0, str(C10))
from joint_source import JointNativeSource


def literal_instrument(model):
    """Reconstruct complete MIX from integer repair and interior CAR parity."""
    entries = {}
    event_columns = 0
    for column, (mask, flux) in enumerate(model.states):
        for seam, (a, b) in enumerate(((0, 1), (1, 2), (2, 0))):
            contrast = ((mask >> a) & 1) - ((mask >> b) & 1)
            for coin in (0, 1):
                current = (contrast - (2*coin-1) * (abs(contrast) % 2)) // 2
                target_mask, parity = mask, 1
                if current:
                    source, target = (a, b) if current > 0 else (b, a)
                    assert (mask >> source) & 1 and not (mask >> target) & 1
                    lo, hi = sorted((source, target))
                    parity = (-1)**sum((mask >> z) & 1 for z in range(lo+1, hi))
                    target_mask ^= (1 << source) | (1 << target)
                target_flux = list(flux)
                target_flux[2*seam] = (target_flux[2*seam] + current) % 6
                target_flux[2*seam+1] = (target_flux[2*seam+1] - current) % 6
                row = model.source.index[(target_mask, tuple(target_flux))]
                for law, value in (("C", current), ("E", contrast)):
                    key = seam, law, coin, value, current
                    rows, cols, values = entries.setdefault(key, ([], [], []))
                    rows.append(row)
                    cols.append(column)
                    values.append(parity / math.sqrt(12))
                event_columns += 1
    return [(key, sparse.csr_matrix((values, (rows, cols)), shape=(384, 384)))
            for key, (rows, cols, values) in sorted(entries.items())], event_columns


def direct_reference(model, ops, rho, duration, rate):
    """Integrate density matrices directly, without a vectorized generator."""
    blocks, packed, offset = [], [], 0
    for ids, *_ in model.spectral.values():
        state = rho[np.ix_(ids, ids)]
        if not np.any(state):
            continue
        dim = len(ids)
        local_h = sparse.csr_matrix(model.hamiltonian[np.ix_(ids, ids)])
        local_ops = [k[ids][:, ids] for _, k in ops]
        blocks.append((ids, dim, slice(offset, offset+dim*dim), local_h, local_ops))
        packed.append(state.reshape(-1))
        offset += dim*dim

    def rhs(_time, vector):
        derivative = np.zeros_like(vector)
        for _ids, dim, segment, h, local_ops in blocks:
            state = vector[segment].reshape(dim, dim)
            out = -1j * (h @ state - (h.T @ state.T).T) - rate * state
            for k in local_ops:
                out += rate * ((k @ state) @ k.T)
            derivative[segment] = out.reshape(-1)
        return derivative

    solution = solve_ivp(rhs, (0, duration), np.concatenate(packed), method="DOP853",
                         rtol=2e-11, atol=2e-13, t_eval=[duration])
    assert solution.success
    out = np.zeros_like(rho, dtype=complex)
    for ids, dim, segment, _, _ in blocks:
        out[np.ix_(ids, ids)] = solution.y[segment, -1].reshape(dim, dim)
    return out, solution.nfev, blocks


def single_event_quadrature(model, rho, duration, rate, blocks, order):
    nodes, weights = leggauss(order)
    out = np.zeros_like(rho, dtype=complex)
    for ids, dim, _, _, local_ops in blocks:
        h = model.hamiltonian[np.ix_(ids, ids)]
        state = rho[np.ix_(ids, ids)]
        block = np.zeros((dim, dim), complex)
        for node, weight in zip(nodes, weights):
            at = duration * (node + 1) / 2
            first = expm(-1j * at * h)
            second = expm(-1j * (duration-at) * h)
            before = first @ state @ first.conj().T
            committed = sum((k @ before) @ k.T for k in local_ops)
            block += weight * second @ committed @ second.conj().T
        out[np.ix_(ids, ids)] = math.exp(-rate*duration) * rate * duration/2 * block
    return out


def run():
    receipt = json.loads((HERE / "COMMON_CLOCK.json").read_text())
    archive = np.load(HERE / receipt["state_archive"])
    model = JointNativeSource()
    rho = archive["input_state"]
    duration = receipt["parameters"]["duration"]
    rate = receipt["parameters"]["native_attempt_rate"]
    ops, event_columns = literal_instrument(model)
    gram = sum(k.T @ k for _, k in ops).toarray()
    kraus_error = 0.0
    for key, k in ops:
        seam, law, coin, value, _ = key
        matching = next(native for meta, native in model.field.operators_by_law["MIX"][seam]
                        if meta["law"] == law and meta["coin"] == coin and meta["record_value"] == value)
        kraus_error = max(kraus_error, float(np.max(abs(k.toarray()-matching/math.sqrt(3)))))
    final, evaluations, blocks = direct_reference(model, ops, rho, duration, rate)
    one24 = single_event_quadrature(model, rho, duration, rate, blocks, 24)
    one32 = single_event_quadrature(model, rho, duration, rate, blocks, 32)
    errors = {"complete_literal_gram": float(np.max(abs(gram-np.eye(384)))),
              "literal_vs_supplied_K": kraus_error,
              "direct_ODE_vs_sparse_semigroup": float(np.max(abs(final-archive["final_state"]))),
              "quadrature24_vs32": float(np.max(abs(one24-one32))),
              "quadrature_vs_count1_generator": float(np.max(abs(one32-archive["count1_state"]))) }
    # Dropping zero-current Kraus operators destroys normalization immediately.
    moving = sum(k.T @ k for key, k in ops if key[-1] != 0).toarray()
    noops_lost = float(np.trace((np.eye(384)-moving) @ rho).real)
    retained = archive["count0_state"] + archive["count1_state"]
    tail_lost = float(1-np.trace(retained).real)
    conditional_mutation = retained / np.trace(retained)
    # A deterministic coherent kick at the mean wait is not the one-event law.
    uhalf = expm(-.5j*duration*model.hamiltonian)
    mid = uhalf @ rho @ uhalf.conj().T
    kick = sum((k @ mid) @ k.T for _, k in ops)
    wrong_one = rate*duration*math.exp(-rate*duration) * uhalf @ kick @ uhalf.conj().T
    midpoint_error = float(np.linalg.norm(wrong_one-one32))
    checks = {
        "all2304_literal_event_columns_and30_complete_marks": event_columns == 2304 and len(ops) == 30,
        "literal_CAR_signs_and_record_updates_match_supplied_K": kraus_error < 1e-14,
        "complete_literal_instrument": errors["complete_literal_gram"] < 1e-14,
        "independent_direct_matrix_ODE_agrees": errors["direct_ODE_vs_sparse_semigroup"] < 2e-10,
        "independent_ordered_time_integral_agrees": errors["quadrature_vs_count1_generator"] < 2e-12 and errors["quadrature24_vs32"] < 2e-12,
        "dropping_native_noops_rejected": noops_lost > .5,
        "discarding_or_renormalizing_count_complement_rejected": tail_lost > .03 and np.linalg.norm(conditional_mutation-final) > 1e-3,
        "replacing_random_wait_by_mean_wait_rejected": midpoint_error > 1e-4,
        "every_primary_claim_check_true": all(receipt["checks"].values()),
    }
    paths = [Path(__file__), HERE / "common_clock.py", HERE / "COMMON_CLOCK.json",
             HERE / receipt["state_archive"], C10 / "joint_source.py", C10 / "charged_reservoir.py",
             C10 / "modular_field_adapter.py"]
    result = {"schema": "oph.common_c10.clock.reference_review.v1",
              "checks": {k: bool(v) for k, v in checks.items()}, "errors": errors,
              "native_event_columns": event_columns, "complete_native_marks": len(ops),
              "direct_ODE_function_evaluations": evaluations,
              "mutations": {"noops_dropped_initial_probability": noops_lost,
                            "tail_dropped_probability": tail_lost,
                            "tail_renormalization_Frobenius_difference": float(np.linalg.norm(conditional_mutation-final)),
                            "mean_wait_substitution_count1_Frobenius_difference": midpoint_error},
              "independence": "Separate literal integer-current and CAR-parity Kraus reconstruction; direct density-matrix DOP853 integration and independent ordered-time Gaussian quadrature. No CommonClock implementation imported. Same author reference check, not an external auditor claim.",
              "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
              "official_rows_changed": []}
    (HERE / "REFERENCE_REVIEW.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))
    if not all(checks.values()):
        raise SystemExit("reference review failed")
    return result


if __name__ == "__main__":
    run()
