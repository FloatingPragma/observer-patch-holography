"""Prepare a regional modular reference using incident ticks, not global time.

Fix a cap and IID uniform native seams, with fair placement coins and one
fixed law C or E. Write I and R for the unread incident and remote channels,
p for the incident seam fraction, and G=p I+(1-p)R. All three channels are
positive self-adjoint Hilbert--Schmidt contractions. On traceless operators
the source certificate gives G <= r Id, r=19079/20000. The channel at the
next incident tick (remote attempts precede its repair) is

    S = p I [Id-(1-p)R]^-1.

For M=Id-(1-p)R, p Id<=M<=Id. Thus B=M^-1/2 p I M^-1/2 is positive and
B<=Id-(1-r)M^-1<=r Id on the traceless space. Since S=M^1/2 B M^-1/2,
||S^n-Dep||_(2->2)<=p^-1/2 r^n. Purification with source Schmidt rank<=66
gives arbitrary-reference half trace distance <=33 p^-1/2 r^n. The code
uses the rational outward bound 33/p times 2^-floor(n/16), since r^16<1/2.
This proves a local-tick preparation bound even though S need not be
self-adjoint. The schedule is exogenous; adaptation to the tape is excluded.

After a fixed local warm-up, retain BOTH outcomes of the next incident tick.
Conditioning that specified tick on movement has uniform-input probability
q=5/33 and reference proportional to its incident seam degree. This is not
waiting for the first successful movement. Native complement repair selects
the conditional expectation of that reference. Its cap modular frequencies
are read from the resulting weights. No seconds, physical boost, passive
clock evolution, or physical field identification is supplied by this proof.

Records may be retained externally. Bounds concern the unread source
marginal (and an initially untouched reference), not its independence from
the generated tape. A region must collect its incident ticks, stop at the
declared count, and switch to the complement-only protocol.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, eye

from ..native_equilibration import cut_completion as cut, construction as equilibrium
from ..observer import native, quantum_mixing as mixing
from ..record_laws import burn_in
from ..record_laws.burn_in import R_UPPER

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "regional_preparation_results.json"
CAP = cut.CAP
SOURCE_PATHS = tuple(sorted(set((
    *cut.SOURCE_PATHS, *equilibrium.SOURCE_PATHS, *burn_in.SOURCE_PATHS,
    "model/native_clock/regional_preparation.py",
    "model/native_equilibration/cut_completion.py",
    "model/native_equilibration/construction.py",
    "model/observer/native.py", "model/observer/quantum_mixing.py",
    "model/observer/QUANTUM_MIXING.md", "model/record_laws/burn_in.py",
    "model/common_source/slow_geometry.py",
    "Lean/Thermodynamics/FiniteConditionalRepair.lean",
    "Lean/Thermodynamics/CapFirstLaw.lean",
))))


def _law(law):
    if type(law) is not str or law not in ("C", "E"):
        raise ValueError("one fixed native law C or E is required")
    return law


def _integer(value, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f"an integer in [{low}, {high}] is required")
    return value


def incident_edges(cap=CAP):
    cap = cut._cap(cap)
    return tuple(e for e in mixing.seams() if set(e).intersection(cap))


def local_bound(ticks, cap=CAP):
    """All-input/reference distance bound; ticks include incident zeros."""
    _integer(ticks, 0, 1000000)
    p = F(len(incident_edges(cap)), 30)
    if not 0 < R_UPPER < 1 or R_UPPER**16 >= F(1, 2):
        raise ArithmeticError("the inherited mixing radius lost its block bound")
    return min(F(1), F(33, 1)/p/2**(ticks//16))


def reference(cap=CAP):
    """Incident-movement reference and its native complement completion."""
    edges = incident_edges(cap)
    basis = native.states(binary=True)
    degrees = tuple(sum(x[a] != x[b] for a, b in edges) for x in basis)
    total = sum(degrees)
    raw = tuple(F(d, total) for d in degrees)
    q = F(total, 2*len(edges)*66)
    completed = [F(0)]*66
    for sector in cut.cut_sectors(cap):
        for row in sector.indices:
            value = sum((raw[i] for i in row), F(0))/len(row)
            for i in row:
                completed[i] = value
    if q != F(5, 33) or sum(completed) != 1 or min(completed) <= 0:
        raise ArithmeticError("local herald completion is not faithful and normalized")
    return raw, tuple(completed), q


def preparation(precision_bits=24, cap=CAP):
    """Finite local-count and complement-attempt budgets, all failures kept.

    If the local warm-up distance is delta, conditioning the next specified
    movement outcome costs at most 2 delta/(q-delta). Complement repair adds
    eta by contractivity. This bound includes initially untouched references.
    For the source marginal, |rho-tau|_op<=delta implies a multiplicative
    radius t=66 delta. A finite complement error eta changes the completed
    faithful reference by radius u=eta/min(sigma). Their product yields
    v=t+u+tu. Normalization places the prepared density between
    (1-v)/(1+v) sigma and (1+v)/(1-v) sigma. Operator monotonicity of log
    bounds the modular generator error by log((1+v)/(1-v))<=2v/(1-v).
    Duhamel gives half diamond distance of modular flows <=|s| times that
    bound. This concerns the actual prepared source marginal's modular
    flow; implementing that flow remains the coherent control protocol.
    """
    _integer(precision_bits, 12, 40)
    cc = cut.cut_certificate(cap, precision_bits=precision_bits)
    p = F(len(incident_edges(cap)), 30)
    blocks = 0
    while F(33, 1)/p/2**blocks > F(1, 2**precision_bits):
        blocks += 1
    ticks = 16*blocks
    delta = local_bound(ticks, cap)
    _, completed, q = reference(cap)
    eta = F(cc["preparation"]["trace_distance_upper"])
    t, u = 66*delta, eta/min(completed)
    v = t+u+t*u
    if delta >= q or v >= 1:
        raise ArithmeticError("finite preparation did not certify a faithful output")
    return {
        "cap": list(cap), "incident_seam_fraction": str(p),
        "local_warmup_ticks_including_zeros": ticks,
        "warmup_reference_uniform_halftrace_upper": str(delta),
        "herald": "movement flag of the NEXT specified incident tick",
        "uniform_input_success_probability": str(q),
        "success_probability_interval": [str(q-delta), str(q+delta)],
        "failure_probability_interval": [str(1-q-delta), str(1-q+delta)],
        "complement_preparation": cc["preparation"],
        "conditional_output_reference_uniform_halftrace_upper": str(min(F(1), 2*delta/(q-delta)+eta)),
        "completed_reference_minimum_eigenvalue": str(min(completed)),
        "relative_pre_normalization_radius": str(v),
        "prepared_modular_generator_operator_error_upper": str(2*v/(1-v)),
        "flow_halfdiamond_upper_per_absolute_modular_parameter": str(2*v/(1-v)),
        "expected_underlying_attempts_through_herald": str(F(ticks+1)/p),
        "remote_wait_tail_per_incident_tick": f"({1-p})^h after h native attempts",
        "global_warmup_attempt_counter_required": False,
        "complement_controller_counts_attempts_and_reports_completion": True,
    }


def seam_counts(seam, law):
    """Twice the actual unread quantum channel on all 4,356 matrix units."""
    _law(law)
    if (type(seam) is not tuple or len(seam) != 2
            or any(type(port) is not int for port in seam) or seam not in mixing.seams()):
        raise ValueError("a native seam is required")
    basis = native.states(binary=True)
    positions = {x: i for i, x in enumerate(basis)}
    rows, cols = [], []
    for coin in (0, 1):
        images = []
        for x in basis:
            out, (_, current) = native.native(x, seam, coin)
            images.append((positions[out], current if law == "C" else x[seam[0]]-x[seam[1]]))
        for i, (a, ri) in enumerate(images):
            for j, (b, rj) in enumerate(images):
                if ri == rj:
                    rows.append(66*a+b)
                    cols.append(66*i+j)
    return coo_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)),
                      shape=(4356, 4356)).tocsr()


def exact_audit(cap=CAP):
    """Check hypotheses exactly, including quantum rather than only bit laws."""
    edges = incident_edges(cap)
    identity = eye(4356, dtype=np.int64, format="csr")
    trace = np.zeros(4356, dtype=np.int64)
    trace[np.arange(66)*67] = 1
    result = {}
    for law in ("C", "E"):
        digest = hashlib.sha256()
        for seam in mixing.seams():
            a = seam_counts(seam, law)
            polynomial = a@(a-2*identity)
            if law == "C":
                polynomial = polynomial@(a-identity)
            if ((a-a.T).nnz or polynomial.nnz or
                    not np.array_equal(a@trace, 2*trace) or
                    not np.array_equal(trace@a, 2*trace)):
                raise ArithmeticError("native local positive self-adjoint channel identity failed")
            # Symmetry and this square-free polynomial place the spectrum in
            # {0,1/2,1} (C) or {0,1} (E), with no floating eigenvalue test.
            rr, cc = a.nonzero()
            digest.update(json.dumps([list(seam), [[int(i), int(j), int(a[i, j])]
                                      for i, j in zip(rr, cc, strict=True)]],
                                     separators=(",", ":")).encode())
        result[law] = {"all_seams": 30, "all_matrix_units_per_seam": 4356,
                       "exact_unital_trace_preserving_positive_self_adjoint": True,
                       "native_seam_channel_sha256": digest.hexdigest()}
    # Independently evaluate both output flags under uniform input using
    # actual primitives, including the placement probability.
    basis = native.states(binary=True)
    positions = {x: i for i, x in enumerate(basis)}
    outputs = [[F(0)]*66 for _ in range(2)]
    for x in basis:
        for seam in edges:
            for coin in (0, 1):
                out, (_, current) = native.native(x, seam, coin)
                outputs[int(current != 0)][positions[out]] += F(1, 66*2*len(edges))
    raw, completed, q = reference(cap)
    if tuple(outputs[1]) != tuple(q*w for w in raw):
        raise ArithmeticError("native local herald differs from its degree formula")
    if any(a+b != F(1, 66) for a, b in zip(*outputs, strict=True)):
        raise ArithmeticError("the complete herald changed the uniform source")
    e, scale = cut.completion_counts(cap)
    # Diagonal exact check, avoiding floats even for the rational completion.
    for sector in cut.cut_sectors(cap):
        for row in sector.indices:
            if any(completed[i] != sum((raw[j] for j in row), F(0))/len(row) for i in row):
                raise ArithmeticError("native tensor completion has the wrong reference")
    if (e@e-scale*e).nnz:
        raise ArithmeticError("native tensor completion is not idempotent")
    return {"laws": result, "both_herald_outputs_retained": True,
            "native_incident_herald_verified": True, "completion_verified": True}


def run():
    raw, completed, q = reference()
    gaps = []
    for sector in cut.cut_sectors():
        weights = sorted(set(completed[row[0]] for row in sector.indices))
        gaps.append({"cap_charge": sector.charge,
                     "positive_modular_gap_log_arguments": sorted({str(b/a) for a in weights for b in weights if b>a})})
    return {
        "schema": "oph.native_regional_reference_preparation.v1", "all_pass": True,
        "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        "audit": exact_audit(), "preparation": preparation(),
        "raw_reference_zero_eigenvalues": raw.count(F(0)),
        "completed_reference_weight_multiplicities": {str(k): v for k, v in sorted(Counter(completed).items())},
        "regional_clock_gaps": gaps,
        "differs_from_global_movement_completed_reference": completed != cut.completed_palm(),
        "contract": {
            "architecture": "binary twelve-port total-two source",
            "laws": "one fixed C or E during local warmup; C/E complement repair",
            "schedule": "IID uniform source seams and fair placements, then uniform complement seams",
            "records": "region collects incident ticks including zeros; all final herald outcomes retained",
            "tape": "unread source marginal; no conditioning on warmup tape and no adaptive feedback",
            "finite_bound_includes_untouched_references": True,
            "passive_modular_dynamics_established": False,
            "laboratory_or_physical_spacetime_identification": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", nargs="?", const=str(OUTPUT))
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    payload = json.dumps(run(), indent=2, sort_keys=True)+"\n"
    if args.check:
        if Path(args.check).read_text(encoding="utf-8") != payload:
            raise SystemExit("regional preparation receipt differs")
        print("regional preparation receipt verified")
    else:
        args.output.write_text(payload, encoding="utf-8", newline="\n")
        print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
