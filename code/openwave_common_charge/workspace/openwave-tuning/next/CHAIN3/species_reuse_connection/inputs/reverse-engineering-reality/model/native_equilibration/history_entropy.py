"""Counted entropy and charge readout on one continuing native history.

The source has binary load two on twelve ports. A bounded, record-adaptive
prefix may cross the protected core. A two-edge census then resolves the core
from E contrasts with certainty, or C currents on a state-independent success
flag. Every native
coin, current/contrast, disturbed endpoint and timeout remains in the instrument.

For a complete terminal tape b, native branches are injective translations.
If b resolves the initial core charge q, its d_b compatible endpoints have
reference path weight w_q 2**(-T_b). Counting endpoints independently therefore
gives, for a state-independent retained flag of probability beta,

 F(p) = H(B|flag) + E_p log(d_B) - E_p[T] log(2) - log(beta),
 F(p)-S(sigma) = h1 delta<m> + h2 delta<c> - D(B_p || B_sigma).

The fixed-moment first variation vanishes; its Hessian is the negative
classical Fisher form of the actual record response. The cost subtraction is
essential even for bounded stopping times. It accounts for retained fair-coin
entropy, not a physical tape erasure. The successful census is an operational
readout, not supplied access to source occupations. Timeout fibres need not
resolve q; their full-instrument counting excess is bounded explicitly.

Completion kernels for different caps are marginals of the SAME conditional
endpoint distribution. Their overlap agreement is a finite consensus result.
Unmeasured cap occupations in that kernel are predictions, not extra records.
No area, proper volume, stress tensor, geometric boost or Einstein stationarity
is identified by this finite population theorem.

The coherent extension below computes the exact recorded-channel response
G on every source matrix unit. The Petz pullback is G entrywise rho; G=1
defines the full recoverable algebra, while G=0 defines entropy-flat output
directions. These are different questions. With the SAME two-moment sigma,
the cost-corrected quantum output entropy is stationary on all coherent
fixed-moment tangents and has Hessian -sum G_ij |X_ij|^2 L(sigma_i,sigma_j).
The all-cap checks intersect this source algebra with the existing E_R
completions. They discharge finite coherent recorded-entropy stationarity;
they supply no geometric variation, area coefficient or stress identification.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np

from . import common_variation as v
from ..native_modular_reference import local_preparation as local
from ..observer import native, quantum_mixing

TREE = ((0, 1), (0, 5))
COMPONENT = (0, 1, 5)


@dataclass(frozen=True)
class Path:
    initial: int
    endpoint: tuple[int, ...]
    tape: tuple[tuple[int, int, int, int], ...]
    chance: F
    terminal: str


def current(record, law):
    """Decode the signed transfer from the declared primitive record."""
    _, _, coin, label = record
    if law == "C":
        return label
    if law == "E":
        return (label-(2*coin-1)*(abs(label) % 2))//2
    raise ValueError("original C or E record law required")


def _step(state, tape, seam, coin, law):
    a, b = seam
    out, (_, transfer) = native.native(state, seam, coin)
    label = transfer if law == "C" else state[a]-state[b]
    return out, (*tape, (a, b, coin, label))


def adaptive_seam(tape, law):
    """A declared causal controller, using records only, not hidden loads."""
    if not tape:
        return (0, 5)
    if len(tape) % 2:
        return (0, 1) if current(tape[-1], law) else (1, 8)
    return (0, 7) if current(tape[-1], law) else (2, 10)


def edge_words(deadline):
    if type(deadline) is not int or not 2 <= deadline <= 16:
        raise ValueError("integer census deadline from 2 through 16 required")
    return tuple(((*((coin,)*(length-1)), 1-coin), "success")
                 for length in range(2, deadline+1) for coin in (0, 1)) + tuple(
                     ((coin,)*deadline, "timeout") for coin in (0, 1))


def prefix_length(tape, law, limit, stop_on_current=False):
    """Replay the declared bounded stopping decision from committed records."""
    if stop_on_current:
        for index, record in enumerate(tape[:limit]):
            if current(record, law):
                return index+1
    return limit


@lru_cache(maxsize=12)
def histories(law="C", prefix=2, deadline=3, stop_on_current=False):
    """Enumerate the complete finite instrument; no terminal branch is dropped."""
    if law not in ("C", "E") or type(prefix) is not int or not 0 <= prefix <= 5:
        raise ValueError("C/E law and a bounded zero-to-five attempt prefix required")
    words = edge_words(deadline) if law == "C" else (((0,), "success"), ((1,), "success"))
    allowed = set(quantum_mixing.seams())
    terminal = []
    for initial, state in enumerate(native.states(binary=True)):
        live, stopped = [(state, ())], []
        for _ in range(prefix):
            continuing = []
            for end, tape in live:
                if stop_on_current and tape and current(tape[-1], law):
                    stopped.append((end, tape))
                    continue
                for coin in (0, 1):
                    seam = adaptive_seam(tape, law)
                    if seam not in allowed:
                        raise ArithmeticError("controller selected a non-source seam")
                    continuing.append(_step(end, tape, seam, coin, law))
            live = continuing
        for end, tape in (*stopped, *live):
            live = [(end, tape)]
            for edge in TREE:
                continuing = []
                for end, tape in live:
                    for word, flag in words:
                        out, record = end, tape
                        for coin in word:
                            out, record = _step(out, record, edge, coin, law)
                        if flag == "timeout":
                            terminal.append(Path(initial, out, record,
                                                 F(1, 2**len(record)), "timeout"))
                        else:
                            continuing.append((out, record))
                live = continuing
            terminal.extend(Path(initial, end, tape, F(1, 2**len(tape)), "success")
                            for end, tape in live)
    return tuple(terminal)


def decode_core(tape, prefix, law="C"):
    """Recover initial core bits from a successful native tape, without a state."""
    position = prefix
    displacement = {i: 0 for i in COMPONENT}
    differences = []
    for edge in TREE:
        if position >= len(tape) or tape[position][:2] != edge:
            raise ValueError("incomplete native census")
        start = position
        if law == "E":
            contrast = tape[start][3]
        else:
            first_coin = tape[start][2]
            while position < len(tape) and tape[position][2] == first_coin:
                if tape[position][:2] != edge:
                    raise ValueError("missing native coin flip")
                position += 1
            if position >= len(tape) or tape[position][:2] != edge:
                raise ValueError("census timeout does not resolve the core")
            first = current(tape[start], law)
            flip = current(tape[position], law)
            direction = 2*first_coin-1
            contrast = first if first else (direction if flip else 0)
        a, b = edge
        differences.append((a, b, contrast-displacement[a]+displacement[b]))
        transfer = sum(current(record, law) for record in tape[start:position+1])
        displacement[a] -= transfer
        displacement[b] += transfer
        position += 1
    if position != len(tape):
        raise ValueError("unexpected records after the completed census")
    potential = {0: 0}
    for a, b, contrast in differences:
        potential[b] = potential[a]-contrast
    # A three-port component cannot be all occupied on the binary total-two
    # source. Its lowest occupation is therefore zero, fixing the offset.
    offset = -min(potential.values())
    initial = {i: potential[i]+offset for i in COMPONENT}
    if any(x not in (0, 1) for x in initial.values()) or sum(initial.values()) > 2:
        raise ValueError("records do not describe the binary total-two source")
    core = [initial[0], initial[1]]
    for record in tape[:prefix]:
        a, b = record[:2]
        transfer = current(record, law)
        for i in (0, 1):
            core[i] -= transfer*((b == i)-(a == i))
    if any(x not in (0, 1) for x in core):
        raise ValueError("prefix records do not recover a binary initial core")
    return tuple(core)


def fibres(paths, cap=()):
    """Group actual histories, optionally refine by a predicted endpoint cap."""
    if tuple(sorted(set(cap))) != cap or any(type(i) is not int or not 0 <= i < 12 for i in cap):
        raise ValueError("sorted distinct cap indices required")
    result = defaultdict(list)
    for path in paths:
        result[(path.tape, tuple(path.endpoint[i] for i in cap))].append(path)
    for row in result.values():
        if len({p.endpoint for p in row}) != len(row) or len({p.chance for p in row}) != 1:
            raise ArithmeticError("a native tape lost its injective, equal-coin fibre")
    return dict(result)


def flag_probability(paths):
    by_input = [F() for _ in range(66)]
    for path in paths:
        by_input[path.initial] += path.chance
    if min(by_input) <= 0 or len(set(by_input)) != 1:
        raise ValueError("retained flag must have positive state-independent probability")
    return by_input[0]


def _sum_forms(*forms):
    result = defaultdict(F)
    for form in forms:
        for prime, coefficient in form.items():
            result[int(prime)] += coefficient
    return v._clean(result)


def source_entropy(weights):
    v.moments(weights)
    result = defaultdict(F)
    for weight in weights:
        if weight:
            v._log_add(result, -weight, weight)
    return v._clean(result)


def ledger(paths, weights, cap=()):
    """Exact prime-log entropy ledger; edge dimensions are counted first."""
    v.moments(weights)
    beta = flag_probability(paths)
    rows = fibres(paths, cap)
    bulk, edge, cost, normalization = (defaultdict(F) for _ in range(4))
    probabilities = {}
    mean_attempts = F()
    for key, row in rows.items():
        probability = sum((weights[p.initial]*p.chance for p in row), F())/beta
        probabilities[key] = probability
        if probability:
            v._log_add(bulk, -probability, probability)
            v._log_add(edge, probability, len(row))
        mean_attempts += probability*len(row[0].tape)
    v._log_add(cost, -mean_attempts, 2)
    v._log_add(normalization, -1, beta)
    return {"probabilities": probabilities, "bulk": v._clean(bulk),
            "edge": v._clean(edge), "coin_cost": v._clean(cost),
            "flag_normalization": v._clean(normalization),
            "functional": _sum_forms(bulk, edge, cost, normalization),
            "flag_probability": beta, "mean_native_attempts": mean_attempts}


def relative_entropy(probabilities, reference):
    result = defaultdict(F)
    if probabilities.keys() != reference.keys():
        raise ValueError("the same record outcomes are required")
    for key, probability in probabilities.items():
        if probability:
            v._log_add(result, probability, probability/reference[key])
    return v._clean(result)


def conditional_cap(row, cap):
    """The common reference's uniform conditional endpoint marginal."""
    if len({v.charges()[p.initial] for p in row}) != 1:
        raise ValueError("counted reference completion needs a charge-resolved fibre")
    result = defaultdict(F)
    for path in row:
        result[tuple(path.endpoint[i] for i in cap)] += F(1, len(row))
    return dict(result)


def value(form):
    return sum(float(coefficient)*math.log(int(prime)) for prime, coefficient in form.items())


@dataclass(frozen=True)
class QuantumBranch:
    tape: tuple
    terminal: str
    inputs: tuple[int, ...]
    outputs: tuple[int, ...]
    weight: F


@dataclass(frozen=True)
class QuantumResponse:
    """Exact data for the complete tape AND disturbed quantum endpoint."""
    law: str
    flag_probability: F
    branches: tuple[QuantumBranch, ...]
    overlap: tuple[tuple[F, ...], ...]
    recovery_classes: tuple[tuple[int, ...], ...]


def quantum_response(paths, *, law):
    """Compute the reference-Petz pullback on all 4,356 native matrix units.

    A validated branch is K_b=sqrt(w_b) U_b P_b, w_b=2**(-T_b)/beta.
    Its exact endpoint inverse is the Petz map for EVERY faithful diagonal
    source reference: the reference factors cancel on each image support.
    An arbitrary trace-preserving replacement completes the decoder outside
    these supports; no actual native output is projected or postselected.
    Thus R_sigma Gamma(X)=G*X, G_ij=sum_b w_b 1_{i,j in domain_b}, where
    * is entrywise multiplication. G is a Gram matrix with diagonal one.
    Its entries equal one precisely within identical domain signatures;
    these classes give the full equality/recovery algebra direct_sum M_C.
    The identity tensors with any untouched reference, without changing G.

    Paths must implement one causal native instrument (possibly a retained
    state-independent flag). Actual endpoints, primitive labels, coin costs
    and record-only decisions are checked rather than trusting a branch list.
    """
    if type(law) is not str or law not in ("C", "E"):
        raise ValueError("original C or E record law required")
    if type(paths) is not tuple or not paths:
        raise ValueError("a nonempty tuple of actual native paths required")
    states = native.states(binary=True)
    positions = {state: i for i, state in enumerate(states)}
    seams = set(quantum_mixing.seams())
    decisions, seen = {}, set()

    def decision(prefix, action):
        if prefix in decisions and decisions[prefix] != action:
            raise ValueError("native control must depend only on committed records")
        decisions[prefix] = action

    for path in paths:
        if (type(path) is not Path or type(path.initial) is not int
                or not 0 <= path.initial < 66 or type(path.tape) is not tuple
                or type(path.terminal) is not str or path.terminal not in ("success", "timeout")
                or type(path.chance) is not F or path.chance != F(1, 2**len(path.tape))):
            raise ValueError("original native path fields and fair-coin probability required")
        if any(type(record) is not tuple or len(record) != 4
               or any(type(x) is not int for x in record)
               or record[:2] not in seams or record[2] not in (0, 1)
               for record in path.tape):
            raise ValueError("a complete original native record required")
        key = path.initial, path.tape, path.terminal
        if key in seen:
            raise ValueError("duplicate native path")
        seen.add(key)
        state, prefix = states[path.initial], ()
        for record in path.tape:
            a, b, coin, label = record
            decision(prefix, ("continue", (a, b)))
            out, (_, transfer) = native.native(state, (a, b), coin)
            expected = transfer if law == "C" else state[a]-state[b]
            if label != expected:
                raise ValueError("record differs from the actual native branch")
            prefix, state = (*prefix, record), out
        decision(prefix, ("terminal", path.terminal))
        if (type(path.endpoint) is not tuple or path.endpoint != state
                or any(type(x) is not int for x in path.endpoint)):
            raise ValueError("disturbed endpoint differs from native execution")
    beta = flag_probability(paths)
    if beta > 1:
        raise ValueError("retained flag probability cannot exceed one")
    rows = fibres(paths)
    branches, signatures = [], [[] for _ in range(66)]
    overlap = [[F() for _ in range(66)] for _ in range(66)]
    for index, row in enumerate(rows.values()):
        row = sorted(row, key=lambda path: path.initial)
        inputs = tuple(path.initial for path in row)
        compatible = []
        for initial, state in enumerate(states):
            for a, b, coin, label in row[0].tape:
                out, (_, transfer) = native.native(state, (a, b), coin)
                if label != (transfer if law == "C" else state[a]-state[b]):
                    break
                state = out
            else:
                compatible.append(initial)
        if inputs != tuple(compatible):
            raise ValueError("retained tapes must include every compatible initial state")
        weight = row[0].chance/beta
        branches.append(QuantumBranch(row[0].tape, row[0].terminal, inputs,
                                      tuple(positions[path.endpoint] for path in row), weight))
        for i in inputs:
            signatures[i].append(index)
            for j in inputs:
                overlap[i][j] += weight
    classes = defaultdict(list)
    for i, signature in enumerate(signatures):
        classes[tuple(signature)].append(i)
    if any(overlap[i][i] != 1 or not 0 <= overlap[i][j] <= 1
           or (overlap[i][j] == 1) != (signatures[i] == signatures[j])
           for i in range(66) for j in range(66)):
        raise ArithmeticError("native quantum response lost its exact sufficiency algebra")
    return QuantumResponse(law, beta, tuple(branches), tuple(map(tuple, overlap)),
                           tuple(tuple(row) for row in classes.values()))


def quantum_region_response(response, cap):
    """Exact visible and recoverable dimensions after the existing E_R.

    On a cap containing the protected core, E_R fixes every sigma(m,c).
    Source relative-entropy equality for Gamma E_R therefore requires BOTH
    rho=E_R(rho) and rho in the tape's recovery algebra. Cap matrix-unit
    supports are disjoint, so checking all their complement copies gives
    the exact intersection and image dimensions, without a floating rank.
    The full cap is a readout only, not an extra native-mixing claim.
    """
    if (type(response) is not QuantumResponse or type(cap) is not tuple
            or any(type(i) is not int or not 0 <= i < 12 for i in cap)
            or tuple(sorted(set(cap))) != cap
            or not set(v.CORE).issubset(cap)):
        raise ValueError("an exact response and a sorted core-containing cap required")
    outside = tuple(i for i in range(12) if i not in cap)
    coordinates = {(tuple(x[i] for i in cap), tuple(x[i] for i in outside)): index
                   for index, x in enumerate(native.states(binary=True))}
    retained = visible = recoverable = 0
    for charge in range(3):
        left = sorted({a for a, _ in coordinates if sum(a) == charge})
        right = sorted({b for a, b in coordinates if sum(a) == charge})
        for a in left:
            for b in left:
                entries = [response.overlap[coordinates[a, c]][coordinates[b, c]] for c in right]
                retained += 1
                visible += any(x > 0 for x in entries)
                recoverable += all(x == 1 for x in entries)
    return {"cap_algebra_dimension": retained, "record_output_image_dimension": visible,
            "flat_cap_tangent_space_dimension": retained-visible,
            "relative_entropy_equality_algebra_dimension": recoverable,
            "common_reference_fixed_by_cap_completion": True}


def quantum_stationarity(response, weights, tangent):
    """Floating evaluation of the analytically exact coherent Hessian.

    For sigma=diag(weights) from the existing faithful two-moment family,
    F(rho)=S(Gamma rho)-E_rho[T]log2-log(beta) obeys
    F(rho)=S(sigma)+h1 delta m+h2 delta c-D(Gamma rho||Gamma sigma).
    Its fixed-moment first derivative vanishes and its Hessian is
    -sum_ij G_ij |X_ij|^2 L(sigma_i,sigma_j). Flat directions are exactly
    G_ij=0; lossless directions instead have G_ij=1. Neither condition is
    a numerical entropy/rank threshold. The parameter is the declared
    dimensionless variation parameter, not an observer or physical clock.
    """
    if type(response) is not QuantumResponse or weights != v.family(*v.moments(weights)):
        raise ValueError("exact response and faithful two-moment family required")

    def original_numeric(value):
        if (value is None or isinstance(value, (bool, np.bool_))
                or np.ma.isMaskedArray(value)):
            raise ValueError("original unmasked numeric tangent entries required")
        if isinstance(value, (tuple, list)):
            for item in value:
                original_numeric(item)
        elif isinstance(value, np.ndarray) and value.dtype.kind not in "iufc":
            raise ValueError("original numeric tangent dtype required")

    original_numeric(tangent)
    raw = np.asarray(tangent)
    if (np.ma.isMaskedArray(tangent) or raw.shape != (66, 66)
            or raw.dtype.kind not in "iufc" or not np.isfinite(raw).all()):
        raise ValueError("a finite full-matrix tangent on the 66-state carrier required")
    x = np.array(raw, dtype=complex, copy=True)
    norm = float(np.linalg.norm(x, "fro"))
    if (np.linalg.norm(x-x.conj().T, "fro") > 1e-12*max(norm, 1.)
            or not math.isfinite(norm)):
        raise ValueError("a Hermitian finite tangent required")
    diagonal = np.diag(x).real
    q = np.asarray(v.charges())
    constraints = (float(diagonal.sum()), float(q@diagonal), float((q*(q-1)/2)@diagonal))
    if any(abs(value) > 1e-12*max(norm, 1.) for value in constraints):
        raise ValueError("normalization and both initial moments must be fixed")
    w = np.asarray(list(map(float, weights)))
    if (not np.isfinite(w).all() or min(w) <= 0
            or not math.isfinite(1/float(min(w)))):
        raise ValueError("reference and entropy curvature must remain faithful and finite in execution precision")
    source_terms, output_terms, loss_terms, norm_terms = [], [], [], []
    flat = lossless = True
    for i, a in enumerate(weights):
        for j, b in enumerate(weights):
            magnitude = float(abs(x[i, j])**2)
            g = response.overlap[i][j]
            lower = float(min(a, b))
            try:
                ratio = float(abs(a-b)/min(a, b))
            except OverflowError:
                ratio = math.inf
            if ratio == 0:
                divided = 1/lower
            elif math.isfinite(ratio):
                # The exact gap may underflow even when a and b round to
                # different binary values. Scale by their relative gap;
                # never divide by float(a-b).
                divided = (math.log1p(ratio)/ratio)/lower
            else:
                divided = abs(math.log(float(a))-math.log(float(b)))/float(abs(a-b))
            source_terms.append(magnitude*divided)
            output_terms.append(float(g)*magnitude*divided)
            # Compute loss directly; subtracting two large Hessians can
            # create cancellation or NaN even on an exactly lossless block.
            loss_terms.append(float(1-g)*magnitude*divided)
            norm_terms.append(float(g)*magnitude)
            if x[i, j] != 0:
                flat &= g == 0
                lossless &= g == 1
    try:
        source, curvature, loss, weighted_norm = map(math.fsum,
            (source_terms, output_terms, loss_terms, norm_terms))
    except OverflowError as error:
        raise ValueError("entropy curvature diagnostics must remain finite in execution precision") from error
    report = {"fixed_moment_first_derivative": -float(diagonal@np.log(w/w[0])),
            "output_entropy_hessian": -curvature, "source_entropy_hessian": -source,
            "information_loss_hessian": loss,
            "visible_HS_norm_squared": weighted_norm,
            "negative_hessian_lower": weighted_norm/float(max(weights)),
            "negative_hessian_upper": weighted_norm/float(min(weights)),
            "exact_structural_flat_direction": flat,
            "exact_structural_lossless_direction": lossless,
            "floating_diagnostics_are_outward_certificates": False}
    if any(type(value) is float and not math.isfinite(value) for value in report.values()):
        raise ValueError("entropy curvature diagnostics must remain finite in execution precision")
    return report


def quantum_response_certificate(paths, *, law):
    """Compact all-source equality/stationarity certificate, not a new packet."""
    import hashlib
    import json

    response = quantum_response(paths, law=law)
    counts = Counter(value for row in response.overlap for value in row)
    class_counts = Counter(map(len, response.recovery_classes))
    region_counts, digest = Counter(), hashlib.sha256()
    for mask in range(1 << 10):
        cap = (*v.CORE, *(i+2 for i in range(10) if mask & (1 << i)))
        region = quantum_region_response(response, cap)
        region_counts[(region["cap_algebra_dimension"], region["record_output_image_dimension"],
                       region["relative_entropy_equality_algebra_dimension"])] += 1
        digest.update(json.dumps([cap, region], sort_keys=True).encode())
    fixed = sum(len(row)**2 for row in response.recovery_classes)
    if fixed != counts[F(1)] or sum(counts.values()) != 4356:
        raise ArithmeticError("native equality algebra omitted a matrix unit")
    touched = tuple(sorted({port for branch in response.branches for record in branch.tape
                            for port in record[:2]}))
    occupation_classes = defaultdict(set)
    for i, state in enumerate(native.states(binary=True)):
        occupation_classes[tuple(state[port] for port in touched)].add(i)
    return {"retained_flag_probability": response.flag_probability,
            "complete_quantum_record_blocks": len(response.branches),
            "all_native_matrix_units_checked": 4356,
            "exact_Petz_overlap_multiplicities": {str(x): n for x, n in sorted(counts.items())},
            "recovery_class_size_multiplicities": {str(x): n for x, n in sorted(class_counts.items())},
            "actual_native_ports_used": touched,
            "recovery_classes_equal_initial_occupation_fibres_on_used_ports":
                {frozenset(row) for row in response.recovery_classes}
                == {frozenset(row) for row in occupation_classes.values()},
            "relative_entropy_equality_algebra_dimension": fixed,
            "fully_destroyed_matrix_units": counts[F(0)],
            "strict_coherent_entropy_curvature_dimension": 4356-counts[F(0)],
            "fixed_moment_Hermitian_tangent_dimension": 4353,
            "strict_fixed_moment_curvature_dimension": 4353-counts[F(0)],
            "diagonal_population_recovery_is_exact": True,
            "same_recovery_for_every_faithful_diagonal_reference": True,
            "arbitrary_untouched_reference_extension": True,
            "source_selected_two_moment_family_stationary": True,
            "core_containing_cap_completions_checked": 1024,
            "cap_intersection_dimension_multiplicities": [
                {"cap_algebra": a, "output_image": b, "equality_algebra": c, "caps": count}
                for (a, b, c), count in sorted(region_counts.items())],
            "exact_cap_intersections_sha256": digest.hexdigest(),
            "physical_metric_area_volume_stress_identified": False}


def run():
    weights = local.reference()
    target = source_entropy(weights)
    results = {}
    for law in ("C", "E"):
        paths = histories(law)
        kept = tuple(p for p in paths if p.terminal == "success")
        beta = flag_probability(kept)
        if flag_probability(paths) != 1 or beta != (F(9, 16) if law == "C" else 1):
            raise ArithmeticError("native census terminal probability changed")
        for path in kept:
            if decode_core(path.tape, 2, law) != native.states(binary=True)[path.initial][:2]:
                raise ArithmeticError("native record core decoder failed")
        rows = fibres(kept)
        observed = ledger(kept, weights)
        if observed["functional"] != target:
            raise ArithmeticError("record-only counted entropy envelope failed")
        # All cap comparisons refine this same history measure. They neither
        # prepare a fresh source nor independently resample an overlap state.
        caps = ((), (0, 1), (0, 1, 5), (0, 1, 2, 5), (0, 1, 5, 8), tuple(range(12)))
        for cap in caps:
            if ledger(kept, weights, cap)["functional"] != target:
                raise ArithmeticError("same-history cap envelope changed")
        complete = ledger(paths, weights)
        excess = _sum_forms(complete["functional"], {k: -x for k, x in target.items()})
        bound = (1-beta)*math.log(float(max(weights)/min(weights)))
        if (law == "C" and not 0 < value(excess) < bound) or (law == "E" and excess):
            raise ArithmeticError("timeout ambiguity control lost its strict excess or bound")
        results[law] = {"terminal_paths": len(paths), "successful_paths": len(kept),
            "record_only_success_fibres": len(rows),
            "success_probability": beta, "timeout_probability": 1-beta,
            "mean_success_attempts": observed["mean_native_attempts"],
            "mean_complete_attempts": complete["mean_native_attempts"],
            "exact_source_entropy": target, "exact_success_functional": observed["functional"],
            "complete_counting_excess": excess,
            "complete_excess_upper_coefficient": 1-beta,
            "complete_excess_upper_log_argument": max(weights)/min(weights),
            "successful_completion_dimensions": sorted({len(row) for row in rows.values()}),
            "same_history_cap_refinements_checked": len(caps),
            "complete_coherent_stationarity": quantum_response_certificate(paths, law=law)}
        if beta < 1:
            results[law]["flagged_coherent_stationarity"] = quantum_response_certificate(kept, law=law)
    stopped = histories("E", prefix=3, stop_on_current=True)
    for path in stopped:
        count = prefix_length(path.tape, "E", 3, True)
        if decode_core(path.tape, count, "E") != native.states(binary=True)[path.initial][:2]:
            raise ArithmeticError("record-stopped prefix decoder failed")
    stopped_ledger = ledger(stopped, weights)
    if stopped_ledger["functional"] != target:
        raise ArithmeticError("record-dependent native stopping cost was lost")
    per_input_cost = [sum((p.chance*len(p.tape) for p in stopped if p.initial == i), F())
                      for i in range(66)]
    h = 24
    alpha = 1-F(1, 2**(h-1))
    finite = v.stationarity_budget("local")
    return {"laws": results, "prefix_attempts": 2, "census_deadline_tested": 3,
        "bounded_C_readout_service": {"deadline": h, "native_attempts_upper": 2+2*h,
            "timeout_probability": 1-alpha**2,
            "mean_attempts_upper": 8,
            "timeout_counting_excess_upper_log_argument": max(weights)/min(weights),
            "timeout_counting_excess_upper_coefficient": 1-alpha**2},
        "unconditional_E_readout_service": {"native_attempts": 4, "timeout_probability": 0,
            "core_census_attempts": 2, "prefix_attempts": 2},
        "record_stopped_E_service": {"prefix_attempts_upper": 3, "total_attempts_upper": 5,
            "conditional_input_mean_attempts_range": (min(per_input_cost), max(per_input_cost)),
            "mean_native_attempts": stopped_ledger["mean_native_attempts"],
            "exact_cost_corrected_functional": stopped_ledger["functional"]},
        "finite_preparation_stationarity": {
            "operator_sandwich_radius": finite["operator_sandwich_radius"],
            "fixed_moment_tangent_error_per_trace_norm_upper":
                finite["entropy_tangent_error_per_trace_norm_upper"],
            "local_reference_native_attempts_upper": finite["native_total_attempts_upper"],
            "complete_flagged_reference_halftrace_upper": finite["complete_flagged_reference_halftrace_upper"],
            "adaptive_record_channel_preserves_the_reference_sandwich": True,
            "same_gradient_bound_for_coherent_record_endpoint_tangents": True,
            "input_dependent_stopping_cost_must_be_varied": True},
        "scope": {"readout_uses_only_declared_native_records": True,
            "prefix_has_record_feedback_and_core_crossings": True,
            "all_timeout_disturbed_outputs_retained": True,
            "cap_completions_share_one_conditional_endpoint_measure": True,
            "unmeasured_cap_predictions_are_additional_records": False,
            "coherent_source_and_untouched_reference_information_identity": True,
            "coherent_stationarity_here_varies_the_source_state": True,
            "coherent_source_fixed_moment_stationarity": True,
            "exact_record_endpoint_quantum_sufficiency_algebra": True,
            "proper_volume_area_stress_dictionary": False}}
