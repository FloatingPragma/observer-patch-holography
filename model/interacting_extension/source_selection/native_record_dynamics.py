#!/usr/bin/env python3
"""Complete native repair records, their conditioned dynamics, and native clocks.

The unchanged r2062 endpoint rule and its exact thirty-seam chart are used.
Uniform seam choice, the C/E record choice, a positive activity tilt, coherent
record access, and any Euclidean continuation are distinguished explicitly.
No particle, mixing, measured energy, or physical-time input is imported.

Receiving API:
    native_instrument(repo, candidate="C")
    doob_instrument(native, z)                 # 0 < z <= 1, z=exp(-s)
    complete_outcomes(instrument, rho)         # every label, even zero ones
    apply(instrument, rho)
    count_distribution(instrument, rho, n)     # actual recorded move count
    filter_implementation(native, z, rho)      # success AND all failures
    native_blocks(native, z)                   # all 4356 operator directions
    record_references(native, rho, repo)       # all chart-transducer outcomes
    two_record_reference_law(native, repo)     # exact full joint raw-record law

This is a conditional native-law extension, not a claim that its tilt or
Euclidean unitary is selected physical dynamics. See PROOF.md.
"""
from __future__ import annotations

import ast
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

from .._source import verify_sources

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PIN = "3cf3c14810671a16e61b7915f69079a1773a12a2"


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def is_zero(matrix):
    return all(sp.simplify(value) == 0 for value in matrix)


@lru_cache(maxsize=2)
def _geometry_from_text(text):
    """Read the verified literal chart without constructing family services."""
    tree = ast.parse(text)
    function = next(node for node in tree.body
                    if isinstance(node, ast.FunctionDef) and node.name == 'raw_frame')
    expression = next(node.value for node in function.body if isinstance(node, ast.Return))
    if not all(isinstance(node, (ast.Tuple, ast.Name, ast.Load)) for node in ast.walk(expression)):
        raise ValueError('source frame does not have the literal chart grammar')
    phi = (1+sp.sqrt(5))/2
    names = dict(zero=sp.Integer(0), one=sp.Integer(1), minus=sp.Integer(-1),
                 phi=phi, negative_phi=-phi)
    frame = sp.ImmutableMatrix(eval(compile(ast.Expression(expression), '<verified chart>', 'eval'),
                                    {'__builtins__': {}}, names))
    if frame.shape != (12, 3):
        raise ValueError('native chart dimension differs')
    projectors = tuple(sp.ImmutableMatrix((frame[j, :].T*frame[j, :]/
                       (frame[j, :]*frame[j, :].T)[0]).applyfunc(sp.simplify))
                       for j in range(12))
    seams = tuple((a, b) for a in range(12) for b in range(a+1, 12)
                  if sp.simplify(sum((frame[a, k]-frame[b, k])**2 for k in range(3))) == 4)
    if (len(seams) != 30 or not is_zero(sum(projectors, sp.zeros(3))-4*sp.eye(3))
            or any(not is_zero(P*P-P) or sp.trace(P) != 1 for P in projectors)):
        raise ArithmeticError('native chart frame identities differ')
    return frame, projectors, seams


def geometry_packet(repo):
    """Verify live source bytes even when the exact chart is cached."""
    repo = Path(repo).resolve()
    name = 'model/common_source/slow_geometry.py'
    hashes = verify_sources(repo, (name,))
    frame, projectors, seams = _geometry_from_text((repo/name).read_text())
    return dict(frame=frame, projectors=projectors, seams=seams, source_hashes=hashes)


@dataclass(frozen=True)
class Instrument:
    labels: tuple
    destinations: np.ndarray   # (all labels, input); -1 = zero column
    amplitudes: np.ndarray     # same shape, zero on invalid columns
    costs: np.ndarray         # recorded |current|, constant on a label
    states: tuple
    seams: tuple
    metadata: dict

    @property
    def dimension(self):
        return len(self.states)


def _rho(value, d):
    a = np.asarray(value, dtype=np.complex128)
    if a.shape != (d, d) or not np.isfinite(a).all():
        raise ValueError("a finite full native density matrix is required")
    if np.max(np.abs(a-a.conj().T)) > 2e-12:
        raise ValueError("density matrix is not Hermitian")
    if np.linalg.eigvalsh(a).min() < -2e-12:
        raise ValueError("density matrix is not positive")
    if abs(np.trace(a)-1) > 2e-12:
        raise ValueError("density matrix must have trace one")
    return a


def native_instrument(repo, candidate="C", *, binary=True):
    """Literal endpoint calls; no assumed transition matrix substitutes for them.

    Source C records current; source E additionally resolves old contrast.
    The two instruments have equal population dynamics and different
    coherence backaction. The binary shell is an exact invariant restriction.
    """
    if candidate not in ("C", "E"):
        raise ValueError("candidate must be C or E")
    if type(binary) is not bool:
        raise ValueError("binary must be an explicit Boolean")
    repo = Path(repo).resolve()
    source_hashes = verify_sources(repo, ("model/observer/native.py", "model/observer/quantum_mixing.py",
                                         "model/observer/source_discrimination/circuit.py"))
    packet = geometry_packet(repo)
    native = _load(repo / "model/observer/native.py", "oph_native_record_source")
    states = tuple(native.states(binary=binary))
    index = {state: k for k, state in enumerate(states)}
    seams = tuple(packet["seams"])
    # Read the actual raw-code convention. Reserved pointer codes remain
    # explicit zero branches even though the invariant shell cannot attain
    # them. This is the exact finite instrument, not a new pulse replay.
    circuit_path = repo / "model/observer/source_discrimination/circuit.py"
    circuit_ast = ast.parse(circuit_path.read_text())
    constants = {n.targets[0].id: ast.literal_eval(n.value)
                 for n in circuit_ast.body if isinstance(n, ast.Assign)
                 and len(n.targets) == 1 and isinstance(n.targets[0], ast.Name)
                 and n.targets[0].id in ("CURRENT_CODES", "RECORD_BITS")}
    current_codes = constants["CURRENT_CODES"]
    record_size = 2**constants["RECORD_BITS"][candidate]
    if current_codes != {0: 0, 1: 1, -1: 2} or record_size not in (4, 8):
        raise ArithmeticError("source raw record grammar changed")
    labels = tuple((e, c, r) for e in range(len(seams))
                   for c in (0, 1) for r in range(record_size))
    label_index = {label: i for i, label in enumerate(labels)}
    dest = np.full((len(labels), len(states)), -1, dtype=np.int64)
    amp = np.zeros_like(dest, dtype=float)
    costs = np.zeros(len(labels), dtype=np.int64)
    observed_costs = defaultdict(set)
    for x, state in enumerate(states):
        for e, seam in enumerate(seams):
            for coin in (0, 1):
                out, record = native.native(state, seam, coin)
                current = record[1]
                contrast = state[seam[0]]-state[seam[1]]
                if out not in index or record != (coin, state[seam[0]]-out[seam[0]]):
                    raise ArithmeticError("literal native record or shell has changed")
                r = current_codes[current]
                if candidate == "E":
                    r += 4*(abs(contrast) % 2)
                k = label_index[e, coin, r]
                dest[k, x] = index[out]
                amp[k, x] = 1/math.sqrt(60)
                observed_costs[k].add(abs(current))
    for k, label in enumerate(labels):
        e, coin, raw = label
        # Reserved current code 3 has no native branch and contributes zero
        # to every count probability; its reported cost is conventionally 0.
        expected = int((raw % 4) in (1, 2))
        if observed_costs[k] and observed_costs[k] != {expected}:
            raise ArithmeticError("record does not determine the selected activity")
        costs[k] = expected
        attained = dest[k, dest[k] >= 0]
        if len(set(map(int, attained))) != len(attained):
            raise ArithmeticError("native branch is not injective on its support")
    if not np.all(np.count_nonzero(amp, axis=0) == 60):
        raise ArithmeticError("not every source event is retained")
    hashes = dict(packet["source_hashes"])
    hashes.update(source_hashes)
    hashes["model/observer/quantum_mixing.py"] = hashlib.sha256(
        (repo / "model/observer/quantum_mixing.py").read_bytes()).hexdigest()
    hashes["model/observer/source_discrimination/circuit.py"] = hashlib.sha256(circuit_path.read_bytes()).hexdigest()
    return Instrument(labels, dest, amp, costs, states, seams, {
        "source_commit": PIN, "source_hashes": hashes,
        "candidate": candidate, "native_quantum_record_convention": "real positive lift",
        "seam_policy": "uniform independent choices on all thirty native seams",
        "seam_policy_uniqueness": "within memoryless seam-only policies under transitive carrier covariance",
        "native_endpoint_rule_unchanged": True,
        "binary_shell_invariant": True, "input_restricted_to_binary_shell": binary,
        "zero_labels_retained": True,
        "raw_pointer_dimension_per_seam_coin": record_size,
        "phase_exact_source_compiler_rerun_in_this_extension": False,
        "physical_state_or_time_selected": False,
    })


def _branch_output(inst, r, rho):
    x = np.flatnonzero(inst.destinations[r] >= 0)
    y = inst.destinations[r, x]
    a = inst.amplitudes[r, x]
    out = np.zeros_like(rho, dtype=np.complex128)
    out[np.ix_(y, y)] = a[:, None] * rho[np.ix_(x, x)] * a.conj()[None, :]
    return out


def complete_outcomes(inst, rho):
    """Subnormalized output states for every label, including impossible labels."""
    rho = _rho(rho, inst.dimension)
    return tuple(_branch_output(inst, r, rho) for r in range(len(inst.labels)))


def apply(inst, rho):
    """Linear channel; permits subnormalized intermediate operators."""
    rho = np.asarray(rho, dtype=np.complex128)
    if rho.shape != (inst.dimension, inst.dimension):
        raise ValueError("wrong native operator dimension")
    out = np.zeros_like(rho)
    for r in range(len(inst.labels)):
        out += _branch_output(inst, r, rho)
    return out


def effects_diagonal(inst):
    # Every native branch is injective on its support, so K* K is diagonal.
    return np.abs(inst.amplitudes)**2


def population_matrix(inst):
    out = np.zeros((inst.dimension, inst.dimension))
    for r in range(len(inst.labels)):
        x = np.flatnonzero(inst.destinations[r] >= 0)
        out[inst.destinations[r, x], x] += abs(inst.amplitudes[r, x])**2
    return out


def distance_orbits(inst):
    if inst.dimension != 66:
        raise ValueError("the distance-orbit Doob construction requires the invariant binary shell")
    d = np.full((12, 12), 99, dtype=int)
    np.fill_diagonal(d, 0)
    for a, b in inst.seams:
        d[a, b] = d[b, a] = 1
    for k in range(12):
        d = np.minimum(d, d[:, k, None]+d[None, k, :])
    orbit = np.asarray([d[tuple(i for i, n in enumerate(state) if n)]
                        for state in inst.states], dtype=int)
    if Counter(map(int, orbit)) != {1: 30, 2: 30, 3: 6}:
        raise ArithmeticError("native pair-distance orbits changed")
    return orbit


def require_native(inst, *, binary=True):
    """Reject tilted/weighted channels in exact equal-event counting APIs."""
    if not isinstance(inst, Instrument):
        raise ValueError("a complete native Instrument is required")
    if binary and inst.dimension != 66:
        raise ValueError("the exact count formula requires the binary native shell")
    if (inst.destinations.ndim != 2 or inst.destinations.dtype.kind not in 'iu'
            or inst.costs.shape != (len(inst.labels),) or inst.costs.dtype.kind not in 'iu'
            or inst.destinations.shape != (len(inst.labels), inst.dimension)):
        raise ValueError("complete native integer destinations and activities are required")
    valid = inst.destinations >= 0
    if (inst.destinations.shape != inst.amplitudes.shape
            or not np.isfinite(inst.amplitudes).all()
            or np.any(inst.destinations[valid] >= inst.dimension)
            or not np.all(np.count_nonzero(valid, axis=0) == 60)
            or not np.all(inst.amplitudes[~valid] == 0)
            or not np.all(inst.amplitudes[valid] == 1/math.sqrt(60))):
        raise ValueError("this formula requires unchanged native equal-event amplitudes")
    if not set(inst.costs.tolist()) <= {0, 1}:
        raise ValueError("native activity must be zero or one")


def quotient_counts(inst):
    """Exact matrices; column = old pair-distance, row = new pair-distance."""
    require_native(inst)
    orbit = distance_orbits(inst)
    tables = {}
    for x in range(inst.dimension):
        q = np.zeros((2, 3), dtype=int)
        for r in range(len(inst.labels)):
            y = inst.destinations[r, x]
            if y >= 0:
                q[inst.costs[r], orbit[y]-1] += 1
        old = int(orbit[x])
        if old in tables and not np.array_equal(tables[old], q):
            raise ArithmeticError("distance is not an equitable native partition")
        tables[old] = q
    return tuple(np.stack([tables[j][a] for j in (1, 2, 3)], axis=1)
                 for a in (0, 1))


def doob_instrument(native, z):
    """Activity-conditioned quantum dynamics, with all native record labels.

    z is explicit and is not identified with public P or record capacity N.
    The derived Kraus factors change the instrument; a complete heralded
    realization from the unchanged native instrument is supplied separately.
    """
    if isinstance(z, (bool, np.bool_)) or not np.isfinite(z) or not 0 < z <= 1:
        raise ValueError("a finite explicit tilt 0 < z <= 1 is required")
    z = float(z)
    orbit = distance_orbits(native)
    stay, move = quotient_counts(native)
    sizes = np.array([30., 30., 6.])
    M = (stay+z*move)/60
    # M acts on orbit probabilities; this conjugation is the orthonormal
    # orbit basis of the symmetric microscopic tilted transition operator.
    B = M * np.sqrt(sizes[None, :]/sizes[:, None])
    if np.max(abs(B-B.T)) > 2e-15:
        raise ArithmeticError("tilted native population operator is not symmetric")
    values, vectors = np.linalg.eigh(B)
    lam, u = float(values[-1]), vectors[:, -1]
    if u.sum() < 0:
        u = -u
    h_orbit = u/np.sqrt(sizes)
    h_orbit /= h_orbit.max()  # harmless Perron-vector normalization
    if z == 1:
        # This is an exact source identity, not a small-eigenvalue clipping.
        lam, h_orbit = 1., np.ones(3)
    h = h_orbit[orbit-1]
    if np.min(h) <= 0:
        raise ArithmeticError("positive Perron effect was not recovered")
    amps = native.amplitudes.copy()
    for r in range(len(native.labels)):
        x = np.flatnonzero(native.destinations[r] >= 0)
        y = native.destinations[r, x]
        amps[r, x] *= np.sqrt(z**int(native.costs[r])*h[y]/h[x]/lam)
    stationary = np.diag(h*h/np.dot(h, h))
    metadata = dict(native.metadata)
    metadata.update({"activity_tilt_z": z, "activity_tilt_s": -math.log(z),
                     "tilt_source_selected": False, "P_or_N_to_tilt_map_derived": False,
                     "instrument_changed_by_conditioning": z != 1,
                     "physical_state_or_time_selected": False})
    driven = Instrument(native.labels, native.destinations.copy(), amps,
                        native.costs.copy(), native.states, native.seams, metadata)
    return {"instrument": driven, "stationary_state": stationary,
            "perron_effect": h, "perron_effect_by_distance": h_orbit,
            "perron_value": lam, "symmetric_orbit_transfer": B,
            "native_orbit_population_transfer": M, "orbit": orbit,
            "herald_success_probability": lam*float(h.min()/h.max())}


def count_distribution(inst, rho, n):
    """Complete move-count output; summing count blocks gives the channel power."""
    if type(n) is not int or n < 0:
        raise ValueError("nonnegative integer attempt count required")
    rho = _rho(rho, inst.dimension)
    blocks = [rho]
    for _ in range(n):
        following = [np.zeros_like(rho) for _ in range(len(blocks)+1)]
        for k, block in enumerate(blocks):
            for r in range(len(inst.labels)):
                following[k+int(inst.costs[r])] += _branch_output(inst, r, block)
        blocks = following
    probabilities = np.array([np.trace(block).real for block in blocks])
    return {"outcomes": list(range(n+1)), "probabilities": probabilities,
            "states": tuple(blocks)}


def filter_implementation(native, z, rho):
    """Operational native+filter implementation; never discards failed records.

    Stages: inverse-Perron input filter, unchanged native repair, recorded
    activity filter, Perron output filter. All stage failures remain in the
    complete instrument. Conditional successful output equals the Doob map.
    """
    rho = _rho(rho, native.dimension)
    result = doob_instrument(native, z)
    h = result["perron_effect"]
    lo, hi = h.min(), h.max()
    pre = np.sqrt(lo/h)
    post = np.sqrt(h/hi)
    pre_fail = np.sqrt(np.maximum(0., 1-lo/h))
    post_fail = np.sqrt(np.maximum(0., 1-h/hi))
    prepared = pre[:, None]*rho*pre[None, :]
    outcomes = [(('pre_failure',), pre_fail[:, None]*rho*pre_fail[None, :])]
    successful = np.zeros_like(rho)
    for r, label in enumerate(native.labels):
        b = _branch_output(native, r, prepared)
        weight = z**int(native.costs[r])
        outcomes.append((('activity_failure', *label), (1-weight)*b))
        outcomes.append((('post_failure', *label),
                         weight*post_fail[:, None]*b*post_fail[None, :]))
        good = weight*post[:, None]*b*post[None, :]
        outcomes.append((('success', *label), good))
        successful += good
    total = sum((np.trace(b).real for _, b in outcomes), 0.)
    alpha = result["herald_success_probability"]
    target = alpha*apply(result["instrument"], rho)
    return {"outcomes": tuple(outcomes), "total_probability": total,
            "success_probability": float(np.trace(successful).real),
            "predicted_success_probability": alpha,
            "successful_output": successful,
            "successful_output_residual": float(np.max(abs(successful-target)))}


def native_blocks(native, z=1.):
    """Complete tilted superoperator in its exact signed-difference blocks.

    Returns (difference, ordered matrix-unit pairs, zero-counts, move-counts).
    Blocks are (zero_counts + z*move_counts)/60. No small eigenvalue or mode
    is dropped. The source's binary restriction has 4356 matrix units.
    """
    require_native(native)
    states = native.states
    groups = defaultdict(list)
    for i, left in enumerate(states):
        for j, right in enumerate(states):
            groups[tuple(a-b for a, b in zip(left, right, strict=True))].append((i, j))
    result = []
    for delta, pairs in sorted(groups.items()):
        positions = {pair: i for i, pair in enumerate(pairs)}
        counts = np.zeros((2, len(pairs), len(pairs)), dtype=np.int64)
        for col, (i, j) in enumerate(pairs):
            for r in range(len(native.labels)):
                a, b = native.destinations[r, i], native.destinations[r, j]
                if a >= 0 and b >= 0:
                    if (a, b) not in positions:
                        raise ArithmeticError("native record changed signed difference")
                    counts[native.costs[r], positions[a, b], col] += 1
        if not all(np.array_equal(v, v.T) for v in counts):
            raise ArithmeticError("tilted source block is not exactly self-adjoint")
        result.append((delta, tuple(pairs), counts[0], counts[1]))
    if sum(len(pairs) for _, pairs, _, _ in result) != 66**2:
        raise ArithmeticError("an operator direction is missing")
    return tuple(result)


def coherent_record_read(inst, psi, record_phases=None, commit_first=False):
    """Every Fourier-pointer outcome of a declared coherent record read.

    Seam and coin choices stay committed. Only each repair pointer is read
    coherently before its own commitment. After commitment, diagonal record
    laws and all count clocks ignore phases. Every pointer outcome remains.
    """
    psi = np.asarray(psi, dtype=np.complex128)
    if psi.shape != (inst.dimension,) or abs(np.vdot(psi, psi)-1) > 2e-12:
        raise ValueError("normalized full native vector required")
    m = len(inst.labels)
    phases = np.ones(m, complex) if record_phases is None else np.asarray(record_phases, complex)
    if phases.shape != (m,) or np.max(abs(abs(phases)-1)) > 2e-12:
        raise ValueError("one unit-modulus phase per full record label required")
    columns = np.zeros((inst.dimension, m), complex)
    for r in range(m):
        x = np.flatnonzero(inst.destinations[r] >= 0)
        columns[inst.destinations[r, x], r] = phases[r]*inst.amplitudes[r, x]*psi[x]
    committed = np.sum(abs(columns)**2, axis=0)
    raw = inst.metadata["raw_pointer_dimension_per_seam_coin"]
    if commit_first:
        # No interference between already committed seam/coin choices.
        probabilities = np.repeat(committed.reshape(60, raw).sum(axis=1)/raw, raw)
    else:
        rotated = np.fft.fft(columns.reshape(inst.dimension, 60, raw), axis=2)/math.sqrt(raw)
        probabilities = np.sum(abs(rotated)**2, axis=0).reshape(m)
    return {"committed_probabilities": committed,
            "pointer_probabilities": probabilities,
            "activity_probabilities": np.array([committed[inst.costs == a].sum() for a in (0, 1)])}


def _reference_geometry(native, repo):
    data = geometry_packet(repo)
    projectors = np.asarray([np.array(P.evalf(30), complex) for P in data["projectors"]])
    frame = np.array(data["frame"].evalf(30), float)
    vectors = frame/np.linalg.norm(frame, axis=1)[:, None]
    occupation = np.asarray(native.states, dtype=np.int64)
    per_state_family = np.einsum("xi,iab->xab", occupation/2, projectors)
    return data, projectors, vectors, occupation, per_state_family


def record_references(native, rho, repo):
    """Complete source-chart transducer following a native recorded repair.

    L_(x,i)=sqrt(n_i(x)/2)|v_i><x| is the specified measurement/preparation
    service; its output label is (native raw record, final load x, selected
    occupied port i). All 66*12 labels, including zeros, remain for every
    raw record. This destroys native coherences. A retained-copy version
    keeps the measured |x> register as an additional output.

    F_r is the reduced family state conditioned only on native record r.
    Conditioning additionally on the transducer's x,i makes its family
    branch pure. This distinction is part of the output, not hidden.
    """
    rho = _rho(rho, native.dimension)
    _, projectors, vectors, occupation, per_state = _reference_geometry(native, repo)
    posterior_diagonals = np.zeros((len(native.labels), native.dimension))
    for r in range(len(native.labels)):
        x = np.flatnonzero(native.destinations[r] >= 0)
        posterior_diagonals[r, native.destinations[r, x]] = (
            abs(native.amplitudes[r, x])**2*np.diag(rho).real[x])
    full_probabilities = posterior_diagonals[:, :, None]*occupation[None, :, :]/2
    reduced_outputs = np.einsum("rx,xab->rab", posterior_diagonals, per_state)
    return {"raw_record_labels": native.labels,
            "transducer_final_native_labels": native.states,
            "transducer_port_labels": tuple(range(12)),
            "complete_probabilities": full_probabilities,
            "family_vectors_for_each_port_outcome": vectors,
            "raw_record_probabilities": posterior_diagonals.sum(axis=1),
            "subnormalized_family_outputs": reduced_outputs,
            "native_chart_projectors": projectors,
            "transducer_is_a_declared_readout_preparation_service": True,
            "full_outcome_family_states_are_pure": True,
            "reduced_raw_record_conditioned_reference_is_autonomous_vacuum": False}


def two_record_reference_law(native, repo):
    """Exact rational joint law and reference numerators for two native events.

    Uses I_66/66 and the unchanged source instrument. Every branch is a
    scaled injection, so a history's probability is its integer support
    count divided by 66*60**2. No probabilities are rounded or truncated.
    """
    require_native(native)
    _, projectors, _, occupation, _ = _reference_geometry(native, repo)
    m, d = native.destinations.shape
    input_support = (native.destinations >= 0).astype(np.int64)
    output_support = np.zeros((m, d), dtype=np.int64)
    for r in range(m):
        x = np.flatnonzero(input_support[r])
        output_support[r, native.destinations[r, x]] = 1
    counts = output_support@input_support.T
    port_counts = np.zeros((m, m, 12), dtype=np.int64)
    for port in range(12):
        final_occupied = np.zeros((m, d), dtype=np.int64)
        for r in range(m):
            x = np.flatnonzero(input_support[r])
            final_occupied[r, x] = occupation[native.destinations[r, x], port]
        port_counts[:, :, port] = output_support@final_occupied.T
    denominator = 66*60**2
    if counts.sum() != denominator or not np.array_equal(port_counts.sum(axis=2), 2*counts):
        raise ArithmeticError("complete two-record Born normalization failed exactly")
    first = input_support.sum(axis=1)
    if not np.array_equal(counts.sum(axis=1), 60*first):
        raise ArithmeticError("second-event marginal changed the first record law")
    # Fixed-time appended records preserve the prefix. The invariant source
    # state also gives the same one-event second marginal.
    if not np.array_equal(counts.sum(axis=0), 60*first):
        raise ArithmeticError("stationary second-event raw law changed")
    subnormalized = np.einsum("rsi,iab->rsab", port_counts/(2*denominator), projectors)
    return {"raw_record_labels": native.labels,
            "probability_numerators": counts, "probability_denominator": denominator,
            "reference_port_numerators": port_counts,
            "subnormalized_family_outputs": subnormalized,
            "history_reference_is_derived_from_complete_posterior": True,
            "zero_probability_histories_retained": True}


def reference_audit(native, repo, archive):
    d = native.dimension
    transducer = record_references(native, np.eye(d)/d, repo)
    p = transducer["raw_record_probabilities"]
    F = transducer["subnormalized_family_outputs"]
    s145 = math.sqrt(145)
    expected = {10: np.array([(20-s145)/50, 1/5, (20+s145)/50]),
                56: np.sort(np.array([(90-s145)/280, 5/14, (90+s145)/280])),
                46: np.sort(np.array([(7-2/math.sqrt(5))/23,
                                      (7+2/math.sqrt(5))/23, 9/23]))}
    spectra = np.zeros((len(p), 3))
    errors = []
    support_counts = np.count_nonzero(native.amplitudes, axis=1)
    group_counts = Counter(map(int, support_counts))
    for r in np.flatnonzero(p > 0):
        spectra[r] = np.linalg.eigvalsh(F[r]/p[r])
        errors.append(float(np.max(abs(spectra[r]-expected[int(support_counts[r])]))))
    history = two_record_reference_law(native, repo)
    counts = history["probability_numerators"]
    selected = counts > 0
    history_states = history["subnormalized_family_outputs"][selected]
    history_prob = counts[selected]/history["probability_denominator"]
    history_spectra = np.linalg.eigvalsh(history_states/history_prob[:, None, None])
    # Only a diagnostic grouping; exact integer port numerators are saved.
    distinct_spectra = sorted(set(tuple(np.round(v, 11)) for v in history_spectra))
    data, P, vectors, occupation, _ = _reference_geometry(native, repo)
    a0, b0 = native.seams[0]
    B0 = np.column_stack((vectors[a0], vectors[b0], np.cross(vectors[a0], vectors[b0])))
    rotation_error = 0.
    for a, b in tuple(native.seams)+tuple((b, a) for a, b in native.seams):
        B = np.column_stack((vectors[a], vectors[b], np.cross(vectors[a], vectors[b])))
        R = B@np.linalg.inv(B0)
        rotation_error = max(rotation_error, float(np.max(abs(R.T@R-np.eye(3)))),
                             abs(float(np.linalg.det(R))-1),
                             float(np.max(abs(R@P[a0]@R.T-P[a]))),
                             float(np.max(abs(R@P[b0]@R.T-P[b]))))
        moved = vectors@R.T
        if max(np.min(np.linalg.norm(vectors-v, axis=1)) for v in moved) > 3e-14:
            raise ArithmeticError("an oriented-edge frame map is not a native carrier rotation")
    # Exact characteristic polynomials for every orbit, independent of a
    # floating eigenvalue grouping. Adjacent source rays overlap by 1/5.
    Pa, Pb = data["projectors"][a0], data["projectors"][b0]
    x = sp.symbols("x")
    exact_matrices = {
        "one_occupied_endpoint": sp.eye(3)/5+sp.Rational(9,20)*Pb-sp.Rational(1,20)*Pa,
        "C_zero_current": (40*sp.eye(3)-9*Pa+Pb)/112,
        "E_even_contrast": (9*sp.eye(3)-2*(Pa+Pb))/23,
    }
    exact_spectra = {
        "one_occupied_endpoint": [(20-sp.sqrt(145))/50, sp.Rational(1,5), (20+sp.sqrt(145))/50],
        "C_zero_current": [(90-sp.sqrt(145))/280, sp.Rational(5,14), (90+sp.sqrt(145))/280],
        "E_even_contrast": [(7-2/sp.sqrt(5))/23, (7+2/sp.sqrt(5))/23, sp.Rational(9,23)]}
    exact_characters = {}
    for label, matrix in exact_matrices.items():
        polynomial = sp.factor((x*sp.eye(3)-matrix).det())
        difference = sp.simplify(polynomial-sp.prod(x-y for y in exact_spectra[label]))
        if difference != 0:
            raise ArithmeticError("exact conditional reference spectrum failed")
        exact_characters[label] = {"characteristic_polynomial": str(polynomial),
                                   "eigenvalues": [str(v) for v in exact_spectra[label]]}
    # Six antipodal ray lines: any three distinct lines are independent.
    # Their 3x3 Gram determinants are 2/5 +/- 2/(5 sqrt(5)), both positive.
    axis_groups = []
    unseen = set(range(12))
    while unseen:
        i = min(unseen)
        group = tuple(j for j in sorted(unseen)
                      if all(sp.simplify(v) == 0 for v in data["projectors"][i]-data["projectors"][j]))
        if len(group) != 2:
            raise ArithmeticError("native antipodal ray pairing changed")
        axis_groups.append(group)
        unseen.difference_update(group)
    used_axes = np.stack([history["reference_port_numerators"][:, :, list(group)].sum(axis=2) > 0
                          for group in axis_groups], axis=2).sum(axis=2)
    exact_rank = np.minimum(3, used_axes)
    rank_histories = {str(k): int(np.count_nonzero((exact_rank == k) & selected)) for k in (1,2,3)}
    rank_probability_numerators = {str(k): int(counts[exact_rank == k].sum()) for k in (1,2,3)}
    if sum(rank_probability_numerators.values()) != history["probability_denominator"]:
        raise ArithmeticError("exact conditional rank partition lost probability")
    candidate = native.metadata["candidate"]
    archive[candidate+"_transducer_complete_probabilities"] = transducer["complete_probabilities"]
    archive[candidate+"_first_record_family_substates"] = F
    archive[candidate+"_first_record_family_spectra"] = spectra
    archive[candidate+"_two_record_probability_numerators"] = counts
    archive[candidate+"_two_record_reference_port_numerators"] = history["reference_port_numerators"]
    archive[candidate+"_two_record_family_spectra"] = history_spectra
    archive[candidate+"_two_record_family_exact_ranks"] = exact_rank
    return {
        "complete_transducer_outcome_count": int(np.prod(transducer["complete_probabilities"].shape)),
        "transducer_total_probability_residual": abs(float(transducer["complete_probabilities"].sum())-1),
        "raw_record_support_size_counts_including_zero_labels": dict(sorted(group_counts.items())),
        "unconditional_family_state_residual_to_I_over_3": float(np.max(abs(F.sum(axis=0)-np.eye(3)/3))),
        "all_nonzero_first_record_reference_spectra_max_residual": max(errors),
        "all_oriented_edge_native_rotation_residual": rotation_error,
        "exact_reference_orbit_spectra": exact_characters,
        "two_record_probability_denominator": history["probability_denominator"],
        "two_record_complete_probability_numerator_sum": int(counts.sum()),
        "two_record_all_record_pairs": int(counts.size),
        "two_record_nonzero_record_pairs": int(np.count_nonzero(counts)),
        "two_record_conditional_spectra_groups_at_11_decimals_diagnostic": len(distinct_spectra),
        "two_record_conditional_eigenvalue_minimum": float(history_spectra.min()),
        "two_record_conditional_reference_exact_rank_counts": rank_histories,
        "two_record_conditional_reference_rank_probability_numerators": rank_probability_numerators,
        "two_record_rank_probability_denominator": history["probability_denominator"],
        "small_negative_float_eigenvalues_are_roundoff_of_exact_positive_sums": True,
        "rank_deficient_histories_must_not_use_a_full_space_modular_log": True,
        "two_record_unconditional_family_residual_to_I_over_3": float(np.max(abs(
            history["subnormalized_family_outputs"].sum(axis=(0,1))-np.eye(3)/3))),
        "full_quantum_outcome_family_states_are_pure": True,
        "faithful_reference_conditions_only_on_native_record_and_traces_transducer_records": True,
        "chart_measurement_preparation_transducer_declared": True,
        "reference_selected_as_an_autonomous_vacuum": False,
    }


def _f(x):
    return float(x)


def audit(repo, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    report = {"schema": "oph.native_record_dynamics.v1", "source_commit": PIN,
              "no_empirical_mass_mixing_or_energy_inputs": True, "candidates": {}}
    archive = {}
    rng = np.random.default_rng(271828)
    for candidate in ("C", "E"):
        native = native_instrument(repo, candidate)
        full_native = native_instrument(repo, candidate, binary=False)
        full_binary = np.array([max(x) == 1 for x in full_native.states])
        transient_counts = np.zeros(78, dtype=int)
        for row in full_native.destinations:
            x = np.flatnonzero(row >= 0)
            transient_counts[x] += ~full_binary[row[x]]
        if not np.array_equal(transient_counts, 50*(~full_binary).astype(int)):
            raise ArithmeticError("the exact full-carrier transient survival law changed")
        d, m = native.dimension, len(native.labels)
        v = rng.normal(size=d)+1j*rng.normal(size=d)
        v /= np.linalg.norm(v)
        rho = np.outer(v, v.conj())
        source_outputs = complete_outcomes(native, rho)
        source_probability = np.array([np.trace(x).real for x in source_outputs])
        stay, move = quotient_counts(native)
        blocks = native_blocks(native)
        summary = {
            "metadata": native.metadata, "dimension": d, "all_record_labels": m,
            "full_native_carrier_dimension": full_native.dimension,
            "full_native_carrier_tp_residual": float(max(abs(effects_diagonal(full_native).sum(axis=0)-1))),
            "full_native_carrier_nonzero_record_labels": int(np.count_nonzero(np.any(full_native.amplitudes,axis=1))),
            "full_native_carrier_all_event_columns": int(np.count_nonzero(full_native.amplitudes)),
            "full_native_double_load_survival_probability_after_n_attempts": "(5/6)**n times initial double-load population",
            "full_native_to_binary_population_law_exact": True,
            "nonzero_record_labels": int(np.count_nonzero(np.any(native.amplitudes, axis=1))),
            "all_native_event_columns": int(np.count_nonzero(native.amplitudes)),
            "tp_residual": _f(max(abs(effects_diagonal(native).sum(axis=0)-1))),
            "unital_residual": _f(np.max(abs(apply(native, np.eye(d))-np.eye(d)))),
            "arbitrary_complex_input_normalization_residual": _f(abs(source_probability.sum()-1)),
            "orbit_sizes": [30, 30, 6], "zero_activity_quotient_counts": stay.tolist(),
            "one_activity_quotient_counts": move.tolist(),
            "operator_block_size_counts": dict(sorted(Counter(len(b[1]) for b in blocks).items())),
            "all_operator_directions": sum(len(b[1]) for b in blocks),
            "zero_tilt_exact_gershgorin_lower_bound_numerator_over_60": int(min(
                (np.diag(a+b)-(np.abs(a+b).sum(axis=1)-np.abs(np.diag(a+b)))).min()
                for _, _, a, b in blocks)),
            "tilts": []}
        if candidate == "C":
            integer_population = np.zeros((66,66),dtype=np.int64)
            for destinations in native.destinations:
                columns = np.flatnonzero(destinations >= 0)
                integer_population[destinations[columns],columns] += 1
            x = sp.symbols("x")
            actual_charpoly = sp.Matrix(integer_population.tolist()).charpoly(x).as_expr()
            expected_charpoly = ((x-60)*(x-54)**5*(x-52)**5*(x-50)**11*(x-48)**5
                *(x*x-110*x+3020)**3*(x*x-100*x+2488)
                *(x*x-100*x+2490)**4*(x*x-100*x+2498)**4
                *(x**3-150*x*x+7480*x-124024)**5)
            if sp.expand(actual_charpoly-expected_charpoly) != 0:
                raise ArithmeticError("exact full native population spectrum changed")
            contraction = (55+sp.sqrt(5))/60
            report["exact_native_preparation"] = {
                "characteristic_polynomial_of_60_times_population_transfer": str(expected_charpoly),
                "exact_nonconstant_quantum_channel_spectral_radius": str(contraction),
                "binary_input_trace_distance_bound_after_n_attempts": "sqrt(65)/2 * ((55+sqrt(5))/60)**n",
                "reference_stable_half_diamond_bound_after_n_attempts": "min(1,33*((55+sqrt(5))/60)**n)",
                "stationary_native_preparation_is_limit_of_complete_recorded_channel": True,
                "retained_tape_capacity_separately_required": True,
            }
        for z in (1., .5, .25):
            driven = doob_instrument(native, z)
            inst, state = driven["instrument"], driven["stationary_state"]
            spectra = np.sort(np.concatenate([np.linalg.eigvalsh((a+z*b)/60)
                                             for _, _, a, b in blocks]))
            lam = driven["perron_value"]
            count = count_distribution(inst, state, 5)
            n = np.arange(6)
            filtered = filter_implementation(native, z, rho)
            probabilities = effects_diagonal(inst)@np.diag(state).real
            mean = float(probabilities@inst.costs)
            row = {
                "z": z, "s": -math.log(z), "perron_value": lam,
                "perron_effect_by_distance": driven["perron_effect_by_distance"].tolist(),
                "stationary_orbit_probabilities": [float(np.trace(state[np.ix_(driven['orbit'] == j, driven['orbit'] == j)]).real)
                                                    for j in (1, 2, 3)],
                "stationary_activity_per_attempt": mean,
                "tp_residual": float(max(abs(effects_diagonal(inst).sum(axis=0)-1))),
                "stationarity_residual": float(np.max(abs(apply(inst, state)-state))),
                "all_operator_eigenvalues_count": len(spectra),
                "tilted_transfer_minimum_eigenvalue": float(spectra[0]),
                "top_eigenvalue_residual": abs(float(spectra[-1])-lam),
                "top_eigenvalue_multiplicity_at_2e-12": int(np.count_nonzero(abs(spectra-lam) < 2e-12)),
                "euclidean_gap_per_attempt": float(-math.log(spectra[-2]/lam)),
                "euclidean_gap_is_a_physical_particle_mass": False,
                "five_attempt_count_probabilities": count["probabilities"].tolist(),
                "five_attempt_count_mean_residual": abs(float(n@count["probabilities"])-5*mean),
                "five_attempt_count_variance": float((n*n)@count["probabilities"]-(n@count["probabilities"])**2),
                "complete_filter_record_labels": len(filtered["outcomes"]),
                "complete_filter_total_probability_residual": abs(filtered["total_probability"]-1),
                "herald_success_probability": filtered["success_probability"],
                "herald_success_probability_residual": abs(filtered["success_probability"]-filtered["predicted_success_probability"]),
                "herald_success_output_residual": filtered["successful_output_residual"],
            }
            summary["tilts"].append(row)
            tag = f"{candidate}_z{str(z).replace('.', '_')}"
            archive[tag+"_amplitudes"] = inst.amplitudes
            archive[tag+"_stationary_state"] = state
            archive[tag+"_all_transfer_eigenvalues"] = spectra
            archive[tag+"_count_probabilities"] = count["probabilities"]
        # Diagonal clocks cannot see an environment phase; a declared
        # coherent read before record commitment can, and every outcome stays.
        psi = np.zeros(d, complex)
        a, b = native.seams[0]
        other = next(i for i in range(12) if i not in (a, b))
        source_state, target_state = [0]*12, [0]*12
        source_state[a] = source_state[other] = 1
        target_state[b] = target_state[other] = 1
        psi[native.states.index(tuple(source_state))] = 1/math.sqrt(2)
        psi[native.states.index(tuple(target_state))] = 1/math.sqrt(2)
        phases = np.ones(m, complex)
        raw_move = 1 if candidate == "C" else 5
        phases[native.labels.index((0, 0, raw_move))] = -1
        ordinary = coherent_record_read(native, psi)
        phased = coherent_record_read(native, psi, phases)
        committed_a = coherent_record_read(native, psi, commit_first=True)
        committed_b = coherent_record_read(native, psi, phases, commit_first=True)
        summary["record_phase_audit"] = {
            "complete_committed_law_difference": float(np.max(abs(ordinary["committed_probabilities"]-phased["committed_probabilities"]))),
            "activity_clock_law_difference": float(np.max(abs(ordinary["activity_probabilities"]-phased["activity_probabilities"]))),
            "complete_coherent_pointer_total_variation": float(.5*np.sum(abs(ordinary["pointer_probabilities"]-phased["pointer_probabilities"]))),
            "commit_first_pointer_law_difference": float(np.max(abs(committed_a["pointer_probabilities"]-committed_b["pointer_probabilities"]))),
            "all_pointer_outcomes": m,
            "seam_and_coin_randomizers_remain_committed": True,
        }
        summary["record_conditioned_geometric_reference"] = reference_audit(native, repo, archive)
        archive[candidate+"_destinations"] = native.destinations
        archive[candidate+"_full78_destinations"] = full_native.destinations
        archive[candidate+"_full78_amplitudes"] = full_native.amplitudes
        archive[candidate+"_full78_states"] = np.asarray(full_native.states,dtype=np.int64)
        archive[candidate+"_labels"] = np.asarray(native.labels, dtype=np.int64)
        archive[candidate+"_costs"] = native.costs
        archive[candidate+"_coherent_pointer_original"] = ordinary["pointer_probabilities"]
        archive[candidate+"_coherent_pointer_phased"] = phased["pointer_probabilities"]
        report["candidates"][candidate] = summary

    # Exact symbolic quotient characteristic and long-time clock cumulants.
    u, z = sp.symbols("u z", real=True)
    p = u**3-(2+8*z)*u**2+(8*z-20*z*z)*u+40*z*z+80*z**3
    du = sp.simplify(-sp.diff(p, z)/sp.diff(p, u))
    d2u = sp.simplify(-(sp.diff(p, z, 2)+2*sp.diff(p, u, z)*du+sp.diff(p, u, 2)*du**2)/sp.diff(p, u))
    uz, uzz = du.subs({z:1, u:10}), d2u.subs({z:1, u:10})
    # s=-log z, so theta_s'=-z lambda_z/lambda and
    # theta_s''=z lambda_z/lambda+z²(lambda_zz/lambda-(lambda_z/lambda)²).
    mean = sp.simplify(uz/60)
    variance_rate = sp.simplify((uz+uzz)/60-(uz/60)**2)
    report["exact_activity_quotient"] = {
        "u_definition": "u=60*lambda-50",
        "characteristic_polynomial": str(p),
        "polynomial_at_z_one": str(sp.factor(p.subs(z,1))),
        "native_mean_activity_per_attempt": str(mean),
        "native_asymptotic_activity_variance_per_attempt": str(variance_rate),
        "iid_variance_would_be": str(sp.simplify(mean*(1-mean))),
        "native_clock_is_iid": bool(variance_rate == mean*(1-mean)),
    }
    # Independent exact Poisson-equation calculation of the same variance.
    M0 = sp.diag(52,50,50)/60
    A = sp.Matrix([[4,4,0],[4,4,10],[0,2,0]])/60
    T = M0+A
    pi = sp.Matrix([sp.Rational(5,11),sp.Rational(5,11),sp.Rational(1,11)])
    ones = sp.ones(1,3)
    centered = A*pi-mean*pi
    correction = (sp.eye(3)-T+pi*ones).inv()*centered
    covariance_sum = sp.simplify((ones*A*correction)[0])
    poisson_variance = sp.simplify(mean*(1-mean)+2*covariance_sum)
    if poisson_variance != variance_rate:
        raise ArithmeticError("independent native clock variance derivations disagree")
    report["exact_activity_quotient"].update({
        "positive_lag_activity_covariance_sum": str(covariance_sum),
        "independent_Poisson_equation_variance": str(poisson_variance),
        "first_lag_activity_covariance": str(sp.simplify((ones*A*A*pi)[0]-mean**2)),
    })
    # Outward rational isolation of the three reported Perron roots. The
    # printed float is compared as its exact binary rational, so its error
    # is measured rather than silently treated as the exact algebraic root.
    perron_certificates = []
    for zvalue in (sp.Rational(1),sp.Rational(1,2),sp.Rational(1,4)):
        intervals = sp.Poly(p.subs(z,zvalue),u).intervals(eps=sp.Rational(1,10**60))
        (low,high), multiplicity = intervals[-1]
        if multiplicity != 1:
            raise ArithmeticError("Perron root is not isolated and simple")
        low_lam, high_lam = (low+50)/60,(high+50)/60
        candidate_row = next(row for row in report["candidates"]["C"]["tilts"] if row["z"] == float(zvalue))
        floating_value = sp.Rational(candidate_row["perron_value"])
        error = max(abs(floating_value-low_lam),abs(floating_value-high_lam))
        perron_certificates.append({
            "z_exact":str(zvalue),"lambda_lower_rational":str(low_lam),
            "lambda_upper_rational":str(high_lam),
            "floating_perron_value_absolute_error_upper_rational":str(error),
            "floating_perron_value_absolute_error_upper_decimal":str(sp.N(error,25)),
            "root_interval_width_less_than_1e_minus_60":bool(high_lam-low_lam < sp.Rational(1,10**60)),
            "certification":"exact rational real-root isolation; no empirical input"})
    report["exact_perron_root_certificates"] = perron_certificates
    report["capacity_theorem"] = {
        "attempt_count_symbol": "n", "record_capacity_symbol": "N_rec=log(D_rec)",
        "complete_native_record_support_bounds": "60**n <= R_n <= min(m**n,66*60**n); m=120 (C),180 (E)",
        "necessary_capacity_for_n_retained_attempts": "N_rec >= n*log(60)",
        "simple_fixed_alphabet_sufficient_capacity": "N_rec >= n*log(m)",
        "prefix_law_capacity_independent_if_k_le_n": True,
        "public_P_or_N_enters_native_transition_or_its_spectrum": False,
        "finite_capacity_changes_native_masses_without_additional_feedback": False,
    }
    report["source_status"] = {
        "native_quantum_record_instrument_reconstructed": True,
        "unique_stationary_native_state_on_invariant_binary_shell": "I_66/66",
        "nontrivial_conditioned_state_and_dynamics_derived_given_one_tilt": True,
        "tilt_or_pressure_to_tilt_identification_source_selected": False,
        "all_filter_failures_retained": True,
        "positive_native_euclidean_transfer_and_logarithm_constructed": True,
        "real_time_continuation_is_the_native_repair_process": False,
        "native_interacting_3plus1_field_or_twelve_particle_masses_recovered": False,
    }
    # Tight, explicit controls on the claimed finite calculation.
    for summary in report["candidates"].values():
        if summary["zero_tilt_exact_gershgorin_lower_bound_numerator_over_60"] <= 0:
            raise ArithmeticError("positive transfer logarithm not certified")
        for row in summary["tilts"]:
            if row["top_eigenvalue_multiplicity_at_2e-12"] != 1:
                raise ArithmeticError("numerical simple Perron check failed")
            if max(row[k] for k in ("tp_residual", "stationarity_residual",
                                    "top_eigenvalue_residual", "five_attempt_count_mean_residual",
                                    "complete_filter_total_probability_residual",
                                    "herald_success_probability_residual", "herald_success_output_residual")) > 3e-12:
                raise ArithmeticError("finite native record audit exceeded its tolerance")
    (output / "results.json").write_text(json.dumps(report, indent=2)+"\n")
    np.savez_compressed(output / "complete_native_instruments.npz", **archive)
    return report


if __name__ == "__main__":
    from .certificate import main
    main()
