"""Finite scalar words compiled to the existing twelve-port source controls.

Every primitive is exp(sum_p c_p D_p) on the full six-level processor.
The logical one/pair operations use its existing two/four-dimensional proper
codes and coherent load/export service. Those services, blank preparation,
analogue precision and quantum readout remain specified apparatus.

Tapes are chronological JSON-compatible dictionaries. Global qubit numbers
are little endian: qubit j is bit j of a statevector index. A pair request
with targets [a,b] packs into the fixed code |a,b> -> |2*a+b>. All scalar
phases and zero-angle primitive requests are retained. There are no arbitrary
logical matrices, free basis changes or discarded service-failure outcomes.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import hashlib
import json
import math
from numbers import Integral, Real
from pathlib import Path

import numpy as np
from scipy.linalg import expm, svd as scipy_svd


ROOT = Path(__file__).resolve().parents[3]
RESPONSE = "code/source_selection_model/response.json"
CONTROL_RECEIPT = "code/m1_source_realization/receipt.json"
KINDS = ("X", "Y", "Z", "I", "H", "J")
_ZERO = (Fraction(0), Fraction(0))


def _add(a, b):
    return a[0]+b[0], a[1]+b[1]


def _multiply(a, b):
    return a[0]*b[0]+5*a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def _field(value):
    if (not isinstance(value, list) or len(value) != 2
            or any(type(x) is not str for x in value)):
        raise ValueError("two rational strings for a+b*sqrt(5) required")
    return tuple(Fraction(x) for x in value)


def _number(value):
    return float(value[0])+math.sqrt(5)*float(value[1])


def _finite(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be a finite real number")
    return value


def _bits(size, maximum):
    if size < 2 or size & (size-1) or size > 1 << maximum:
        raise ValueError(f"a power-of-two size from 2 through {1 << maximum} required")
    return size.bit_length()-1


def _real_array(value, ndim, name):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot be masked")
    raw = np.asarray(value)
    if raw.ndim != ndim or raw.dtype.kind not in "iuf" or not np.all(np.isfinite(raw)):
        raise ValueError(f"{name} must be a finite real {ndim}-dimensional array")
    return raw.astype(float, copy=True)


@lru_cache(maxsize=1)
def _source():
    """Check the retained exact control identities before using their pulses.

    Arithmetic here is in Q(sqrt(5)); no inverse problem fits a desired gate.
    The existing independent source checker separately reconstructs the D_p
    from their geometric formula. Tests use that independent route as well.
    """
    response = json.loads((ROOT/RESPONSE).read_text(encoding="utf-8"))
    retained = json.loads((ROOT/CONTROL_RECEIPT).read_text(encoding="utf-8"))
    blocks_list = response["generators"]
    if len(blocks_list) != 12:
        raise ValueError("the source must have exactly twelve generators")
    exact = []
    numerical = []
    for blocks in blocks_list:
        if (len(blocks) != 2 or any(len(block) != 3 or any(len(row) != 3 for row in block)
                                   for block in blocks)):
            raise ValueError("two three-dimensional source blocks required")
        full = [[(_ZERO, _ZERO) for _ in range(6)] for _ in range(6)]
        for b, block in enumerate(blocks):
            for i, row in enumerate(block):
                for j, cell in enumerate(row):
                    if len(cell) != 4:
                        raise ValueError("complex source field entry required")
                    full[3*b+i][3*b+j] = (_field(cell[:2]), _field(cell[2:]))
        exact.append(full)
        numerical.append(np.array([[complex(_number(r), _number(im)) for r, im in row]
                                   for row in full]))
    targets = {}
    for i in range(3):
        target = np.zeros((6, 6), complex)
        target[i, i] = 1j
        targets[f"phase_{i}"] = target
    for i, j in ((0, 1), (0, 2), (1, 2)):
        for kind in ("real", "imaginary"):
            target = np.zeros((6, 6), complex)
            target[i, j], target[j, i] = ((1, -1) if kind == "real" else (1j, 1j))
            targets[f"{kind}_{i}{j}"] = target
    targets["relative_phase"] = np.diag([1j, 1j, 1j, 0, 0, 0])
    rows = retained["evidence"]["controls"]["rows"]
    if [row["name"] for row in rows] != list(targets):
        raise ValueError("source control catalogue changed")
    coefficients = {}
    for row in rows:
        c = [_field(x) for x in row["coefficients"]]
        if len(c) != 12:
            raise ValueError("twelve integrated port coefficients required")
        target = targets[row["name"]]
        for i in range(6):
            for j in range(6):
                parts = []
                for component in (0, 1):
                    result = _ZERO
                    for coefficient, generator in zip(c, exact, strict=True):
                        result = _add(result, _multiply(coefficient, generator[i][j][component]))
                    parts.append(result)
                z = target[i, j]
                if parts != [(Fraction(int(z.real)), Fraction(0)),
                             (Fraction(int(z.imag)), Fraction(0))]:
                    raise ValueError(f"exact source control mismatch: {row['name']}")
        coefficients[row["name"]] = c
    selected = {
        "X": coefficients["imaginary_01"],
        "Y": coefficients["real_01"],
        "Z": [(a[0]-b[0], a[1]-b[1]) for a, b in
              zip(coefficients["phase_0"], coefficients["phase_1"], strict=True)],
        "I": [_add(a, b) for a, b in
              zip(coefficients["phase_0"], coefficients["phase_1"], strict=True)],
        "J": coefficients["relative_phase"],
    }
    numeric_coefficients = {kind: np.array(list(map(_number, values)))
                            for kind, values in selected.items()}
    numeric_coefficients["H"] = (numeric_coefficients["X"]+numeric_coefficients["Z"])/math.sqrt(2)
    generators = np.array(numerical)
    for array in [generators, *numeric_coefficients.values()]:
        array.flags.writeable = False
    return generators, numeric_coefficients, selected


def source_summary():
    """Serializable source commitments and exact base-control coefficients."""
    generators, coefficients, exact = _source()
    return {
        "source_sha256": {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                          for path in (RESPONSE, CONTROL_RECEIPT)},
        "coefficient_field": "Q(sqrt(5)); [a,b] means a+b*sqrt(5)",
        "exact_base_coefficients": {kind: [[str(a), str(b)] for a, b in row]
                                    for kind, row in exact.items()},
        "H_coefficients": "(X+Z)/sqrt(2)",
        "exact_control_rows_checked": 10,
        "source_generators": int(len(generators)),
        "coherent_proper_code_dimensions": [2, 4],
        "coherent_code_gluing_supplied": True,
        "full_algebra_transfer_has_retained_failure_outcome": True,
        "native_generator_antihermitian_residual": float(max(
            np.linalg.norm(g+g.conj().T, 2) for g in generators)),
        "maximum_base_coefficient": float(max(np.max(abs(c)) for c in coefficients.values())),
    }


def integrated_controls(kind, angle):
    """The twelve actual integrated response coefficients for one primitive."""
    if kind not in KINDS:
        raise ValueError("unknown source pulse kind")
    angle = _finite(angle, "pulse angle")
    return angle*_source()[1][kind]


@lru_cache(maxsize=32768)
def _native(kind, angle):
    coefficients = integrated_controls(kind, angle)
    generator = np.tensordot(coefficients, _source()[0], axes=(0, 0))
    full = expm(generator)
    full.flags.writeable = False
    return full


def native_pulse(kind, angle):
    """Execute the full M6 source response, without inserting a logical gate."""
    if kind not in KINDS:
        raise ValueError("unknown source pulse kind")
    return _native(kind, _finite(angle, "pulse angle")).copy()


def transfer_kraus(d, coherent=True):
    """Complete M6 -> code+failure export, including every complement branch.

    The receiver's previous contents are overwritten and the sender is reset
    by the declared service. Those fixed outputs are omitted from these
    rectangular matrices. Valid-code identity preserves arbitrary references;
    the dephasing control supplies one retained Kraus branch per code value.
    """
    if type(d) is not int or d not in (2, 4) or type(coherent) is not bool:
        raise ValueError("a declared d=2 or d=4 code and Boolean coherence required")
    kraus = []
    if coherent:
        success = np.zeros((d+1, 6), complex)
        success[:d, :d] = np.eye(d)
        kraus.append(success)
    else:
        for j in range(d):
            branch = np.zeros((d+1, 6), complex)
            branch[j, j] = 1
            kraus.append(branch)
    for j in range(d, 6):
        branch = np.zeros((d+1, 6), complex)
        branch[d, j] = 1
        kraus.append(branch)
    return kraus


def pulse(kind, angle, targets):
    """Construct one request; zero angle still has a load and export."""
    if kind not in KINDS:
        raise ValueError("unknown source pulse kind")
    expected = 2 if kind == "J" else 1
    if (not isinstance(targets, (list, tuple)) or len(targets) != expected
            or any(isinstance(j, (bool, np.bool_)) or not isinstance(j, Integral) or j < 0
                   for j in targets) or len(set(targets)) != len(targets)):
        raise ValueError("distinct nonnegative target qubits of the native code required")
    return {"kind": kind, "targets": [int(j) for j in targets], "angle": _finite(angle, "pulse angle")}


def _validate(tape, qubits=None):
    if not isinstance(tape, (list, tuple)):
        raise ValueError("a chronological native tape required")
    out = []
    for item in tape:
        if type(item) is not dict or set(item) != {"kind", "targets", "angle"}:
            raise ValueError("native pulse schema mismatch")
        event = pulse(item["kind"], item["angle"], item["targets"])
        if qubits is not None and max(event["targets"]) >= qubits:
            raise ValueError("pulse target outside the requested register")
        out.append(event)
    return out


def inverse(tape):
    return [pulse(event["kind"], -event["angle"], event["targets"])
            for event in reversed(_validate(tape))]


def remap(tape, targets):
    if (not isinstance(targets, (list, tuple)) or len(set(targets)) != len(targets)
            or any(isinstance(j, (bool, np.bool_)) or not isinstance(j, Integral) or j < 0
                   for j in targets)):
        raise ValueError("distinct nonnegative register destinations required")
    return [pulse(event["kind"], event["angle"], [targets[j] for j in event["targets"]])
            for event in _validate(tape, len(targets))]


def cnot(control, target):
    """Phase-exact iH, -CZ, iH word from the existing relative-phase source."""
    return [pulse("H", math.pi/2, [target]), pulse("J", -math.pi, [control, target]),
            pulse("H", math.pi/2, [target])]


def pauli_rotation(label, angle):
    """Compile exp(i*angle*P); label[j] is the Pauli at little-endian qubit j."""
    if not isinstance(label, (str, list, tuple)) or not label or any(p not in "IXYZ" for p in label):
        raise ValueError("a nonempty tensor-Pauli label required")
    angle = _finite(angle, "Pauli angle")
    support = [j for j, p in enumerate(label) if p != "I"]
    if not support:
        return [pulse("I", angle, [0])]
    basis = []
    for j, p in enumerate(label):
        if p == "X":
            basis.append(pulse("Y", -math.pi/4, [j]))
        elif p == "Y":
            basis.append(pulse("X", math.pi/4, [j]))
    pivot = support[-1]
    parity = [event for j in support[:-1] for event in cnot(j, pivot)]
    return inverse(basis)+parity+[pulse("Z", angle, [pivot])]+inverse(parity)+basis


def compile_diagonal(phases):
    """Compile diag(exp(i*phases)) by its complete commuting Walsh expansion."""
    values = _real_array(phases, 1, "diagonal phases")
    bits = _bits(len(values), 8)
    word = []
    for mask in range(len(values)):
        coefficient = math.fsum(float(x)*(-1)**((index & mask).bit_count())
                                for index, x in enumerate(values))/len(values)
        label = ["Z" if mask >> j & 1 else "I" for j in range(bits)]
        word.extend(pauli_rotation(label, coefficient))
    return word


def _gray_rotation(a, b, theta, bits):
    """exp(i theta Y_ab), with ordered Gray-neighbor coordinates (a,b)."""
    difference = a ^ b
    if not difference or difference & (difference-1):
        raise ValueError("a Gray-neighbor pair required")
    target = difference.bit_length()-1
    orientation = (-1)**((a >> target) & 1)
    word = []
    for mask in range(1 << bits):
        if mask & difference:
            continue
        label = ["Y" if j == target else "Z" if mask >> j & 1 else "I" for j in range(bits)]
        coefficient = orientation*theta*(-1)**((a & mask).bit_count())/(1 << (bits-1))
        word.extend(pauli_rotation(label, coefficient))
    return word


def compile_real_orthogonal(matrix):
    """Compile a complete real orthogonal basis on at most four qubits.

    Gray order is applied to both rows and pivot columns; the remaining
    diagonal signs and all fixed-bit projectors are actual emitted gates.
    """
    value = _real_array(matrix, 2, "orthogonal basis")
    bits = _bits(value.shape[0], 4)
    if value.shape != (len(value), len(value)) or np.linalg.norm(value.T@value-np.eye(len(value)), 2) > 2e-11:
        raise ValueError("a real orthogonal square basis required")
    gray = [j ^ (j >> 1) for j in range(len(value))]
    residual = value[np.ix_(gray, gray)].copy()
    rotations = []
    for column in range(len(value)-1):
        for row in range(len(value)-1, column, -1):
            x, y = residual[row-1, column], residual[row, column]
            # The schedule includes empty pivots. Keeping their zero-angle
            # requests avoids changing service counts with floating QR zeros.
            theta = 0. if x == 0 and y == 0 else math.atan2(y, x)
            c, s = math.cos(theta), math.sin(theta)
            old = residual[[row-1, row], :].copy()
            residual[row-1], residual[row] = c*old[0]+s*old[1], -s*old[0]+c*old[1]
            rotations.append((gray[row-1], gray[row], theta))
    signs = np.where(np.diag(residual) >= 0, 1., -1.)
    if np.linalg.norm(residual-np.diag(signs), 2) > 3e-11:
        raise ArithmeticError("Gray elimination failed to retain the complete basis")
    phases = np.zeros(len(value))
    for j, sign in enumerate(signs):
        phases[gray[j]] = 0 if sign == 1 else math.pi
    word = compile_diagonal(phases)
    for a, b, theta in reversed(rotations):
        word.extend(_gray_rotation(a, b, -theta, bits))
    return word


def _compile_real_state_tree(value, bits):
    """Binary amplitude tree with commuting Walsh-Y factors at each level.

    Higher bits label the prepared prefix. A level splits its remaining norm
    between target zero and one; the final level retains signed amplitudes.
    Every conditional rotation is compiled from native tensor-Pauli words,
    including all zero coefficients, in a dimension-dependent fixed schedule.
    """
    word = []
    for target in reversed(range(bits)):
        controls = bits-target-1
        angles = []
        for prefix in range(1 << controls):
            block = value[prefix*(1 << (target+1)):(prefix+1)*(1 << (target+1))]
            half = 1 << target
            x = float(block[0]) if target == 0 else float(np.linalg.norm(block[:half]))
            y = float(block[1]) if target == 0 else float(np.linalg.norm(block[half:]))
            angles.append(-math.atan2(y, x) if x != 0 or y != 0 else 0.)
        for mask in range(1 << controls):
            coefficient = math.fsum(theta*(-1)**((prefix & mask).bit_count())
                                    for prefix, theta in enumerate(angles))/(1 << controls)
            label = ["I"]*bits
            label[target] = "Y"
            for control in range(controls):
                if mask >> control & 1:
                    label[target+1+control] = "Z"
            word.extend(pauli_rotation(label, coefficient))
    return word


def compile_real_state(vector):
    """Prepare a normalized real state from |0>, using at most twelve qubits.

    This emits a state-preparation word, not an arbitrary state replacement.
    Its action on all other inputs is the same explicitly compiled unitary.
    Through six qubits the existing Gray word is preserved. Larger vectors
    use an O(bits*2**bits) binary-tree word with commuting Walsh-Y factors.
    The admitted maximum is 4096 real amplitudes and 135171 emitted pulses.
    """
    value = _real_array(vector, 1, "state amplitudes")
    bits = _bits(len(value), 12)
    if abs(np.linalg.norm(value)-1) > 2e-12:
        raise ValueError("a normalized real state required")
    if bits > 6:
        word = _compile_real_state_tree(value, bits)
        if len(word) != 3*((bits-1)*(1 << bits)+1) or len(word) > 135171:
            raise ArithmeticError("the bounded amplitude-tree source schedule changed")
        return word
    gray = [j ^ (j >> 1) for j in range(len(value))]
    residual = value[gray].copy()
    rotations = []
    for row in range(len(value)-1, 0, -1):
        x, y = residual[row-1], residual[row]
        # Retain empty-subtree rotations and their proper-code services.
        theta = 0. if x == 0 and y == 0 else math.atan2(y, x)
        residual[row-1], residual[row] = math.hypot(x, y), 0
        rotations.append((gray[row-1], gray[row], theta))
    word = []
    for a, b, theta in reversed(rotations):
        word.extend(_gray_rotation(a, b, -theta, bits))
    return word


def _apply(state, gate, targets, qubits):
    axes = [qubits-1-j for j in targets]
    axes += [axis for axis in range(qubits) if axis not in axes]
    axes += [qubits]
    arranged = state.reshape([2]*qubits+[state.shape[1]]).transpose(axes)
    updated = (gate@arranged.reshape(len(gate), -1)).reshape(arranged.shape)
    return updated.transpose(np.argsort(axes)).reshape(state.shape)


def execute(tape, qubits, state=None):
    """Execute all native pulses on every input column, or on supplied data.

    Contraction to the valid proper code is exact for these block-preserving
    generators. The full export channel's complement outcomes are exposed by
    transfer_kraus and tested independently; no state is renormalized here.
    State columns may retain arbitrary untouched reference entanglement.
    """
    if (isinstance(qubits, (bool, np.bool_)) or not isinstance(qubits, Integral)
            or not 1 <= qubits <= 16):
        raise ValueError("an integer register size from one through sixteen required")
    qubits = int(qubits)
    tape = _validate(tape, qubits)
    vector = False
    if state is None:
        if qubits > 10:
            raise ValueError("complete unitary execution is limited to ten qubits")
        out = np.eye(1 << qubits, dtype=complex)
    else:
        if np.ma.isMaskedArray(state):
            raise ValueError("state cannot be masked")
        out = np.asarray(state)
        if (out.ndim not in (1, 2) or out.shape[0] != 1 << qubits
                or out.dtype.kind not in "iufc" or not np.all(np.isfinite(out))):
            raise ValueError("finite statevector or reference-column array on this register required")
        vector = out.ndim == 1
        out = np.array(out, dtype=complex, copy=True).reshape(1 << qubits, -1)
    for event in tape:
        full = _native(event["kind"], event["angle"])
        dimension = 1 << len(event["targets"])
        # The exact source identities imply invariant code. This guard checks
        # the numerical backend; it never deletes a nonzero physical branch.
        if np.linalg.norm(full[dimension:, :dimension], 2) > 2e-11:
            raise ArithmeticError("native pulse leaked out of its declared proper code")
        out = _apply(out, full[:dimension, :dimension], event["targets"], qubits)
    return out[:, 0] if vector else out


def stabilize_executed_unitary(matrix):
    """Remove bounded floating drift from an already executed unitary word.

    The caller supplies the raw contraction of source pulses, never a target
    matrix. If M=L diag(s) R*, its polar unitary is L R*. This is a numerical
    representation of the same ideally unitary subroutine, not a physical
    channel, error-correction service, new primitive, or state normalization.
    The raw executor is unchanged. Its full-phase target checks remain useful
    because projection cannot repair an incorrect unitary gate or phase.

    Fixed mandatory guards bound the input defect, actual matrix change and
    restored defect. They reject a materially nonunitary request. Returned
    norms describe floating computations; they are not interval certificates
    or estimates of physical gate precision. No state or outcome is used.
    """
    if np.ma.isMaskedArray(matrix):
        raise ValueError("executed unitary cannot be masked")
    raw = np.asarray(matrix)
    if (raw.ndim != 2 or raw.shape[0] != raw.shape[1] or raw.dtype.kind not in "iufc"
            or not np.all(np.isfinite(raw))):
        raise ValueError("a finite square executed-unitary matrix required")
    _bits(raw.shape[0], 10)
    raw = np.array(raw, dtype=complex, copy=True)
    identity = np.eye(len(raw))
    input_limit = 1e-10
    correction_limit = 5.1e-11
    output_limit = 2e-13
    defect = float(np.linalg.norm(raw.conj().T@raw-identity, 2))
    if not math.isfinite(defect) or defect > input_limit:
        raise ArithmeticError("executed-unitary defect exceeds the fixed rounding guard")
    try:
        left, _, right = np.linalg.svd(raw, full_matrices=False)
    except np.linalg.LinAlgError:
        # Some LAPACK builds fail to converge in the divide-and-conquer
        # driver on these almost exactly degenerate singular values. Use
        # the independent QR driver on the same executed matrix. The input,
        # correction and output guards are identical for either backend.
        left, _, right = scipy_svd(raw, full_matrices=False,
                                   lapack_driver="gesvd", check_finite=False)
    corrected = left@right
    correction = float(np.linalg.norm(corrected-raw, 2))
    restored = float(np.linalg.norm(corrected.conj().T@corrected-identity, 2))
    if (not math.isfinite(correction) or correction > correction_limit
            or not math.isfinite(restored) or restored > output_limit):
        raise ArithmeticError("executed-unitary polar correction exceeds its fixed rounding guard")
    # In exact arithmetic, ||polar(M)-M|| <= delta/(1+sqrt(1-delta)).
    # This stable expression avoids cancellation for tiny delta. The value
    # is diagnostic, since delta and the SVD are themselves floating results.
    polar_distance_bound = defect/(1+math.sqrt(1-defect))
    return corrected, {
        "method": "guarded polar projection of the executed source-word matrix",
        "dimension": len(raw),
        "raw_unitarity_defect": defect,
        "projection_change": correction,
        "restored_unitarity_defect": restored,
        "exact_polar_distance_formula_evaluated": polar_distance_bound,
        "maximum_raw_unitarity_defect": input_limit,
        "maximum_projection_change": correction_limit,
        "maximum_restored_unitarity_defect": output_limit,
        "physical_error_certificate": False,
    }


def resources(tape):
    """Retained word counts, including services for all zero-angle requests."""
    tape = _validate(tape)
    pair = sum(event["kind"] == "J" for event in tape)
    unique = {(event["kind"], event["angle"]) for event in tape}
    leakage = max((np.linalg.norm(_native(kind, angle)[(4 if kind == "J" else 2):,
                                                           :(4 if kind == "J" else 2)], 2)
                   for kind, angle in unique), default=0.)
    return {
        "native_pulses": len(tape),
        "single_code_pulses": len(tape)-pair,
        "pair_code_pulses": pair,
        "zero_angle_requests": sum(event["angle"] == 0 for event in tape),
        "identity_phase_pulses": sum(event["kind"] == "I" for event in tape),
        "proper_code_load_events": len(tape),
        "proper_code_export_events": len(tape),
        "proper_code_service_events": 2*len(tape),
        "integrated_processor_norm_upper": math.fsum(abs(event["angle"]) for event in tape),
        "maximum_work_qubits": max((max(event["targets"]) for event in tape), default=-1)+1,
        "maximum_numerical_valid_code_leakage": float(leakage),
        "service_failure_outcomes_retained": True,
        "valid_code_failure_probability_in_exact_arithmetic": 0,
        "initialization_storage_readout_physical_costs_included": False,
    }
