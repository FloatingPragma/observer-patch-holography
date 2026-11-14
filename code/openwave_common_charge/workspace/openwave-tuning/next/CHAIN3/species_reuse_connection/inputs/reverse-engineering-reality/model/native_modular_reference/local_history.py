"""One locally prepared reference, source clocks and continuing native reads.

The two probes act on the SAME mixed carrier.  Exact two-address rotations
fix every other address; they do not silently assume a selected pure input.
Each analyzer is followed by a complete actual C/E coin/current instrument.
The second rotation's sign is controlled by the first visible current.
Native code export has a retained complement abort, with the work erased on
that terminal branch.  Coherent code services remain specified architecture.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..native_interface import construction as interface
from ..compiled_native_interface import construction as compiled
from ..observer.source_scalar import compiler
from ..native_clock import regional_preparation
from . import common_caps as common, compatibility, full_flow as full, local_preparation as local

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name("local_history_results.json")
PAIRS = ((328, 73), (0, 73))
SEAMS = ((1, 7), (1, 8))
PARAMETER = F(1)
ROTATION_GATE = F(1, 10**8)
WORD_GATE = F(1, 10**9)


def _probe(index, sign=1):
    if type(index) is not int or index not in (0, 1) or type(sign) is not int or sign not in (-1, 1):
        raise ValueError("one of two probes and an integer sign required")
    return index, sign


def pair_word(index, sign):
    """Full controlled two-level rotation, not an unconditional XOR mixer."""
    _probe(index, sign)
    a, b = PAIRS[index]
    difference = a ^ b
    fan = [event for bit in range(1, 9) if difference >> bit & 1
           for event in compiler.cnot(0, bit)]
    # Both start addresses are even. The reversible fan sends b to a^1.
    return [*fan, *compiler._gray_rotation(a, a ^ 1, sign*math.pi/4, 9),
            *compiler.inverse(fan)]


@lru_cache(maxsize=4)
def _rotation(index, sign):
    actual = compiler.execute(pair_word(index, sign), 9)
    return actual.dtype.str, actual.tobytes()


def rotation(index, sign):
    _probe(index, sign)
    dtype, data = _rotation(index, sign)
    return np.frombuffer(data, dtype=np.dtype(dtype)).reshape(512, 512)


def ideal_rotation(index, sign):
    _probe(index, sign)
    a, b = PAIRS[index]
    target = np.eye(512, dtype=complex)
    target[np.ix_((a, b), (a, b))] = np.array([[1, sign], [-sign, 1]])/math.sqrt(2)
    return target


def certify_rotation(index, sign, array=None):
    """Rational Frobenius enclosure of every original binary array entry."""
    _probe(index, sign)
    actual = rotation(index, sign) if array is None else array
    if type(actual) is not np.ndarray:
        raise ValueError("an original base array required")
    actual = full.exact._array(actual, (512, 512))
    lo, hi = full.exact._sqrt_interval(F(1, 2))
    middle, radius = (lo+hi)/2, (hi-lo)/2
    a, b = PAIRS[index]
    special = {(a, a): middle, (b, b): middle, (a, b): sign*middle, (b, a): -sign*middle}
    square = F()
    for i, row in enumerate(actual):
        for j, value in enumerate(row):
            re, im = full.exact._component(value)
            wanted = special.get((i, j), F(int(i == j)))
            square += (re-wanted)**2+im*im
    return full.exact._sqrt_interval(square)[1]+2*radius


@lru_cache(maxsize=1)
def ratios():
    values = [F(1)]*512
    for code, weight in zip(full.support()[0], local.reference(), strict=True):
        values[code] = 60*weight
    return tuple(values)


def word(s=PARAMETER):
    return common.diagonal_word(s, ratios())


@lru_cache(maxsize=2)
def _clock(s):
    array = compiler.execute(word(s), 9)
    return array.dtype.str, array.tobytes()


def clock(s=PARAMETER):
    s = full.parameter(s)
    dtype, data = _clock(s)
    return np.frombuffer(data, dtype=np.dtype(dtype)).reshape(512, 512)


@lru_cache(maxsize=1)
def _source_budget():
    for index in (0, 1):
        for sign in (-1, 1):
            full.corner.report._gate(certify_rotation(index, sign), ROTATION_GATE, "selected native pair rotation")
    full.corner.report._gate(common.certify_diagonal(PARAMETER, ratios(), clock()), WORD_GATE,
                             "locally prepared common clock")
    # Four rotations and two clock words, including both feedback choices.
    error, gain = F(), F(1)
    for eps in (ROTATION_GATE, WORD_GATE, ROTATION_GATE)*2:
        error = (1+eps)**2*error+eps*(2+eps)/2
        gain *= (1+eps)**2
    pulses = sum(2*len(pair_word(i, 1))+len(word()) for i in (0, 1))
    programs = [word(), *(pair_word(i, sign) for i in (0, 1) for sign in (-1, 1))]
    return tuple({"complete_history_half_diamond_upper": error, "channel_norm_upper": gain,
            "maximum_source_pulses": pulses, "maximum_proper_code_service_events": 2*pulses,
            "continuing_native_read_attempts": 2, "native_export_checks": 2,
            "source_address_register_qubits": 9,
            "terminal_flag_and_native_read_record_bits_upper": 13,
            "source_word_position_counter_bits": max(len(p) for p in programs).bit_length(),
            "emitted_program_canonical_utf8_bytes": sum(len(full.corner.canonical(p).encode()) for p in programs),
            "host_interpreter_or_physical_coherent_service_workspace_in_budget": False,
            "projective_work_bit_readouts": 0, "selected_configuration_initialization": False}.items())


def source_budget():
    """Fresh certificate data; callers cannot mutate later bounds or costs."""
    return dict(_source_budget())


def _density(density):
    if (type(density) is not np.ndarray or density.dtype.kind not in "fc"
            or density.ndim != 2 or density.shape[0] != density.shape[1]
            or density.shape[0] not in (66, 132) or not np.isfinite(density).all()):
        raise ValueError("finite native density with at most two reference levels required")
    try:
        with np.errstate(over="raise", invalid="raise"):
            hermitian = float(np.linalg.norm(density-density.conj().T, "nuc"))
            eig = np.linalg.eigvalsh((density+density.conj().T)/2)
            trace = float(np.trace(density).real)
    except (FloatingPointError, np.linalg.LinAlgError) as exc:
        raise ValueError("finite density validation required") from exc
    negative = float(-sum(x for x in eig if x < 0))
    if (not all(math.isfinite(x) for x in (hermitian, negative, trace))
            or not np.isfinite(eig).all()
            or hermitian+negative+max(0., trace-1.) > 1e-13 or trace < 0):
        raise ValueError("positive Hermitian subnormalized density required")
    return density


@lru_cache(maxsize=16)
def _stage(index, sign, raw, s):
    if type(raw) is not bool:
        raise ValueError("a Boolean raw-array switch required")
    s = full.parameter(s)
    r = rotation if raw else ideal_rotation
    u = clock(s) if raw else np.diag([complex(float(w))**(-1j*float(s)) for w in ratios()])
    array = r(index, sign)@u@r(index, -sign)
    return np.frombuffer(array.tobytes(), dtype=array.dtype).reshape(512, 512)


def _read(density, index, sign, law, raw, s):
    """Complete export and current branches, retaining the abort reference."""
    ref = density.shape[0]//66
    codes = np.array(full.support()[0])
    t = _stage(index, sign, raw, s)[:, codes]
    outside = np.array([i for i in range(512) if i not in set(codes)])
    native_t = np.kron(t[codes], np.eye(ref))
    after = native_t@density@native_t.conj().T
    escaped = np.kron(t[outside], np.eye(ref))
    abort = (escaped@density@escaped.conj().T).reshape(len(outside), ref, len(outside), ref)
    abort = np.einsum("aras->rs", abort)
    outputs = {}
    for branch in interface.fixed_maps(law, SEAMS[index]):
        k = np.kron(branch.matrix, np.eye(ref))
        outputs[(branch.coin, branch.label, branch.current)] = k@after@k.conj().T
    return outputs, abort


def execute(density, law, *, raw=True, s=PARAMETER, feedback=True):
    """Two complete continuing instruments; no re-preparation between reads.

    Fine C/E labels are retained, including every zero current. On a native
    export abort the work is discarded but any external reference survives.
    """
    density = _density(density)
    if type(law) is not str or law not in ("C", "E") or type(feedback) is not bool:
        raise ValueError("one fixed native law and Boolean feedback switch required")
    if type(raw) is not bool:
        raise ValueError("a Boolean raw-array switch required")
    s = full.parameter(s)
    if raw and s != PARAMETER:
        raise ValueError("raw histories are certified at the retained parameter only")
    first, abort = _read(density, 0, 1, law, raw, s)
    result = {("abort", 0): abort}
    for record, rho in first.items():
        sign = -1 if feedback and record[2] < 0 else 1
        second, abort = _read(rho, 1, sign, law, raw, s)
        result[("abort", 1, *record)] = abort
        for later, final in second.items():
            result[("read", *record, *later)] = final
    if any(not np.isfinite(rho).all() for rho in result.values()):
        raise ArithmeticError("complete continuing output is nonfinite")
    return result


def prepared_history(law, *, raw=True):
    """Finite source preparation from basis input zero, including all flags.

    The arbitrary-input/reference claim comes from the channel theorem, not
    this one population witness. No successful state is normalized here.
    """
    maps = local.population_maps(law=law)
    result = {("preparation", flag): np.diag(matrix[:, 0]).astype(complex)
              for flag, matrix in maps.items() if flag != "success"}
    result.update(execute(np.diag(maps["success"][:, 0]).astype(complex), law, raw=raw))
    return result


def target_history(law):
    """Complete ideal reference target, with exactly the same flag alphabet."""
    sigma = np.diag(list(map(float, local.reference()))).astype(complex)
    raw_reference = regional_preparation.reference(common.CORE)[0]
    failure = np.diag([float(F(1, 66)-F(5, 33)*w) for w in raw_reference]).astype(complex)
    result = {("preparation", "herald_failure"): failure,
              ("preparation", "warmup_timeout"): np.zeros((66, 66), complex),
              ("preparation", "completion_timeout"): np.zeros((66, 66), complex)}
    result.update(execute(float(F(5, 33))*sigma, law, raw=False))
    return result


def distance(left, right):
    if type(left) is not dict or type(right) is not dict or not left or set(left) != set(right):
        raise ValueError("identical complete record alphabets required")
    result = 0.
    for key in left:
        a, b = left[key], right[key]
        if (any(type(x) is not np.ndarray or x.dtype.kind not in "fc" or x.ndim != 2
                or x.shape[0] != x.shape[1] or not x.size or not np.isfinite(x).all() for x in (a, b))
                or a.shape != b.shape):
            raise ValueError("matching original finite quantum output blocks required")
        try:
            with np.errstate(over="raise", invalid="raise"):
                result += float(np.linalg.norm(a-b, "nuc"))/2
        except (FloatingPointError, np.linalg.LinAlgError) as exc:
            raise ValueError("finite complete quantum distance required") from exc
        if not math.isfinite(result):
            raise ValueError("finite complete quantum distance required")
    return result


def probabilities(outputs):
    return {key: float(np.trace(value).real) for key, value in outputs.items()}


def exact_scope():
    states = common.regional.source()[0].states
    codes = full.support()[0]
    sigma = local.reference()
    for index, (a, b) in enumerate(PAIRS):
        i, j = codes.index(a), codes.index(b)
        if sigma[j]/sigma[i] != F(25, 8):
            raise ArithmeticError("source-selected probe frequency changed")
        if any(states[i][p] != states[j][p] for p in range(12) if p not in common.CAPS[index]):
            raise ArithmeticError("selected probe changes its cap exterior")
    # Pair selection additionally fixes an exterior configuration. It is a
    # controlled global operation, not a bare cap algebra element.
    return {"same_reference_for_both_probes": True,
            "gap_log_argument": "25/8", "same_continuing_carrier": True,
            "source_control_requires_exterior_address_selection": True,
            "strictly_cap_local_coherent_control_derived": False,
            "coherent_load_export_services_supplied": True,
            "physical_time_or_atomic_identification": False}


def _compact_outputs(outputs):
    return [{"record": list(key), "probability": value}
            for key, value in probabilities(outputs).items()]


def run():
    budget = source_budget()
    sigma = np.diag(list(map(float, local.reference()))).astype(complex)
    rows = []
    preparation_budget = local.budget()
    total_bound = (budget["channel_norm_upper"]*F(preparation_budget["complete_flagged_reference_halftrace_upper"])
                   +F(5, 33)*budget["complete_history_half_diamond_upper"])
    for law in ("C", "E"):
        ideal, raw = execute(sigma, law, raw=False), execute(sigma, law)
        error = distance(raw, ideal)
        zero = execute(sigma, law, raw=False, s=F())
        control = distance(ideal, zero)
        if error > float(budget["complete_history_half_diamond_upper"]) or control <= 1e-4:
            raise ArithmeticError("source history fails its bound or nontrivial clock control")
        prepared = prepared_history(law)
        complete_error = distance(prepared, target_history(law))
        if complete_error > float(total_bound):
            raise ArithmeticError("finite preparation and complete history exceed the joined bound")
        rows.append({"law": law, "complete_output": _compact_outputs(prepared),
                     "complete_preparation_and_history_halftrace_diagnostic": complete_error,
                     "complete_preparation_and_history_mass_diagnostic": sum(probabilities(prepared).values()),
                     "raw_total_mass_diagnostic": sum(probabilities(raw).values()),
                     "complete_quantum_halftrace_error_diagnostic": error,
                     "zero_clock_complete_quantum_distance_diagnostic": control})
    paths = sorted(set((*common.SOURCE_PATHS, *interface.SOURCE_PATHS, *local.SOURCE_PATHS,
        "model/native_modular_reference/local_history.py",
        "model/native_modular_reference/local_preparation.py",
        "model/native_modular_reference/compatibility.py",
        "model/native_modular_reference/common_caps.py",
        "model/native_modular_reference/full_flow.py")))
    return {"schema": "oph.native_modular_reference.local_history.v1", "all_pass": True,
            "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
            "source_words": {"clock": hashlib.sha256(full.corner.canonical(word()).encode()).hexdigest(),
                "rotations": [hashlib.sha256(full.corner.canonical(pair_word(i, sign)).encode()).hexdigest()
                              for i in (0, 1) for sign in (-1, 1)]},
            "parameter": str(PARAMETER), "source_budget": full.preparation._compact(budget),
            "all_reference_compatibility": compatibility.exact_certificate(),
            "complete_preparation_and_history_halftrace_upper": full.compiled._outward(total_bound),
            "source_clock_vs_actual_prepared_modular_halfdiamond_upper": full.compiled._outward(
                WORD_GATE*(2+WORD_GATE)/2+abs(PARAMETER)*F(
                    preparation_budget["prepared_modular_generator_operator_error_upper"])),
            "reference_preparation_budget": preparation_budget, "scope": exact_scope(), "histories": rows}


def compare_receipt(supplied, actual):
    full.corner.compare_receipt(supplied, actual)
    for receipt in (supplied, actual):
        for row in receipt["histories"]:
            values = [record["probability"] for record in row["complete_output"]]
            if (any(type(p) not in (int, float) or not math.isfinite(p) or p < 0 or p > 1+2e-10
                    for p in values) or abs(sum(values)-1) > 2e-10
                    or abs(sum(values)-row["complete_preparation_and_history_mass_diagnostic"]) > 2e-12):
                raise ValueError("complete nonnegative probability law required")
            for key, value in row.items():
                if key.endswith("_diagnostic") and (type(value) not in (int, float)
                        or not math.isfinite(value) or value < 0):
                    raise ValueError("nonnegative finite history diagnostic required")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path); group.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = strict_json.loads(full.preparation.canonical(run()))
    if args.check:
        supplied = strict_json.load(args.check)
        compare_receipt(supplied, result)
        print("Locally prepared common-clock continuing history verified.")
    elif args.write:
        with args.write.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(full.corner.canonical(result)+"\n")
    else:
        print(full.corner.canonical(result))


if __name__ == "__main__":
    main()
