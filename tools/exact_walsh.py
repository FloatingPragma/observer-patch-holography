"""Exact radix-two Walsh arithmetic for the declared native control compiler.

For N=2**k, H_N[i,j]=(-1)**popcount(i&j).  Splitting an index by
its highest bit gives H_(2N)=[[H_N,H_N],[H_N,-H_N]].  Induction
therefore proves both the butterfly algorithm and H_N**2=N I.  Every
operation below is integer or Fraction addition/subtraction: no numerical
approximation, spectrum fit, or computed-result cache is introduced.

The optional source adapter replaces three pinned compiler implementations.
It changes an arithmetic algorithm, so any execution using it MUST be
reported as an arithmetic-adapted replay, never as metadata-only.  All
native gates, zero coefficients, address ordering, exact binary64 conversion
checks and inverse-transform checks are retained.  The existing finite
compiler proof (model/native_effect_synthesis/PROOF.md, Walsh expansion)
provides the control identity; the induction above changes only its cost.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
from pathlib import Path
import types

SYNTHESIS_PATH = "model/native_effect_synthesis/synthesis.py"
SYNTHESIS_SHA256 = "559ea14f06b7d30ab017b2113feeeebc4d90a517bf85257fd3251d098cb4e372"
FUNCTIONS = ("_walsh", "_state_word", "_phase_word")


def transform(values):
    """Return H_N values exactly, in the original binary index ordering.

    A nonempty power-of-two vector of exact int/Fraction scalars is required.
    The empty-list compatibility of the original private `_walsh` helper is
    handled explicitly by its adapter, not added to this mathematical domain.
    """
    rows = list(values)
    length = len(rows)
    if not length or length & (length-1) or any(type(x) not in (int, Fraction) for x in rows):
        raise ValueError("a nonempty power-of-two exact integer/Fraction vector required")
    width = 1
    while width < length:
        for start in range(0, length, 2*width):
            for offset in range(width):
                left, right = rows[start+offset], rows[start+offset+width]
                rows[start+offset], rows[start+offset+width] = left+right, left-right
        width *= 2
    return tuple(rows)


def replacements(source):
    """Create exact replacements; keep source._walsh visible to mutation tests."""
    F, iv = source.F, source.iv

    def _walsh(values):
        length = len(values)
        if length & (length-1) or any(type(v) is not int for v in values):
            raise ValueError('a power-of-two list of integer angle numerators required')
        if not length:
            return ()
        return tuple(iv.binary_exact(value, iv.AQ*length) for value in transform(values))

    def _state_word(weights):
        if len(weights) != 512 or sum(weights, F()) != 1 or min(weights) < 0:
            raise ValueError('all512 exact nonnegative squared amplitudes required')
        word = []; parameters = []
        for bit in reversed(range(9)):
            controls = 8-bit; angles = []
            for prefix in range(1 << controls):
                block = weights[prefix*(1 << (bit+1)):(prefix+1)*(1 << (bit+1))]
                angles.append(-iv.split_angle(sum(block[:1 << bit], F()),
                                               sum(block[1 << bit:], F())))
            coefficients = source._walsh(angles)
            for mask, coefficient in enumerate(coefficients):
                label = ['I']*9; label[bit] = 'Y'
                for control in range(controls):
                    if mask >> control & 1:
                        label[bit+1+control] = 'Z'
                word.extend(source.compiler.pauli_rotation(label, coefficient))
            # Exact inverse transform of the actual returned coefficients.
            # This still rejects a substituted or corrupted forward compiler.
            reconstructed = transform(tuple(F(c) for c in coefficients))
            for prefix, angle in enumerate(angles):
                if reconstructed[prefix] != F(angle, iv.AQ):
                    raise ArithmeticError('a state-tree angle changed during native compilation')
            parameters.extend(angles)
        if len(word) != 12291 or len(parameters) != 511:
            raise ArithmeticError('the complete nine-bit amplitude-tree schedule changed')
        return word

    def _phase_word(phases):
        if len(phases) != 512:
            raise ValueError('the complete work-register phase list is required')
        coefficients = source._walsh(phases); word = []
        for mask, coefficient in enumerate(coefficients):
            label = ['Z' if mask >> bit & 1 else 'I' for bit in range(9)]
            word.extend(source.compiler.pauli_rotation(label, coefficient))
        reconstructed = transform(tuple(F(c) for c in coefficients))
        for index, phase in enumerate(phases):
            if reconstructed[index] != F(phase, iv.AQ):
                raise ArithmeticError('the emitted diagonal word lost a phase')
        if len(word) != 11270:
            raise ArithmeticError('the complete nine-bit diagonal schedule changed')
        return word

    return dict(zip(FUNCTIONS, (_walsh, _state_word, _phase_word), strict=True))


class SynthesisAdapter:
    """Authenticate source bytes and original code before reversible substitution."""
    def __init__(self, source, root):
        self.source = source
        self.path = Path(root).resolve() / SYNTHESIS_PATH
        if (Path(source.__file__).resolve() != self.path
                or hashlib.sha256(self.path.read_bytes()).hexdigest() != SYNTHESIS_SHA256):
            raise ValueError("native synthesis implementation requires renewed review")
        compiled = compile(self.path.read_bytes(), source.__file__, "exec")
        codes = {code.co_name: code for code in compiled.co_consts if isinstance(code, types.CodeType)}
        self.originals = {name: getattr(source, name) for name in FUNCTIONS}
        for name, function in self.originals.items():
            if (type(function) is not types.FunctionType or function.__code__ != codes[name]
                    or function.__globals__ is not source.__dict__):
                raise ValueError("native synthesis function was already substituted: " + name)
        self.adapted = replacements(source)
        self.installed = False

    def install(self):
        if self.installed:
            raise ValueError("arithmetic adapter is already installed")
        for name, function in self.adapted.items():
            setattr(self.source, name, function)
        self.installed = True

    def verify(self):
        if (not self.installed
                or hashlib.sha256(self.path.read_bytes()).hexdigest() != SYNTHESIS_SHA256
                or any(getattr(self.source, name) is not function for name, function in self.adapted.items())):
            raise ValueError("arithmetic implementation changed during replay")

    def restore(self):
        for name, function in self.originals.items():
            setattr(self.source, name, function)
        self.installed = False

    def evidence(self):
        return {"source": SYNTHESIS_PATH, "source_sha256": SYNTHESIS_SHA256,
                "functions": list(FUNCTIONS), "algorithm_changed": True,
                "algorithm": "exact integer/Fraction radix-two Walsh butterfly",
                "proof": "H_(2N)=[[H_N,H_N],[H_N,-H_N]], hence H_N^2=N I by induction",
                "source_file_changed": False, "computed_result_cache": False,
                "native_gate_schedule_or_exact_inverse_checks_removed": False}
