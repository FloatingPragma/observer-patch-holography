"""Exact uncertainty observables for the six compensated Higgs probes.

These are finite, even, charge-compensated quadratures. They are not a
position/momentum identification or a finite representation of the CCR.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse
import sympy as sp

from ..interacting_extension._source import ROOT, verify_sources
from ..interacting_extension.higgs_join.native_inputs import SOURCE_NAMES
from ..interacting_extension.higgs_join.interacting_carrier import load_operators
from ..interacting_extension.recorded_response import partial_isometry

HERE = Path(__file__).resolve().parent
HIGGS = ROOT / "model/interacting_extension/higgs_join"


@dataclass(frozen=True)
class Probe:
    """A signed partial permutation in the complete orthonormal carrier."""

    name: str
    neutral_dimension: int
    charged_dimension: int
    entries: tuple[tuple[int, int, int], ...]

    def __post_init__(self):
        if any(type(n) is not int or n <= 0 for n in
               (self.neutral_dimension, self.charged_dimension)):
            raise ValueError("positive complete sector dimensions required")
        rows, columns = set(), set()
        for i, j, sign in self.entries:
            if (type(i) is not int or type(j) is not int or type(sign) is not int
                    or not 0 <= i < self.charged_dimension
                    or not 0 <= j < self.neutral_dimension or sign not in (-1, 1)
                    or i in rows or j in columns):
                raise ValueError("signed partial permutation required")
            rows.add(i)
            columns.add(j)

    @property
    def dimension(self):
        return self.neutral_dimension + self.charged_dimension

    def creation(self):
        rows = [i for i, _, _ in self.entries]
        columns = [j for _, j, _ in self.entries]
        values = [sign for _, _, sign in self.entries]
        return sparse.csr_matrix((values, (rows, columns)),
                                 shape=(self.charged_dimension, self.neutral_dimension),
                                 dtype=complex)

    def exact_operators(self):
        """Full sparse matrices over Q(i), with no eigensolver or rounding."""
        x, y, z, support = {}, {}, {}, {}
        for i, j, sign in self.entries:
            k = self.neutral_dimension + i
            x[j, k] = x[k, j] = sign
            y[j, k], y[k, j] = -sp.I * sign, sp.I * sign
            z[j, j], z[k, k] = 1, -1
            support[j, j] = support[k, k] = 1
        return {name: sp.SparseMatrix(self.dimension, self.dimension, data)
                for name, data in (("X", x), ("Y", y), ("D", z), ("support", support))}

    def projectors(self, axis="X"):
        """The complete -1, 0, +1 projective instrument, including blocking."""
        if axis not in ("X", "Y"):
            raise ValueError("quadrature must be X or Y")
        b = partial_isometry(self.creation(), self.neutral_dimension,
                             self.charged_dimension)
        phase = 1 if axis == "X" else -1j
        a = sparse.bmat([[None, phase * b.getH()],
                         [phase.conjugate() * b, None]], format="csr")
        support = sparse.block_diag((b.getH() @ b, b @ b.getH()), format="csr")
        unit = sparse.eye(self.dimension, format="csr", dtype=complex)
        return {-1: (support - a) / 2, 0: unit - support, 1: (support + a) / 2}

    def outcomes(self, state, axis="X"):
        """Unnormalized poststates; columns may be arbitrary reference amplitudes."""
        state = np.asarray(state, dtype=complex)
        if (state.ndim not in (1, 2) or state.shape[0] != self.dimension
                or not np.isfinite(state).all()):
            raise ValueError("finite complete-carrier vector or purification columns required")
        return {label: effect @ state for label, effect in self.projectors(axis).items()}


def _probe_from_rows(name, rows, neutral_norms, charged_norms):
    entries = []
    for i, j, value in rows:
        if (type(i) is not int or type(j) is not int
                or not 0 <= i < len(charged_norms) or not 0 <= j < len(neutral_norms)):
            raise ValueError("exact probe coordinate outside the complete carrier")
        numerator = F(value)
        if numerator * numerator != charged_norms[i] * neutral_norms[j]:
            raise ValueError("native normalized probe is not a signed unit entry")
        entries.append((i, j, 1 if numerator > 0 else -1))
    return Probe(name, len(neutral_norms), len(charged_norms), tuple(entries))


def source_bindings():
    bindings = verify_sources(ROOT, SOURCE_NAMES)
    upstream = json.loads((HIGGS / "results.json").read_text())["source_sha256"]
    names = ("joined_basis.json", "joined_exact_operators.json", "joined_operator_receipt.json",
             "joined_operators.npz", "interacting_carrier.py", "native_inputs.py", "source_algebra.py")
    for name in names:
        path = HIGGS / name
        relative = path.relative_to(ROOT).as_posix()
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if upstream.get(relative) != digest:
            raise ValueError("Higgs source differs from its checked receipt: " + name)
        bindings[relative] = digest
    paths = [HIGGS / "results.json", ROOT / "model/interacting_extension/_source.py",
             ROOT / "model/interacting_extension/source_manifest.json",
             ROOT / "model/interacting_extension/recorded_response.py",
             ROOT / "Lean/EventAlgebra/Robertson.lean",
             HERE / "uncertainty.py", HERE / "test_uncertainty.py", HERE / "UNCERTAINTY.md"]
    for path in paths:
        bindings[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(sorted(bindings.items()))


def load_probes():
    """Bind the exact probes to the existing complete numeric Higgs producer."""
    source_bindings()
    basis = json.loads((HIGGS / "joined_basis.json").read_text())
    packet = json.loads((HIGGS / "joined_exact_operators.json").read_text())
    if packet["families"] != 3 or [basis[str(q)]["dimension"] for q in (0, 1)] != [438, 1596]:
        raise ValueError("the complete three-family Higgs carrier is required")
    norms = {q: [F(row["raw_norm_squared"]) for row in basis[str(q)]["labels"]]
             for q in (0, 1)}
    if any(n <= 0 for values in norms.values() for n in values):
        raise ValueError("positive exact basis norms required")
    numeric, _ = load_operators(HIGGS)
    probes = []
    for family in range(3):
        for spin in range(2):
            name = f"addition_{family}_{spin}"
            probe = _probe_from_rows(name, packet["operators"][name], norms[0], norms[1])
            difference = probe.creation() - numeric[name]
            difference.eliminate_zeros()
            if difference.nnz:
                raise ArithmeticError("exact probe differs from the executed Higgs operator")
            active = {j for _, j, _ in probe.entries}
            # The stored basis fixes each charged-lepton mode occupation.
            expected = {j for j, label in enumerate(basis["0"]["labels"])
                        if 2 * family + spin not in label["e_at0"]}
            if active != expected:
                raise ArithmeticError("neutral blocking differs from 1-n_j")
            probes.append(probe)
    return tuple(probes)


def exact_audit(probe):
    operators = probe.exact_operators()
    x, y, d, support = (operators[k] for k in ("X", "Y", "D", "support"))
    unit = sp.SparseMatrix.eye(probe.dimension)
    identities = [x.H - x, y.H - y, d.H - d, x*x - support, y*y - support,
                  d*d - support, support*support - support,
                  x*y - y*x - 2*sp.I*d, x*y + y*x]
    for a in (x, y):
        effects = [(support - a)/2, unit - support, (support + a)/2]
        identities.extend([sum(effects, sp.SparseMatrix.zeros(probe.dimension)) - unit])
        identities.extend(e.H - e for e in effects)
        identities.extend(e*e - e for e in effects)
        identities.extend(e*f for i, e in enumerate(effects) for f in effects[i+1:])
    if any(matrix.nnz() for matrix in identities):
        raise ArithmeticError("exact compensated-probe identity failed")
    rank = len(probe.entries)
    return {"name": probe.name, "neutral_dimension": probe.neutral_dimension,
            "charged_dimension": probe.charged_dimension, "rank_each_signed_outcome": rank,
            "rank_blocked_outcome": probe.dimension - 2*rank,
            "neutral_blocked_dimension": probe.neutral_dimension - rank,
            "exact_rational_complex_identities_checked": len(identities),
            "neutral_tracial_probabilities_minus_zero_plus":
                [str(F(rank, 2*probe.neutral_dimension)),
                 str(1-F(rank, probe.neutral_dimension)),
                 str(F(rank, 2*probe.neutral_dimension))]}


def certificate():
    before = source_bindings()
    rows = [exact_audit(probe) for probe in load_probes()]
    if source_bindings() != before:
        raise ValueError("uncertainty source changed during replay")
    return {"schema": "oph.compensated_probe_uncertainty.v1", "source_sha256": before,
            "probes": rows,
            "neutral_state_theorem": {"domain": "all neutral densities, including arbitrary references",
                "active_probability": "w=Tr(rho b^dagger b)=1-<n_j>",
                "means_X_Y": "0", "variances_X_Y": "w",
                "commutator_expectation": "2 i w", "Robertson_product": "Delta X Delta Y=w",
                "probabilities_minus_zero_plus": ["w/2", "1-w", "w/2"]},
            "scope": {"exact_finite_probe_algebra": True, "all_blocked_outcomes_retained": True,
                "physical_position_momentum_identified": False, "Planck_scale_derived": False,
                "new_native_compiler_for_projective_instrument": False,
                "interacting_relativistic_scattering_derived": False}}


if __name__ == "__main__":
    print(json.dumps(certificate(), sort_keys=True, indent=2, allow_nan=False))
