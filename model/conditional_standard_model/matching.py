"""Exact finite checks for a supplied source-to-field sector embedding.

The residual controls spectral localization in the supplied source generator.
It supplies no topological classification, continuum limit or exact source
level multiplicities. The generators must use the same energy reference
and clock. See PROOF.md for the spectral-transfer argument.
"""
from __future__ import annotations

from collections.abc import Mapping

import sympy as s


def _matrix(value, name):
    matrix = s.ImmutableMatrix(value)
    if not matrix.rows or not matrix.cols:
        raise ValueError(f"{name} must be nonempty")
    if any(x.free_symbols or x.has(s.Float) or x.is_finite is not True for x in matrix):
        raise ValueError(f"{name} requires finite exact, fully specified entries")
    return matrix


def _zero(matrix):
    return all(s.simplify(x) == 0 for x in matrix)


def _hermitian(value, name):
    matrix = _matrix(value, name)
    if matrix.rows != matrix.cols or not _zero(matrix - matrix.H):
        raise ValueError(f"{name} must be Hermitian")
    return matrix


def check_intertwiner(source_generator, field_generator, embedding, *, charges=None):
    """Check the entire rectangular residual, including leakage from the image.

    ``charges`` maps labels to (source charge, field charge). These are
    conserved Hermitian observables. Their intertwining must be exact.
    The returned Frobenius residual squared is an exact upper bound on
    the squared operator error for every normalized field state.
    """
    hs = _hermitian(source_generator, "source generator")
    he = _hermitian(field_generator, "field generator")
    j = _matrix(embedding, "embedding")
    if j.shape != (hs.rows, he.rows):
        raise ValueError("embedding dimensions must match both generators")
    if not _zero(j.H*j - s.eye(he.rows)):
        raise ValueError("embedding must be an exact isometry")
    if charges is None:
        charges = {}
    if not isinstance(charges, Mapping):
        raise ValueError("charges must map labels to source/field observables")
    for name, pair in charges.items():
        if not isinstance(name, str) or not name or len(pair) != 2:
            raise ValueError("each charge needs a nonempty label and two observables")
        qs = _hermitian(pair[0], f"source charge {name}")
        qe = _hermitian(pair[1], f"field charge {name}")
        if qs.shape != hs.shape or qe.shape != he.shape:
            raise ValueError("charge dimensions must match their generators")
        if not _zero(qs*hs - hs*qs) or not _zero(qe*he - he*qe):
            raise ValueError(f"{name} is not a conserved charge")
        if not _zero(qs*j - j*qe):
            raise ValueError(f"{name} does not intertwine")
    residual = hs*j - j*he
    error_squared = s.simplify(s.trace(residual.H*residual))
    if error_squared.is_nonnegative is not True:
        raise ValueError("could not verify a nonnegative exact residual norm")
    return {
        "source_dimension": hs.rows,
        "field_dimension": he.rows,
        "isometry_exact": True,
        "conserved_charges_intertwined": sorted(charges),
        "generator_intertwining_exact": _zero(residual),
        "generator_residual_frobenius_squared": str(error_squared),
        "spectral_distance_upper": str(s.sqrt(error_squared)),
        "source_spectrum_exhausted": False,
        "topological_sector_origin_proved": False,
        "continuum_limit_proved": False,
    }


def example_certificate():
    """A rational control distinguishes compression from full intertwining."""
    j = s.ImmutableMatrix([[1, 0], [0, 1], [0, 0]])
    he = s.diag(1, 2)
    exact = s.diag(1, 2, 5)
    leaking = s.ImmutableMatrix([[1, 0, s.Rational(1, 10)],
                                 [0, 2, 0], [s.Rational(1, 10), 0, 5]])
    charges = {"control_charge": (s.diag(1, -1, 1), s.diag(1, -1))}
    return {
        "role": "finite spectral-transfer control, not an OPH defect witness",
        "exact": check_intertwiner(exact, he, j, charges=charges),
        "leakage_control": check_intertwiner(leaking, he, j, charges=charges),
        "leakage_control_compression_equals_field_generator": _zero(j.H*leaking*j - he),
    }
