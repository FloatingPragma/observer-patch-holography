"""Validation for supplied symmetric neutrino response kernels.

The family-response and holonomy producers normalize a supplied symmetric
kernel ``K`` as ``D^-1/2 K D^-1/2`` with ``D = diag(K)``.  That step is
defined only for a finite symmetric matrix with a strictly positive diagonal,
and the normalized matrix must remain positive semidefinite.  The ``-log``
defect weights are defined only for strictly positive normalized entries.

These checks reject such inputs.  They never substitute a floor for a
nonpositive diagonal, a nonpositive normalized entry, or a non-PSD kernel.
"""

from __future__ import annotations

import numpy as np


# Relative spectral allowance for a supplied kernel that is positive
# semidefinite up to binary64 roundoff.
_PSD_RELATIVE_ALLOWANCE = 1.0e-12


def validated_symmetric_kernel(value, name: str = "K_core_majorana_sym") -> np.ndarray:
    """Return a finite, exactly symmetric, PSD kernel with a positive diagonal."""
    kernel = np.asarray(value, dtype=float)
    if kernel.shape != (3, 3):
        raise ValueError(f"{name} must be a 3x3 matrix")
    if not np.all(np.isfinite(kernel)):
        raise ValueError(f"{name} must be finite")
    if not np.array_equal(kernel, kernel.T):
        raise ValueError(f"{name} must be exactly symmetric")
    if np.any(np.diagonal(kernel) <= 0.0):
        raise ValueError(f"{name} must have a strictly positive diagonal")
    scale = float(np.max(np.abs(kernel)))
    if float(np.min(np.linalg.eigvalsh(kernel))) < -_PSD_RELATIVE_ALLOWANCE * scale:
        raise ValueError(f"{name} must be positive semidefinite")
    return kernel


def normalized_response(kernel: np.ndarray, name: str = "E_nu") -> np.ndarray:
    """Apply the symmetric normalization ``D^-1/2 K D^-1/2`` without flooring."""
    diagonal = np.sqrt(np.diagonal(kernel))
    inverse = np.diag(1.0 / diagonal)
    normalized = inverse @ kernel @ inverse
    if not np.all(np.isfinite(normalized)):
        raise ValueError(f"{name} is outside binary64 range after symmetric normalization")
    return normalized


def defect_weights(normalized: np.ndarray, name: str = "E_nu") -> np.ndarray:
    """Return ``-log(E_nu)``, refusing entries whose logarithm is undefined."""
    if np.any(normalized <= 0.0):
        raise ValueError(f"{name} needs strictly positive entries for -log defect weights")
    return -np.log(normalized)
