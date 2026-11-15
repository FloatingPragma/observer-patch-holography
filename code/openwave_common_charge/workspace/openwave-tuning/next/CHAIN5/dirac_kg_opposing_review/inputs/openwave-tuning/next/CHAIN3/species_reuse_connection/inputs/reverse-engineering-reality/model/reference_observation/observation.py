"""Normalized reference transpose maps, without a retraction assertion.

The finite proof is in README.md. These bounded floating matrix controls
neither certify spectral intervals nor supply a native implementation of
the reference square roots. No state or outcome is postselected.
"""
from itertools import combinations
import numpy as np
from ..field_regions.regions import _groups, _modes, regional_basis


WHITENING_TOLERANCE = 1e-10


def matrix(value, shape):
    if np.ma.isMaskedArray(value):
        raise ValueError('masked matrices are not admitted')
    a = np.asarray(value)
    if a.shape != shape or a.dtype.kind not in 'ficu' or not np.isfinite(a).all():
        raise ValueError('finite matrix on the complete declared domain required')
    return np.asarray(a, complex)


def frozen(value):
    return np.frombuffer(value.tobytes(), dtype=value.dtype).reshape(value.shape)


def powers(density):
    """Principal powers of a declared faithful state; no eigenvalue clipping."""
    a = np.asarray(density)
    if a.ndim != 2 or a.shape[0] != a.shape[1] or not 1 <= a.shape[0] <= 512:
        raise ValueError('bounded square density of dimension 1..512 required')
    a = matrix(density, a.shape)
    if np.linalg.norm(a-a.conj().T, 2) > 1e-12 or abs(np.trace(a)-1) > 1e-12:
        raise ValueError('Hermitian unit-trace reference required')
    eigenvalues, vectors = np.linalg.eigh(a)
    if eigenvalues[0] <= 0:
        raise ValueError('reference must be faithful on the declared support')
    return ((vectors*np.sqrt(eigenvalues))@vectors.conj().T,
            (vectors/np.sqrt(eigenvalues))@vectors.conj().T)


class ReferenceObservation:
    """UCP map from source observables to restricted observables.

    A declared CPTP state restriction T with Kraus matrices K_a and faithful
    rho defines E(X)=T(rho^.5 X rho^.5) in the rho_T^(-.5) sandwich.
    Its adjoint is the reference Petz recovery, not ordinary state tracing.
    """
    def __init__(self, reference, restriction):
        root, _ = powers(reference)
        size = root.shape[0]
        if not isinstance(restriction, (tuple, list)) or not restriction:
            raise ValueError('nonempty finite Kraus family required')
        first = np.asarray(restriction[0])
        if first.ndim != 2 or first.shape[1] != size or not 1 <= first.shape[0] <= 512:
            raise ValueError('bounded Kraus output dimension required')
        shape = first.shape
        kraus = tuple(matrix(k, shape) for k in restriction)
        if np.linalg.norm(sum(k.conj().T@k for k in kraus)-np.eye(size), 2) > 1e-12:
            raise ValueError('complete trace-preserving restriction required')
        rho = matrix(reference, (size, size))
        reduced = sum(k@rho@k.conj().T for k in kraus)
        _, inverse = powers(reduced)
        whitened = tuple(inverse@k@root for k in kraus)
        if not all(np.isfinite(k).all() for k in whitened):
            raise ArithmeticError('reference whitening produced nonfinite operators')
        unital_error = np.linalg.norm(
            sum(k@k.conj().T for k in whitened)-np.eye(shape[0]), 2)
        reference_error = np.linalg.norm(
            sum(k.conj().T@reduced@k for k in whitened)-rho, 2)
        if max(unital_error, reference_error) > WHITENING_TOLERANCE:
            raise ArithmeticError('reference whitening exceeds the floating residual budget')
        self.source_reference = frozen(rho)
        self.target_reference = frozen(reduced)
        self.kraus = tuple(frozen(k) for k in whitened)
        self.source_dimension, self.target_dimension = size, shape[0]

    def observe(self, observable, reference_dimension=1):
        if type(reference_dimension) is not int or reference_dimension < 1:
            raise ValueError('positive integer untouched-reference dimension required')
        r = reference_dimension
        x = matrix(observable, (self.source_dimension*r,)*2)
        operators = (np.kron(k, np.eye(r)) for k in self.kraus)
        return sum(k@x@k.conj().T for k in operators)

    def recover(self, density, reference_dimension=1):
        """Trace-preserving adjoint; exact recovery is asserted only for rho_T."""
        if type(reference_dimension) is not int or reference_dimension < 1:
            raise ValueError('positive integer untouched-reference dimension required')
        r = reference_dimension
        y = matrix(density, (self.target_dimension*r,)*2)
        operators = (np.kron(k, np.eye(r)) for k in self.kraus)
        return sum(k.conj().T@y@k for k in operators)


def car_restriction(modes, kept_modes, *, regional_input=False):
    """Existing signed CAR trace as explicit Kraus maps on N=2 or N<=2.

    This bounded diagnostic retains every regional number sector. Choosing
    an output with zero-support sectors requires an explicit support map;
    ReferenceObservation will not silently delete them to force faithfulness.
    """
    modes = _modes(modes)
    kept = _modes(kept_modes, ambient=modes)
    if type(regional_input) is not bool:
        raise ValueError('explicit regional-input Boolean required')
    if not 2 <= len(modes) <= 12:
        raise ValueError('bounded CAR restriction requires 2..12 modes')
    basis = regional_basis(modes) if regional_input else tuple(combinations(modes, 2))
    size = len(regional_basis(kept))
    out = []
    for group in _groups(basis, kept).values():
        k = np.zeros((size, len(basis)), complex)
        for source, local, sign in group:
            k[local, source] = sign
        out.append(frozen(k))
    return tuple(out)
