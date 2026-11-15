"""Operational regions of the declared finite common-source experiment.

The tensor localization H = H_scalar1 tensor H_scalar0 tensor H_gauge is an
apparatus assumption. Scalar sites carry commuting full matrix algebras. The
four-dimensional gauge factor is the ENTIRE Gauss-constrained charged/link
sector: its two occupation bits are species coordinates, not two independent
spatial subsystems. Neutral scalar operations preserve that physical sector.

Normalized partial trace supplies completely positive conditional expectations
onto regional algebras. It is not a multiplicative restriction of the legacy
FiniteCausalObserverNet contract. The actual Hamiltonian determines which of
these supplied regions interact, via its tensor-centered components and
double commutators. This yields a two-site operational interpretation, without
a metric, three spatial dimensions, physical clock or continuum assertion.
"""
from __future__ import annotations

from itertools import combinations
from numbers import Integral

import numpy as np

from . import scalar_coupling as scalar


REGIONS = ("scalar_site0", "scalar_site1", "gauge_matter_link")
FACTOR_ORDER = ("scalar_site1", "scalar_site0", "gauge_matter_link")


def _region(region):
    if isinstance(region, str):
        region = (region,)
    if not isinstance(region, (tuple, list)) or any(name not in FACTOR_ORDER for name in region):
        raise ValueError("a region is a tuple of declared factor names")
    if len(set(region)) != len(region):
        raise ValueError("a regional factor cannot occur twice")
    return tuple(name for name in FACTOR_ORDER if name in region)


def _matrix(value, size, name="operator"):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} must not be masked")
    raw = np.asarray(value)
    if raw.shape != (size, size) or raw.dtype.kind not in "fiuc":
        raise ValueError(f"{name} must be a numeric {size} by {size} matrix")
    result = np.asarray(raw, dtype=complex)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result


def _layout(d, region):
    d = scalar._dimension(d)
    region = _region(region)
    dims = (d, d, 4)
    kept = tuple(i for i, name in enumerate(FACTOR_ORDER) if name in region)
    rest = tuple(i for i in range(3) if i not in kept)
    return dims, kept, rest


def _size(dims, axes):
    return int(np.prod([dims[i] for i in axes], dtype=int))


def embed_operator(d, region, operator):
    """Embed the region's algebra by tensoring identities on its complement.

    A multi-factor regional matrix uses the order inherited from FACTOR_ORDER,
    independently of the order in which the caller lists its factor names.
    """
    dims, kept, rest = _layout(d, region)
    operator = _matrix(operator, _size(dims, kept), "regional operator")
    order = kept + rest
    joint = np.kron(operator, np.eye(_size(dims, rest)))
    tensor = joint.reshape(tuple(dims[i] for i in order) * 2)
    back = tuple(order.index(i) for i in range(3))
    return tensor.transpose(back + tuple(i + 3 for i in back)).reshape(4*d*d, 4*d*d)


def scalar_operator(d, site, operator):
    """Neutral scalar-site observable; it acts trivially on charged matter/link."""
    if isinstance(site, (bool, np.bool_)) or not isinstance(site, Integral) or site not in (0, 1):
        raise ValueError("scalar site must be 0 or 1")
    return embed_operator(d, f"scalar_site{site}", operator)


def restriction_kraus(d, region):
    """Rectangular Kraus operators for X -> Tr_complement(X)/dim(complement).

    Sum K_b K_b* = I_region proves unitality; the Kraus form proves complete
    positivity, including a spectator reference. This observable restriction
    is distinct from a density-matrix marginal, which has no normalization
    denominator. Embedding its result gives the full-algebra conditional
    expectation used below.
    """
    dims, kept, rest = _layout(d, region)
    keep_size, rest_size = _size(dims, kept), _size(dims, rest)
    indices = np.arange(4*d*d).reshape(dims).transpose(kept + rest).reshape(keep_size, rest_size)
    maps = []
    for column in range(rest_size):
        matrix = np.zeros((keep_size, 4*d*d), complex)
        matrix[np.arange(keep_size), indices[:, column]] = 1 / np.sqrt(rest_size)
        maps.append(matrix)
    return tuple(maps)


def _partial_trace(d, matrix, region):
    dims, kept, rest = _layout(d, region)
    matrix = _matrix(matrix, 4*d*d)
    keep_size, rest_size = _size(dims, kept), _size(dims, rest)
    order = kept + rest
    tensor = matrix.reshape(dims + dims).transpose(order + tuple(i + 3 for i in order))
    return np.trace(tensor.reshape(keep_size, rest_size, keep_size, rest_size), axis1=1, axis2=3)


def conditional_expectation(d, operator, region):
    """Trace-preserving unital CP projection onto a supplied regional algebra.

    E_R(X) = embed_R(Tr_Rc X / dim(Rc)). It fixes regional observables and
    satisfies E_R(A X B)=A E_R(X) B for A,B in that algebra. In general
    E_R(XY) != E_R(X)E_R(Y); it must not fill a *-homomorphism field.
    """
    dims, _, rest = _layout(d, region)
    return embed_operator(d, region, _partial_trace(d, operator, region) / _size(dims, rest))


def regional_density(d, state_or_density, region):
    """Physical partial trace of a state or an unnormalized retained branch.

    No branch probability is divided out. The returned trace is the input
    weight and Tr(rho_R A)=Tr(rho embed_R(A)). Restricting the state does not
    assert autonomous regional evolution under the interacting Hamiltonian.
    """
    _layout(d, region)
    if np.ma.isMaskedArray(state_or_density):
        raise ValueError("a state must not be masked")
    raw = np.asarray(state_or_density)
    if raw.dtype.kind not in "fiuc" or not np.all(np.isfinite(raw)):
        raise ValueError("a finite numeric state is required")
    if raw.shape == (4*d*d,):
        density = np.outer(raw, raw.conj())
    else:
        density = _matrix(raw, 4*d*d, "density")
    if np.linalg.norm(density-density.conj().T, 2) > 1e-11:
        raise ValueError("density must be Hermitian")
    if np.linalg.eigvalsh(density)[0] < -1e-11:
        raise ValueError("density must be positive")
    weight = float(np.trace(density).real)
    if not -1e-11 <= weight <= 1+1e-11:
        raise ValueError("a normalized state or subnormalized branch is required")
    return _partial_trace(d, density, region)


def interaction_components(d, hamiltonian):
    """Unique tensor-centered decomposition relative to the supplied factors.

    For S subset factors, H_S=sum_(T subset S) (-1)^|S-T| E_T(H).
    Sum_S H_S=H, and each nonempty component has zero conditional expectation
    after any one of its factors is removed. Thus nonzero pair components
    are actual cross-region terms, not merely local terms written on a pair.
    """
    scalar._dimension(d)
    hamiltonian = _matrix(hamiltonian, 4*d*d, "Hamiltonian")
    if np.linalg.norm(hamiltonian-hamiltonian.conj().T, 2) > 1e-11:
        raise ValueError("Hamiltonian must be Hermitian")
    subsets = [group for size in range(4) for group in combinations(REGIONS, size)]
    expectations = {group: conditional_expectation(d, hamiltonian, group) for group in subsets}
    return {group: sum(((-1)**(len(group)-len(part)) * expectations[part]
                       for part in subsets if set(part).issubset(group)),
                      np.zeros_like(hamiltonian)) for group in subsets}


def commutator(left, right):
    return left @ right - right @ left


def interaction_graph(d, coupling=scalar.COUPLING, hamiltonian=None):
    """Compute interaction support and witnesses; distances are not assigned.

    The default is the same Hamiltonian executed by scalar_coupling. A supplied
    alternate Hamiltonian is useful for deletion and basis-change controls.
    Graph support uses all centered components, not one possibly vanishing
    probe. The named double commutators are additional witnesses for the
    default construction. Classification tolerances are floating diagnostics.
    """
    scalar._dimension(d)
    hamiltonian = scalar.ideal_hamiltonian(d, coupling) if hamiltonian is None else hamiltonian
    pieces = interaction_components(d, hamiltonian)
    norms = {group: float(np.linalg.norm(value, 2)) for group, value in pieces.items()}
    tolerance = 1e-11
    edges = [list(pair) for pair in combinations(REGIONS, 2)
             if any(norm > tolerance and set(pair).issubset(group) for group, norm in norms.items())]
    number = np.diag(np.arange(d))
    probes = {"scalar_site0": scalar_operator(d, 0, number),
              "scalar_site1": scalar_operator(d, 1, number),
              "gauge_matter_link": embed_operator(d, "gauge_matter_link",
                  np.kron(np.eye(2), np.array([[0., 1.], [1., 0.]])))}
    witnesses = [{"regions": list(pair),
                  "norm": float(np.linalg.norm(commutator(
                      commutator(hamiltonian, probes[pair[0]]), probes[pair[1]]), 2))}
                 for pair in combinations(REGIONS, 2)]
    return {"nodes": list(REGIONS), "edges": edges,
            "centered_component_norms": [{"region": list(group), "norm": norm}
                                          for group, norm in norms.items()],
            "double_commutator_witnesses": witnesses,
            "support_tolerance": tolerance,
            "factorization_and_localization_supplied": True,
            "charged_endpoint_tensor_factorization_claimed": False,
            "metric_or_spatial_dimension_recovered": False}


def occupation_embedding(coarse=2, fine=4):
    """Corner isometry retaining occupation labels, not a unital algebra map."""
    coarse, fine = scalar._dimension(coarse), scalar._dimension(fine)
    if coarse > fine:
        raise ValueError("occupation embedding must increase the cutoff")
    result = np.zeros((4*fine*fine, 4*coarse*coarse))
    for site1 in range(coarse):
        for site0 in range(coarse):
            for gauge in range(4):
                result[gauge+4*(site0+fine*site1), gauge+4*(site0+coarse*site1)] = 1
    return result


def spatial_receipt(d, steps=16, time=1., coupling=scalar.COUPLING):
    """Inspect the source-executed common state with the supplied region maps."""
    state, metadata = scalar.source_evolve(d, steps, time, coupling)
    regional = []
    for name in REGIONS:
        reduced = regional_density(d, state, name)
        regional.append({"region": name, "matrix_dimension": len(reduced),
                         "weight": float(np.trace(reduced).real),
                         "purity": float(np.trace(reduced @ reduced).real),
                         "density": [[[float(z.real), float(z.imag)] for z in row] for row in reduced]})
    return {"dimension": d, "factor_order": list(FACTOR_ORDER),
            "model_time": float(time), "steps": steps,
            "regional_states": regional, "graph": interaction_graph(d, coupling),
            "source_execution": metadata,
            "restriction": "unital CP normalized partial trace; not multiplicative",
            "gauge_localization": "one complete shared matter/link sector; no charged endpoint tensor split",
            "physical_clock_or_spatial_continuum_claimed": False}
