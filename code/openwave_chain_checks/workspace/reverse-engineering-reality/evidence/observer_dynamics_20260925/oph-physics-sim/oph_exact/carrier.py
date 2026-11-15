"""The exact twelve-port icosahedral carrier and its canonical repair law.

The carrier is the oriented ``(12, 30, 20)`` incidence of
``oph_fpe.dynamics.self_readback_repair_closure.ORIENTED_BASE_FACES``.  Its
sixty proper rotations are derived combinatorially from the oriented faces and
coincide with the coordinate group of ``oph_fpe.core.icosahedral``.  The
antipodal (inverse-port) pairing is graph distance three.

The canonical repair law is the seam-mean retraction of the canonical repair
law document: on seam ``e = {i, j}`` the two endpoint readings are replaced by
their mean, ``E_e = I - b_e b_e^T / 2`` with ``b_e = e_i - e_j``.  Uniform
scheduling over the thirty seams gives the repair mean ``T = I - L/60`` with
``L = 5I - A`` the seam Laplacian.  The centered infinite-response kernel
``C_n = Q T^{2n} Q`` normalizes to the intrinsic Gram ``G = 4 P_slow`` of rank
three, where ``P_slow`` is the projector onto the ``5 - sqrt(5)`` Laplacian
band.  All of this is exact over ``Q(sqrt 5)``; floating arrays are provided
for dynamics at scale.

The integer law is the nearest-agreement retraction of the bounded
self-readback closure: ``(x_i, x_j) -> (floor(s/2), ceil(s/2))`` with ``s`` the
endpoint total, the odd tie placed by the directed seam label.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
from typing import Iterable, Sequence

import numpy as np

from oph_fpe.dynamics.self_readback_repair_closure import (
    ORIENTED_BASE_FACES,
    exact_reference_edges,
)

PORT_COUNT = 12
SEAM_COUNT = 30
FACE_COUNT = 20
ROTATION_COUNT = 60
SCHEDULE_DENOMINATOR = 60  # T = I - L/60


# --------------------------------------------------------------------------
# Q(sqrt 5) arithmetic: a number is (a, b) meaning a + b*sqrt(5), a, b in Q.
# --------------------------------------------------------------------------

QS5 = tuple[Fraction, Fraction]


def q5(a: int | Fraction, b: int | Fraction = 0) -> QS5:
    return (Fraction(a), Fraction(b))


def q5_add(x: QS5, y: QS5) -> QS5:
    return (x[0] + y[0], x[1] + y[1])


def q5_sub(x: QS5, y: QS5) -> QS5:
    return (x[0] - y[0], x[1] - y[1])


def q5_mul(x: QS5, y: QS5) -> QS5:
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def q5_inv(x: QS5) -> QS5:
    norm = x[0] * x[0] - 5 * x[1] * x[1]
    if norm == 0:
        raise ZeroDivisionError("zero element of Q(sqrt5)")
    return (x[0] / norm, -x[1] / norm)


def q5_sign(x: QS5) -> int:
    """Exact sign of a + b sqrt5."""

    a, b = x
    if b == 0:
        return (a > 0) - (a < 0)
    if a >= 0 and b > 0:
        return 1
    if a <= 0 and b < 0:
        return -1
    d = a * a - 5 * b * b
    return ((d > 0) - (d < 0)) * (1 if a > 0 else -1)


def q5_float(x: QS5) -> float:
    return float(x[0]) + float(x[1]) * 5.0**0.5


def q5_str(x: QS5) -> list[str]:
    return [str(x[0]), str(x[1])]


GOLDEN: QS5 = (Fraction(1, 2), Fraction(1, 2))  # (1 + sqrt5)/2
SQRT5: QS5 = (Fraction(0), Fraction(1))


def q5_matmul(x: Sequence[Sequence[QS5]], y: Sequence[Sequence[QS5]]) -> list[list[QS5]]:
    n, k, m = len(x), len(y), len(y[0])
    out = [[q5(0) for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for t in range(k):
            xit = x[i][t]
            if xit == (0, 0):
                continue
            row = y[t]
            for j in range(m):
                out[i][j] = q5_add(out[i][j], q5_mul(xit, row[j]))
    return out


def q5_matrix_to_float(x: Sequence[Sequence[QS5]]) -> np.ndarray:
    return np.array([[q5_float(v) for v in row] for row in x], dtype=float)


# --------------------------------------------------------------------------
# Incidence
# --------------------------------------------------------------------------


@lru_cache(maxsize=1)
def seams() -> tuple[tuple[int, int], ...]:
    """The thirty unoriented seams ``(i, j)`` with ``i < j``."""

    return exact_reference_edges()


@lru_cache(maxsize=1)
def oriented_faces() -> tuple[tuple[int, int, int], ...]:
    return tuple(tuple(int(v) for v in face) for face in ORIENTED_BASE_FACES)


@lru_cache(maxsize=1)
def adjacency() -> np.ndarray:
    a = np.zeros((PORT_COUNT, PORT_COUNT), dtype=np.int64)
    for i, j in seams():
        a[i, j] = 1
        a[j, i] = 1
    a.setflags(write=False)
    return a


@lru_cache(maxsize=1)
def laplacian() -> np.ndarray:
    lap = 5 * np.eye(PORT_COUNT, dtype=np.int64) - adjacency()
    lap.setflags(write=False)
    return lap


@lru_cache(maxsize=1)
def graph_distance() -> np.ndarray:
    a = adjacency().astype(bool)
    dist = np.full((PORT_COUNT, PORT_COUNT), -1, dtype=np.int64)
    np.fill_diagonal(dist, 0)
    reach = np.eye(PORT_COUNT, dtype=bool)
    for step in range(1, PORT_COUNT):
        reach_next = reach | (reach.astype(np.int64) @ a.astype(np.int64) > 0)
        newly = reach_next & ~reach
        dist[newly] = step
        reach = reach_next
        if reach.all():
            break
    dist.setflags(write=False)
    return dist


@lru_cache(maxsize=1)
def antipode() -> tuple[int, ...]:
    """The unique fixed-point-free inverse-port pairing, graph distance three."""

    dist = graph_distance()
    pairing = []
    for p in range(PORT_COUNT):
        far = [q for q in range(PORT_COUNT) if dist[p, q] == 3]
        if len(far) != 1:
            raise AssertionError("distance-three partner is not unique")
        pairing.append(far[0])
    pairing_t = tuple(pairing)
    if any(pairing_t[pairing_t[p]] != p or pairing_t[p] == p for p in range(PORT_COUNT)):
        raise AssertionError("antipodal pairing is not a fixed-point-free involution")
    return pairing_t


def _cyclic(face: Sequence[int]) -> tuple[int, ...]:
    face = tuple(face)
    index = face.index(min(face))
    return face[index:] + face[:index]


@lru_cache(maxsize=1)
def incidence_automorphisms() -> tuple[tuple[int, ...], ...]:
    """All 120 permutations of the ports preserving the seam set."""

    seam_set = set(seams())
    found = []
    # Backtracking over port images with adjacency consistency.
    adj = [set() for _ in range(PORT_COUNT)]
    for i, j in seam_set:
        adj[i].add(j)
        adj[j].add(i)

    def extend(assignment: list[int], used: set[int]) -> None:
        p = len(assignment)
        if p == PORT_COUNT:
            found.append(tuple(assignment))
            return
        for image in range(PORT_COUNT):
            if image in used:
                continue
            consistent = True
            for q in range(p):
                if (q in adj[p]) != (assignment[q] in adj[image]):
                    consistent = False
                    break
            if consistent:
                assignment.append(image)
                used.add(image)
                extend(assignment, used)
                assignment.pop()
                used.discard(image)

    extend([], set())
    if len(found) != 120:
        raise AssertionError("expected 120 incidence automorphisms")
    return tuple(found)


@lru_cache(maxsize=1)
def rotations() -> tuple[tuple[int, ...], ...]:
    """The sixty proper rotations: automorphisms preserving oriented faces."""

    faces = {_cyclic(face) for face in oriented_faces()}
    proper = tuple(
        perm
        for perm in incidence_automorphisms()
        if all(_cyclic((perm[a], perm[b], perm[c])) in faces for a, b, c in oriented_faces())
    )
    if len(proper) != ROTATION_COUNT:
        raise AssertionError("expected sixty proper rotations")
    return proper


def permutation_matrix(perm: Sequence[int]) -> np.ndarray:
    m = np.zeros((PORT_COUNT, PORT_COUNT), dtype=np.int64)
    for old, new in enumerate(perm):
        m[new, old] = 1
    return m


# --------------------------------------------------------------------------
# Canonical repair law
# --------------------------------------------------------------------------


def seam_mean_matrix(seam: tuple[int, int]) -> np.ndarray:
    """``E_e = I - b_e b_e^T / 2`` as a float matrix."""

    i, j = seam
    b = np.zeros(PORT_COUNT)
    b[i] = 1.0
    b[j] = -1.0
    return np.eye(PORT_COUNT) - 0.5 * np.outer(b, b)


def seam_mean_matrix_exact(seam: tuple[int, int]) -> list[list[Fraction]]:
    i, j = seam
    m = [[Fraction(int(r == c)) for c in range(PORT_COUNT)] for r in range(PORT_COUNT)]
    half = Fraction(1, 2)
    m[i][i] -= half
    m[j][j] -= half
    m[i][j] += half
    m[j][i] += half
    return m


def repair_mean() -> np.ndarray:
    """``T = I - L/60`` (float)."""

    return np.eye(PORT_COUNT) - laplacian().astype(float) / SCHEDULE_DENOMINATOR


def repair_mean_exact() -> list[list[Fraction]]:
    lap = laplacian()
    return [
        [Fraction(int(r == c)) - Fraction(int(lap[r, c]), SCHEDULE_DENOMINATOR) for c in range(PORT_COUNT)]
        for r in range(PORT_COUNT)
    ]


def apply_seam_mean(x: np.ndarray, seam: tuple[int, int]) -> np.ndarray:
    """Return the state after the canonical mean repair on one seam."""

    i, j = seam
    y = np.array(x, dtype=float, copy=True)
    m = 0.5 * (y[i] + y[j])
    y[i] = m
    y[j] = m
    return y


def integer_nearest_agreement(xi: int, xj: int, *, ceiling_to_first: bool) -> tuple[int, int]:
    """Integer law: preserve the total, land in the nearest agreement shell.

    For an even total both endpoints receive ``s/2``.  For an odd total the
    ceiling goes to the first endpoint when ``ceiling_to_first`` is true and
    to the second otherwise; the two placements are the two directed labels of
    the seam.
    """

    s = int(xi) + int(xj)
    low, high = s // 2, -((-s) // 2)
    if ceiling_to_first:
        return high, low
    return low, high


def mismatch_potential(x: np.ndarray) -> float:
    """``Phi(x) = sum_seams (x_i - x_j)^2`` on one carrier."""

    x = np.asarray(x, dtype=float)
    return float(sum((x[i] - x[j]) ** 2 for i, j in seams()))


# --------------------------------------------------------------------------
# Spectral structure over Q(sqrt 5)
# --------------------------------------------------------------------------

LAPLACIAN_BANDS: tuple[tuple[QS5, int], ...] = (
    (q5(0), 1),
    (q5(5, -1), 3),  # 5 - sqrt5 : the slow band
    (q5(6), 5),
    (q5(5, 1), 3),  # 5 + sqrt5
)


def _q5_identity() -> list[list[QS5]]:
    return [[q5(int(r == c)) for c in range(PORT_COUNT)] for r in range(PORT_COUNT)]


def _q5_from_int_matrix(m: np.ndarray) -> list[list[QS5]]:
    return [[q5(int(m[r, c])) for c in range(PORT_COUNT)] for r in range(PORT_COUNT)]


@lru_cache(maxsize=None)
def band_projector_exact(band_index: int) -> tuple[tuple[QS5, ...], ...]:
    """Exact spectral projector of the seam Laplacian onto one band.

    ``P_k = prod_{m != k} (L - mu_m I) / (mu_k - mu_m)``.
    """

    mu_k = LAPLACIAN_BANDS[band_index][0]
    lap = _q5_from_int_matrix(laplacian())
    result = _q5_identity()
    for m, (mu_m, _) in enumerate(LAPLACIAN_BANDS):
        if m == band_index:
            continue
        factor = [
            [q5_sub(lap[r][c], mu_m) if r == c else lap[r][c] for c in range(PORT_COUNT)]
            for r in range(PORT_COUNT)
        ]
        result = q5_matmul(result, factor)
        scale = q5_inv(q5_sub(mu_k, mu_m))
        result = [[q5_mul(v, scale) for v in row] for row in result]
    return tuple(tuple(row) for row in result)


def slow_band_projector_exact() -> tuple[tuple[QS5, ...], ...]:
    return band_projector_exact(1)


def slow_band_projector() -> np.ndarray:
    return q5_matrix_to_float(slow_band_projector_exact())


def intrinsic_gram_exact() -> tuple[tuple[QS5, ...], ...]:
    """``G = 4 P_slow``: entries in ``{1, 1/sqrt5, -1/sqrt5, -1}`` by distance."""

    p = slow_band_projector_exact()
    return tuple(tuple(q5_mul(q5(4), v) for v in row) for row in p)


def intrinsic_gram() -> np.ndarray:
    return 4.0 * slow_band_projector()


def centered_response_kernel(n: int) -> np.ndarray:
    """``C_n = Q T^{2n} Q`` with ``Q = I - J/12`` (float)."""

    q = np.eye(PORT_COUNT) - np.ones((PORT_COUNT, PORT_COUNT)) / PORT_COUNT
    t = repair_mean()
    return q @ np.linalg.matrix_power(t, 2 * n) @ q


def normalized_response_kernel(n: int) -> np.ndarray:
    """``12 C_n / tr C_n``; converges to ``4 P_slow`` as ``n`` grows.

    Computed from the exact band projectors with the band powers taken
    relative to the slow band, so the value is accurate at every ``n``; the
    matrix-power form ``centered_response_kernel`` loses relative precision
    once ``tr C_n`` falls near the floating-point floor.
    """

    damping = band_damping_factors()
    slow = damping["slow_3"]
    kernel = np.zeros((PORT_COUNT, PORT_COUNT))
    trace = 0.0
    for band_index, key in ((1, "slow_3"), (2, "middle_5"), (3, "fast_3")):
        share = (damping[key] / slow) ** (2 * n)
        projector = q5_matrix_to_float(band_projector_exact(band_index))
        kernel += share * projector
        trace += share * LAPLACIAN_BANDS[band_index][1]
    return PORT_COUNT * kernel / trace


def band_damping_factors() -> dict[str, float]:
    """Per-band eigenvalues of ``T``: ``1 - mu/60``."""

    return {
        "constant": 1.0,
        "slow_3": q5_float(q5_sub(q5(1), q5_mul(q5(Fraction(1, 60)), LAPLACIAN_BANDS[1][0]))),
        "middle_5": 1.0 - 6.0 / 60.0,
        "fast_3": q5_float(q5_sub(q5(1), q5_mul(q5(Fraction(1, 60)), LAPLACIAN_BANDS[3][0]))),
    }


def a5_character_on_ports(perm: Sequence[int]) -> int:
    """Trace of the permutation on ``R^12`` (number of fixed ports)."""

    return sum(int(perm[p] == p) for p in range(PORT_COUNT))


@dataclass(frozen=True)
class CarrierSummary:
    ports: int
    seams: int
    faces: int
    rotations: int
    antipode: tuple[int, ...]
    laplacian_bands: tuple[tuple[list[str], int], ...]
    damping: dict[str, float]

    @staticmethod
    def build() -> "CarrierSummary":
        return CarrierSummary(
            ports=PORT_COUNT,
            seams=len(seams()),
            faces=len(oriented_faces()),
            rotations=len(rotations()),
            antipode=antipode(),
            laplacian_bands=tuple((q5_str(mu), mult) for mu, mult in LAPLACIAN_BANDS),
            damping=band_damping_factors(),
        )


def check_carrier() -> dict[str, object]:
    """Exact self-checks; every value must be true."""

    lap = laplacian()
    p_slow = slow_band_projector_exact()
    p_float = q5_matrix_to_float(p_slow)
    idempotent = np.allclose(p_float @ p_float, p_float, atol=1e-12)
    trace_three = abs(np.trace(p_float) - 3.0) < 1e-12
    eig_ok = np.allclose(lap.astype(float) @ p_float, q5_float(LAPLACIAN_BANDS[1][0]) * p_float, atol=1e-12)
    gram = intrinsic_gram()
    dist = graph_distance()
    expected = {0: 1.0, 1: 5 ** -0.5, 2: -(5 ** -0.5), 3: -1.0}
    gram_by_distance = all(
        abs(gram[p, q] - expected[int(dist[p, q])]) < 1e-12 for p in range(PORT_COUNT) for q in range(PORT_COUNT)
    )
    limit = normalized_response_kernel(200)
    limit_ok = np.allclose(limit, gram, atol=1e-9)
    rank_three = int(np.linalg.matrix_rank(gram, tol=1e-9)) == 3
    rot = rotations()
    closed = all(tuple(a[b[p]] for p in range(PORT_COUNT)) in set(rot) for a in rot[:6] for b in rot)
    anti = antipode()
    anti_commutes = all(tuple(perm[anti[p]] for p in range(PORT_COUNT)) == tuple(anti[perm[p]] for p in range(PORT_COUNT)) for perm in rot)
    return {
        "slow_projector_idempotent": bool(idempotent),
        "slow_projector_trace_three": bool(trace_three),
        "slow_projector_is_laplacian_eigenprojector": bool(eig_ok),
        "gram_entries_by_distance": bool(gram_by_distance),
        "normalized_kernel_limit_equals_gram": bool(limit_ok),
        "gram_rank_three": bool(rank_three),
        "rotation_sample_closed": bool(closed),
        "antipode_commutes_with_rotations": bool(anti_commutes),
    }


__all__ = [
    "PORT_COUNT",
    "SEAM_COUNT",
    "FACE_COUNT",
    "ROTATION_COUNT",
    "QS5",
    "q5",
    "q5_add",
    "q5_sub",
    "q5_mul",
    "q5_inv",
    "q5_sign",
    "q5_float",
    "q5_str",
    "GOLDEN",
    "SQRT5",
    "seams",
    "oriented_faces",
    "adjacency",
    "laplacian",
    "graph_distance",
    "antipode",
    "incidence_automorphisms",
    "rotations",
    "permutation_matrix",
    "seam_mean_matrix",
    "seam_mean_matrix_exact",
    "repair_mean",
    "repair_mean_exact",
    "apply_seam_mean",
    "integer_nearest_agreement",
    "mismatch_potential",
    "LAPLACIAN_BANDS",
    "band_projector_exact",
    "slow_band_projector_exact",
    "slow_band_projector",
    "intrinsic_gram_exact",
    "intrinsic_gram",
    "centered_response_kernel",
    "normalized_response_kernel",
    "band_damping_factors",
    "a5_character_on_ports",
    "CarrierSummary",
    "check_carrier",
]
