"""Exact federation of twelve-port carriers on the geodesic icosahedral tower.

Lane L1 of the exact-carrier package.  ``N = 20 * 4^L`` carriers sit on the
cells of the geodesic icosahedral tower at level ``L``.  Every carrier runs
the canonical seam-mean law on its thirty icosahedral seams
(``oph_exact.carrier``).  Two gluing modes are declared:

``isolated``
    no inter-carrier seams; the control, on which every federation quantity
    reduces to the single-carrier objects of ``oph_exact.carrier``.
``port_pair``
    the production convention of
    ``oph_fpe.core.screen_ports.assign_echosahedral_ports``: cells are
    adjacent through the cell-dual graph (three edge neighbours per cell),
    and each endpoint of a dual edge is routed to the port whose icosahedral
    direction, read in the cell's local tangent frame, has the largest inner
    product with the neighbour direction (collisions repaired by an exact
    assignment).  Every routed port pair is one inter-carrier seam.  This is a
    declared convention; the source-derived gluing is the open (M1) item of
    the theory and the receipt measures what this convention does.

Dynamics.  The asynchronous schedule draws seams uniformly over all seams
(intra and inter) with a seeded numpy generator; one sweep is ``|S|``
attempts.  Every accepted move is the seam-mean retraction (both endpoints
replaced by their mean).  The descent functional is ``V(x) = sum_p x_p^2``:
the seam mean is the orthogonal projection onto the seam equalizer, so
``V`` drops by exactly ``(x_i - x_j)^2 / 2`` per non-wait move.  The mismatch
potential ``Phi(x) = sum_seams (x_i - x_j)^2`` is the termination potential
(``Phi = 0`` is consensus); it is not monotone along single moves and the
receipt records an exact witness of that.  The synchronous expectation
operator is ``T_fed = I - L_fed / D`` with ``D = 2|S|/N`` (one attempt per
carrier per tick); on the isolated mode ``D = 60`` and every carrier block is
``carrier.repair_mean() = I - L_ico/60`` exactly.

The integer law uses ``carrier.integer_nearest_agreement`` with the odd-total
tie placed by a fair coin (A3).  Its terminal class is the balanced class
(every reading in ``{q, q+1}`` per component), which is one orbit of the
neutral unit swaps; the canonical quotient is the multiset of readings per
component.

Three float engines (pure Python reference, layered numpy, native C through
ctypes) apply identical IEEE operations in a dependency-respecting order and
produce bit-identical states; the receipt is engine independent.  Exact
arithmetic on the small rungs uses dyadic integers ``x_p = n_p / 2^K``.

Receipt: ``data/exact/federation_canonical_mean_receipt.json`` written by
``python -m oph_exact.federation --write`` and checked by ``--check``.  The
independent verifier is ``oph_exact/verify_federation_independent.py``.
"""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.csgraph as csgraph
import scipy.sparse.linalg as sparse_linalg

from oph_exact import carrier
from oph_fpe.core.icosahedral import (
    build_geodesic_icosahedral_tower,
    geodesic_icosahedral_patch_arrays,
)
from oph_fpe.core.screen_ports import assign_echosahedral_ports

SCHEMA = "oph.exact.federation-canonical-mean.v1"
REPO_ROOT = Path(__file__).resolve().parents[1]
RECEIPT_PATH = REPO_ROOT / "data" / "exact" / "federation_canonical_mean_receipt.json"
PRODUCER_PATH = Path(__file__).resolve()
VERIFIER_PATH = PRODUCER_PATH.parent / "verify_federation_independent.py"
TEST_PATH = REPO_ROOT / "tests" / "test_exact_federation.py"
CARRIER_PATH = PRODUCER_PATH.parent / "carrier.py"
PRODUCTION_LAW_REPORT = "runs/e6_64k_dense_20260820/finite_consensus_replay_report.json"

PORTS = carrier.PORT_COUNT
GLUINGS = ("isolated", "port_pair")
RECEIPT_LEVELS = (0, 1, 2, 3)
EXACT_LEVELS = (0, 1)
SCHEDULES = 16
KERNEL_STEPS = (1, 5, 30, 100, 300)
KERNEL_SAMPLE = 12
FLOAT_PHI_THRESHOLD = 1e-18
EXACT_PHI_THRESHOLD_DENOMINATOR = 10**18  # Phi <= 1 / 10^18 exactly
EXACT_TRIGGER_FACTOR = 1.1  # exact Phi is evaluated once the float companion is below factor * threshold
LATTICE_SNAP_MARGIN = 0.25  # a component-lattice snap is unambiguous only when every residual |m_c x_p - q_p| is below this
LOAD_SEED_BASE = 20260909
SCHEDULE_SEED_BASE = 909000
LOAD_MAX = 5
MAX_SWEEPS_DEFAULT = 400_000
PHI_TRACE_ROUND = 12
FLOAT_ROUND = 12


# --------------------------------------------------------------------------
# Canonical JSON and hashing
# --------------------------------------------------------------------------


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False) + "\n"


def sha256_of(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("ascii")).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _round(value: float, digits: int = FLOAT_ROUND) -> float:
    out = float(round(float(value), digits))
    return 0.0 if out == 0.0 else out


def _sig(value: float, digits: int = 6) -> float:
    """Round to ``digits`` significant figures (for quantities near the float floor)."""

    value = float(value)
    if value == 0.0 or not math.isfinite(value):
        return value
    scale = digits - 1 - int(math.floor(math.log10(abs(value))))
    return float(round(value, scale))


# --------------------------------------------------------------------------
# Federation geometry
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class PortGraph:
    """Seam endpoints and the neighbour structure shared by the move engines."""

    a: np.ndarray
    b: np.ndarray
    indptr: np.ndarray
    indices: np.ndarray
    degree: np.ndarray
    neighbour_lists: tuple[tuple[int, ...], ...]
    neighbour_table: np.ndarray  # (ports, max degree), -1 padded, CSR order

    @staticmethod
    def build(a: np.ndarray, b: np.ndarray, ports: int) -> "PortGraph":
        a = np.ascontiguousarray(a, dtype=np.int64)
        b = np.ascontiguousarray(b, dtype=np.int64)
        adjacency = sparse.coo_matrix((np.ones(a.size), (a, b)), shape=(ports, ports))
        adjacency = (adjacency + adjacency.T).tocsr()
        adjacency.sort_indices()
        indptr = np.ascontiguousarray(adjacency.indptr, dtype=np.int64)
        indices = np.ascontiguousarray(adjacency.indices, dtype=np.int64)
        degree = np.diff(indptr).astype(np.int64)
        width = int(degree.max()) if degree.size else 0
        table = np.full((ports, width), -1, dtype=np.int64)
        lists = []
        for port in range(ports):
            row = indices[indptr[port] : indptr[port + 1]]
            table[port, : row.size] = row
            lists.append(tuple(int(v) for v in row))
        return PortGraph(a=a, b=b, indptr=indptr, indices=indices, degree=degree, neighbour_lists=tuple(lists), neighbour_table=table)


@dataclass(frozen=True)
class ExactFederation:
    """``N`` carriers on the tower cells with a declared gluing."""

    level: int
    gluing: str
    carriers: int
    seam_a: np.ndarray
    seam_b: np.ndarray
    graph: PortGraph
    intra_count: int
    inter_pairs: np.ndarray  # (|inter|, 4): cell_a, port_a, cell_b, port_b
    component_of_port: np.ndarray
    component_sizes: np.ndarray
    laplacian: sparse.csr_matrix
    pentagonal_cell: np.ndarray
    gluing_receipt: dict[str, Any]

    @property
    def ports(self) -> int:
        return PORTS * self.carriers

    @property
    def seams(self) -> int:
        return int(self.seam_a.size)

    @property
    def inter_count(self) -> int:
        return int(self.inter_pairs.shape[0])

    @property
    def mean_denominator(self) -> Fraction:
        """``D = 2|S|/N``: ``T_fed = I - L_fed/D``."""

        return Fraction(2 * self.seams, self.carriers)

    @property
    def components(self) -> int:
        return int(self.component_sizes.size)

    def component_means(self, x: np.ndarray) -> np.ndarray:
        totals = np.bincount(self.component_of_port, weights=np.asarray(x, dtype=float), minlength=self.components)
        return totals / self.component_sizes

    def mismatch_potential(self, x: np.ndarray) -> float:
        x = np.asarray(x, dtype=float)
        d = x[self.seam_a] - x[self.seam_b]
        return float(np.dot(d, d))

    def descent_functional(self, x: np.ndarray) -> float:
        x = np.asarray(x, dtype=float)
        return float(np.dot(x, x))

    def apply_synchronous(self, y: np.ndarray) -> np.ndarray:
        """``T_fed y = y - (L_fed y)/D``."""

        return y - (self.laplacian @ y) / float(self.mean_denominator)

    def terminal_quotient_hash(self, x: np.ndarray) -> tuple[str, float]:
        """Snap every reading to its component lattice ``(1/m_c) Z`` and hash.

        Returns the hash and the maximal lattice residual ``|m_c x_p - q_p|``.
        """

        x = np.asarray(x, dtype=float)
        m = self.component_sizes[self.component_of_port].astype(float)
        scaled = m * x
        q = np.rint(scaled)
        residual = float(np.max(np.abs(scaled - q))) if scaled.size else 0.0
        payload = {
            "canonicalizer": "component_lattice_snap",
            "component_sizes": self.component_sizes.tolist(),
            "component_of_port": self.component_of_port.tolist(),
            "q": q.astype(np.int64).tolist(),
        }
        return sha256_of(payload), residual

    def expected_terminal_quotient_hash(self, loads: np.ndarray) -> str:
        totals = np.bincount(self.component_of_port, weights=np.asarray(loads, dtype=float), minlength=self.components)
        q = np.rint(totals[self.component_of_port]).astype(np.int64)
        payload = {
            "canonicalizer": "component_lattice_snap",
            "component_sizes": self.component_sizes.tolist(),
            "component_of_port": self.component_of_port.tolist(),
            "q": q.tolist(),
        }
        return sha256_of(payload)

    def integer_quotient_hash(self, x: np.ndarray) -> str:
        """The multiset of readings per component (components ordered by lowest port)."""

        x = np.asarray(x, dtype=np.int64)
        order = np.argsort(self.component_of_port, kind="stable")
        comps = self.component_of_port[order]
        vals = x[order]
        starts = np.flatnonzero(np.r_[True, comps[1:] != comps[:-1]])
        stops = np.r_[starts[1:], comps.size]
        entries = []
        for start, stop in zip(starts, stops):
            block = np.sort(vals[start:stop])
            values, counts = np.unique(block, return_counts=True)
            entries.append([int(stop - start), [[int(v), int(c)] for v, c in zip(values, counts)]])
        return sha256_of({"canonicalizer": "component_multiset", "components": entries})

    def expected_integer_quotient_hash(self, loads: np.ndarray) -> str:
        totals = np.bincount(self.component_of_port, weights=np.asarray(loads, dtype=float), minlength=self.components)
        entries = []
        for c in range(self.components):
            m = int(self.component_sizes[c])
            total = int(round(float(totals[c])))
            q, r = divmod(total, m)
            multiset = []
            if m - r:
                multiset.append([q, m - r])
            if r:
                multiset.append([q + 1, r])
            entries.append([m, multiset])
        return sha256_of({"canonicalizer": "component_multiset", "components": entries})

    def mean_minimum_descent(self, loads: np.ndarray) -> float:
        """``V`` at the component mean: ``sum_c S_c^2 / m_c``, exact in integer arithmetic before the final rounding.

        The integer loads make every component total an integer; squaring in float64 would
        lose exactness once a total exceeds ``2^26.5`` (level ten and beyond).
        """

        totals = np.bincount(self.component_of_port, weights=np.asarray(loads, dtype=float), minlength=self.components)
        exact = sum(Fraction(int(round(float(t))) ** 2, int(m)) for t, m in zip(totals.tolist(), self.component_sizes.tolist()))
        return float(exact)

    def integer_minimum_descent(self, loads: np.ndarray) -> int:
        totals = np.bincount(self.component_of_port, weights=np.asarray(loads, dtype=float), minlength=self.components)
        v_min = 0
        for c in range(self.components):
            m = int(self.component_sizes[c])
            total = int(round(float(totals[c])))
            q, r = divmod(total, m)
            v_min += (m - r) * q * q + r * (q + 1) * (q + 1)
        return int(v_min)


def _intra_seams(carriers: int) -> tuple[np.ndarray, np.ndarray]:
    local = np.asarray(carrier.seams(), dtype=np.int64)
    base = (np.arange(carriers, dtype=np.int64) * PORTS)[:, None]
    a = (base + local[None, :, 0]).ravel()
    b = (base + local[None, :, 1]).ravel()
    return a, b


def _pentagonal_cells(level: int) -> np.ndarray:
    faces = build_geodesic_icosahedral_tower(level).levels[level].faces
    return np.any(np.asarray(faces) < 12, axis=1)


@lru_cache(maxsize=8)
def build_federation(level: int, gluing: str) -> ExactFederation:
    if gluing not in GLUINGS:
        raise ValueError(f"unknown gluing {gluing!r}; declared gluings are {GLUINGS}")
    points, left, right = geodesic_icosahedral_patch_arrays(level, patch_basis="cells")
    carriers = int(points.shape[0])
    if carriers != 20 * 4**level:
        raise AssertionError("cell count differs from 20 * 4^L")
    intra_a, intra_b = _intra_seams(carriers)
    antipode = carrier.antipode()
    if gluing == "isolated":
        inter = np.zeros((0, 4), dtype=np.int64)
        gluing_receipt = {
            "mode": "isolated",
            "description": "no inter-carrier seams; every carrier is its own component",
            "inter_seam_count": 0,
        }
    else:
        port_map = assign_echosahedral_ports(left, right, carriers, points=points)
        if port_map.routing_mode != "icosahedral_directional_assignment":
            raise AssertionError("production geometric routing was not selected")
        if port_map.overflow_count != 0:
            raise AssertionError("production routing overflowed a twelve-slot carrier")
        lp = port_map.left_port.astype(np.int64)
        rp = port_map.right_port.astype(np.int64)
        inter = np.stack([left.astype(np.int64), lp, right.astype(np.int64), rp], axis=1)
        slots = np.concatenate([left * PORTS + lp, right * PORTS + rp])
        if np.unique(slots).size != slots.size:
            raise AssertionError("a carrier port is glued more than once")
        alignment = np.asarray(port_map.directional_alignment, dtype=float)
        usage = np.bincount(np.concatenate([lp, rp]), minlength=PORTS)
        antipodal = int(np.sum(np.asarray([antipode[int(p)] for p in lp], dtype=np.int64) == rp))
        gluing_receipt = {
            "mode": "port_pair",
            "convention": "oph_fpe.core.screen_ports.assign_echosahedral_ports(left, right, N, points=cell_centres)",
            "routing_mode": port_map.routing_mode,
            "description": (
                "cells adjacent through the cell-dual graph of "
                "geodesic_icosahedral_patch_arrays(level, patch_basis='cells') (three edge "
                "neighbours per cell); each dual-edge endpoint is routed to the port whose "
                "icosahedral template direction, read in the cell's local tangent frame "
                "(tangent_x = normalize(reference x normal) with reference (0,0,1), or (1,0,0) "
                "when |normal_z| > 0.9; tangent_y = normal x tangent_x), has the largest "
                "float32 inner product with the unit vector from the cell centre to the "
                "neighbour centre; duplicate choices inside one cell are repaired by an exact "
                "maximum-alignment assignment; each routed pair is one inter-carrier seam"
            ),
            "declared_convention": True,
            "source_derived": False,
            "inter_seam_count": int(inter.shape[0]),
            "glued_ports_per_carrier": 3,
            "port_usage_histogram": usage.tolist(),
            "antipodal_consistent_pairs": antipodal,
            "antipodal_consistent_fraction": _round(antipodal / inter.shape[0]),
            "alignment_min": _round(float(alignment.min())),
            "alignment_mean": _round(float(alignment.mean())),
            "local_frame_hash": port_map.local_frame_hash,
            "port_pairs_sha256": sha256_of(inter.tolist()),
        }
        if level <= 1:
            gluing_receipt["port_pairs"] = inter.tolist()
    seam_a = np.concatenate([intra_a, inter[:, 0] * PORTS + inter[:, 1]]).astype(np.int64)
    seam_b = np.concatenate([intra_b, inter[:, 2] * PORTS + inter[:, 3]]).astype(np.int64)
    ports = PORTS * carriers
    adjacency = sparse.coo_matrix(
        (np.ones(seam_a.size, dtype=float), (seam_a, seam_b)), shape=(ports, ports)
    )
    adjacency = (adjacency + adjacency.T).tocsr()
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    laplacian = (sparse.diags(degree) - adjacency).tocsr()
    count, labels = csgraph.connected_components(adjacency, directed=False)
    # Relabel components by lowest port index so the labelling is canonical.
    first_port = np.full(count, ports, dtype=np.int64)
    np.minimum.at(first_port, labels, np.arange(ports, dtype=np.int64))
    relabel = np.empty(count, dtype=np.int64)
    relabel[np.argsort(first_port, kind="stable")] = np.arange(count, dtype=np.int64)
    labels = relabel[labels].astype(np.int64)
    sizes = np.bincount(labels, minlength=count).astype(np.int64)
    seam_a.setflags(write=False)
    seam_b.setflags(write=False)
    return ExactFederation(
        level=level,
        gluing=gluing,
        carriers=carriers,
        seam_a=seam_a,
        seam_b=seam_b,
        graph=PortGraph.build(seam_a, seam_b, ports),
        intra_count=int(intra_a.size),
        inter_pairs=inter,
        component_of_port=labels,
        component_sizes=sizes,
        laplacian=laplacian,
        pentagonal_cell=_pentagonal_cells(level),
        gluing_receipt=gluing_receipt,
    )


def initial_loads(level: int, ports: int) -> np.ndarray:
    """Seeded integer loads in ``{0..5}``; identical for every gluing at a level."""

    rng = np.random.default_rng(LOAD_SEED_BASE + level)
    return rng.integers(0, LOAD_MAX + 1, size=ports).astype(np.int64)


def schedule_seed(level: int, gluing: str, index: int) -> int:
    return SCHEDULE_SEED_BASE + 1000 * level + 100 * GLUINGS.index(gluing) + index


# --------------------------------------------------------------------------
# Native kernel (optional, bit-identical to the reference engines)
# --------------------------------------------------------------------------

_KERNEL_SOURCE = r"""
#include <stdint.h>

/* Mean law.  counts = {waits, descent violations, seam-sum raising moves,
   first raising attempt, first raising seam}; witness = {x_i, x_j, N_i, N_j,
   deg_i, deg_j, delta_phi} of the first raising move. */
void oph_fed_mean(double *x, const int64_t *a, const int64_t *b,
                  const int64_t *indptr, const int64_t *indices,
                  const int64_t *seq, int64_t n, double *ledger,
                  int64_t *counts, double *witness)
{
    int64_t waits = 0, violations = 0, phi_up = 0;
    double led = 0.0;
    for (int64_t t = 0; t < n; ++t) {
        const int64_t s = seq[t];
        const int64_t i = a[s];
        const int64_t j = b[s];
        const double xi = x[i];
        const double xj = x[j];
        const double d = xi - xj;
        if (d == 0.0) { ++waits; continue; }
        const double half_sq = 0.5 * (d * d);
        if (!(half_sq > 0.0)) ++violations;
        led += half_sq;
        double ni = 0.0, nj = 0.0;
        for (int64_t k = indptr[i]; k < indptr[i + 1]; ++k) ni += x[indices[k]];
        for (int64_t k = indptr[j]; k < indptr[j + 1]; ++k) nj += x[indices[k]];
        const double degi = (double)(indptr[i + 1] - indptr[i]);
        const double degj = (double)(indptr[j + 1] - indptr[j]);
        const double delta = -d * ((degi * xi - ni) - (degj * xj - nj)) + (d * d) * ((degi + degj + 2.0) * 0.25);
        if (delta > 0.0) {
            if (phi_up == 0 && counts[3] < 0) {
                counts[3] = t;
                counts[4] = s;
                witness[0] = xi; witness[1] = xj; witness[2] = ni; witness[3] = nj;
                witness[4] = degi; witness[5] = degj; witness[6] = delta;
            }
            ++phi_up;
        }
        const double m = 0.5 * (xi + xj);
        x[i] = m;
        x[j] = m;
    }
    ledger[0] += led;
    counts[0] += waits;
    counts[1] += violations;
    counts[2] += phi_up;
}

/* Integer law.  counts = {descents, swaps, waits, unit transfers, decrement
   identity violations}.  Returns the index of the first move after which
   V equals v_min, or -1. */
int64_t oph_fed_integer(int64_t *x, const int64_t *a, const int64_t *b,
                        const int64_t *seq, const int64_t *coin, int64_t n,
                        int64_t *counts, int64_t *v, int64_t v_min)
{
    int64_t first = -1;
    int64_t vv = v[0];
    int64_t descents = 0, swaps = 0, waits = 0, transfers = 0, violations = 0;
    for (int64_t t = 0; t < n; ++t) {
        const int64_t s = seq[t];
        const int64_t i = a[s];
        const int64_t j = b[s];
        const int64_t xi = x[i];
        const int64_t xj = x[j];
        const int64_t tot = xi + xj;
        const int64_t lo = tot >> 1;
        const int64_t hi = tot - lo;
        const int64_t ni = coin[t] ? hi : lo;
        const int64_t nj = tot - ni;
        const int64_t d = xi - xj;
        if (d == 0) { ++waits; continue; }
        if (d == 1 || d == -1) {
            if (ni == xi) { ++waits; continue; }
            ++swaps;
        } else {
            const int64_t ad = d < 0 ? -d : d;
            const int64_t change = ni * ni + nj * nj - xi * xi - xj * xj;
            if (change != -((ad * ad - (ad & 1)) / 2)) ++violations;
            transfers += ad >> 1;
            ++descents;
            vv += change;
        }
        x[i] = ni;
        x[j] = nj;
        if (first < 0 && vv == v_min) first = t;
    }
    counts[0] += descents;
    counts[1] += swaps;
    counts[2] += waits;
    counts[3] += transfers;
    counts[4] += violations;
    v[0] = vv;
    return first;
}
"""

_NATIVE: dict[str, Any] = {"lib": None, "tried": False, "error": None}


def _kernel_dir() -> Path:
    root = os.environ.get("OPH_EXACT_KERNEL_DIR")
    return Path(root) if root else Path(tempfile.gettempdir()) / "oph_exact_kernels"


def native_kernel() -> Any | None:
    """Compile (once) and load the native kernel; ``None`` when unavailable.

    Compiled with ``-ffp-contract=off`` so every floating operation matches
    the Python and numpy engines bit for bit.
    """

    if _NATIVE["tried"]:
        return _NATIVE["lib"]
    _NATIVE["tried"] = True
    if os.environ.get("OPH_EXACT_DISABLE_NATIVE"):
        _NATIVE["error"] = "disabled by OPH_EXACT_DISABLE_NATIVE"
        return None
    try:
        digest = hashlib.sha256(_KERNEL_SOURCE.encode("utf-8")).hexdigest()[:16]
        directory = _kernel_dir()
        directory.mkdir(parents=True, exist_ok=True)
        suffix = ".dll" if sys.platform == "win32" else ".so"
        library = directory / f"oph_fed_kernel_{digest}{suffix}"
        if not library.exists():
            compiler = shutil.which("cc") or shutil.which("clang") or shutil.which("gcc")
            if compiler is None:
                _NATIVE["error"] = "no C compiler on PATH"
                return None
            source = directory / f"oph_fed_kernel_{digest}.c"
            source.write_text(_KERNEL_SOURCE, encoding="utf-8")
            staging = directory / f"oph_fed_kernel_{digest}.{os.getpid()}{suffix}"
            subprocess.run(
                [compiler, "-O2", "-ffp-contract=off", "-shared", "-fPIC", "-o", str(staging), str(source)],
                check=True,
                capture_output=True,
            )
            os.replace(staging, library)
        lib = ctypes.CDLL(str(library))
        i64p = ctypes.POINTER(ctypes.c_int64)
        f64p = ctypes.POINTER(ctypes.c_double)
        lib.oph_fed_mean.argtypes = [f64p, i64p, i64p, i64p, i64p, i64p, ctypes.c_int64, f64p, i64p, f64p]
        lib.oph_fed_mean.restype = None
        lib.oph_fed_integer.argtypes = [i64p, i64p, i64p, i64p, i64p, ctypes.c_int64, i64p, i64p, ctypes.c_int64]
        lib.oph_fed_integer.restype = ctypes.c_int64
        _NATIVE["lib"] = lib
        return lib
    except Exception as error:  # pragma: no cover - environment dependent
        _NATIVE["error"] = f"{type(error).__name__}: {error}"
        return None


def resolve_engine(engine: str, ports: int | None = None) -> str:
    """``auto``: native when available, else the sequential loop on small
    federations (layering has no benefit below a few thousand ports) and the
    layered numpy engine above."""

    if engine == "auto":
        if native_kernel() is not None:
            return "native"
        return "python" if (ports is not None and ports < 4096) else "numpy"
    if engine == "native" and native_kernel() is None:
        raise RuntimeError(f"native kernel unavailable: {_NATIVE['error']}")
    if engine not in ("native", "numpy", "python"):
        raise ValueError(f"unknown engine {engine!r}")
    return engine


# --------------------------------------------------------------------------
# Move engines: mean law (float)
# --------------------------------------------------------------------------


@dataclass
class MeanTally:
    """Counters of one float mean-law engine call (sequential semantics)."""

    waits: int = 0
    ledger: float = 0.0
    descent_violations: int = 0
    phi_raising_moves: int = 0
    first_raise: dict[str, Any] | None = None


def _seam_sum_change(d: float, xi: float, xj: float, ni: float, nj: float, degi: float, degj: float) -> float:
    """``Phi(E_e x) - Phi(x)`` for seam ``(i, j)`` with neighbour sums ``ni, nj``.

    Identical operation order in every engine: ``-d ((deg_i x_i - N_i) -
    (deg_j x_j - N_j)) + d^2 (deg_i + deg_j + 2) / 4``.
    """

    return -d * ((degi * xi - ni) - (degj * xj - nj)) + (d * d) * ((degi + degj + 2.0) * 0.25)


def mean_moves_python(x: np.ndarray, graph: "PortGraph", seq: np.ndarray, tally: MeanTally) -> None:
    """Reference sequential engine."""

    xl = x.tolist()
    al = graph.a.tolist()
    bl = graph.b.tolist()
    nbrs = graph.neighbour_lists
    waits = 0
    ledger = 0.0
    violations = 0
    phi_up = 0
    for t, s in enumerate(seq.tolist()):
        i = al[s]
        j = bl[s]
        xi = xl[i]
        xj = xl[j]
        d = xi - xj
        if d == 0.0:
            waits += 1
            continue
        half_sq = 0.5 * (d * d)
        if not half_sq > 0.0:
            violations += 1
        ledger += half_sq
        ni = 0.0
        for k in nbrs[i]:
            ni += xl[k]
        nj = 0.0
        for k in nbrs[j]:
            nj += xl[k]
        degi = float(len(nbrs[i]))
        degj = float(len(nbrs[j]))
        delta = _seam_sum_change(d, xi, xj, ni, nj, degi, degj)
        if delta > 0.0:
            if phi_up == 0 and tally.first_raise is None:
                tally.first_raise = {"index": t, "seam": s, "x_i": xi, "x_j": xj, "N_i": ni, "N_j": nj, "deg_i": degi, "deg_j": degj, "delta_phi": delta}
            phi_up += 1
        m = 0.5 * (xi + xj)
        xl[i] = m
        xl[j] = m
    x[:] = xl
    tally.waits += waits
    tally.ledger += ledger
    tally.descent_violations += violations
    tally.phi_raising_moves += phi_up


def _dependency_layers(ea: np.ndarray, eb: np.ndarray, ports: int, table: np.ndarray | None = None):
    """Yield index arrays of pairwise-disjoint seams in dependency order.

    A seam is ready when it is the earliest unapplied seam at both of its
    ports; ready seams are pairwise disjoint and every earlier conflicting
    seam has been applied, so applying them together equals the sequential
    order exactly (bit for bit, since disjoint moves touch disjoint memory).
    With a neighbour ``table`` the criterion extends to the closed
    neighbourhoods of both ports (no earlier unapplied seam touches a
    neighbour), so every neighbour sum read from the pre-layer state is the
    sequential one.
    """

    remaining = np.arange(ea.size, dtype=np.int64)
    while remaining.size:
        ra = ea[remaining]
        rb = eb[remaining]
        r = remaining.size
        positions = np.arange(r, dtype=np.int64)
        ports_cat = np.concatenate([ra, rb])
        pos_cat = np.concatenate([positions, positions])
        order = np.lexsort((pos_cat, ports_cat))
        sorted_ports = ports_cat[order]
        starts = np.r_[True, sorted_ports[1:] != sorted_ports[:-1]]
        first = np.full(ports, r, dtype=np.int64)
        first[sorted_ports[starts]] = pos_cat[order][starts]
        ready = (first[ra] == positions) & (first[rb] == positions)
        if table is not None:
            for endpoints in (ra, rb):
                rows = table[endpoints]
                valid = rows >= 0
                f = first[np.where(valid, rows, 0)]
                ready &= np.all(~valid | (f >= positions[:, None]), axis=1)
        yield remaining[ready]
        remaining = remaining[~ready]


def _chunk_size(ports: int, chunk: int | None) -> int:
    if chunk is not None:
        return int(chunk)
    return int(min(16384, max(64, ports // 4)))


def _neighbour_sums(x: np.ndarray, table: np.ndarray, ports_of_moves: np.ndarray) -> np.ndarray:
    """Neighbour sums in neighbour-table (CSR) order, slot by slot, as the sequential engines."""

    rows = table[ports_of_moves]
    valid = rows >= 0
    values = np.where(valid, x[np.where(valid, rows, 0)], 0.0)
    total = np.zeros(ports_of_moves.size)
    for k in range(rows.shape[1]):
        total = total + values[:, k]
    return total


def mean_moves_numpy(x: np.ndarray, graph: "PortGraph", seq: np.ndarray, tally: MeanTally, chunk: int | None = None) -> None:
    """Layered vectorized engine; bit-identical to the reference engine."""

    ports = x.size
    a = graph.a
    b = graph.b
    degree = graph.degree.astype(float)
    table = graph.neighbour_table
    chunk = _chunk_size(ports, chunk)
    for start in range(0, seq.size, chunk):
        ids = seq[start : start + chunk]
        ea = a[ids]
        eb = b[ids]
        for layer in _dependency_layers(ea, eb, ports, table):
            la = ea[layer]
            lb = eb[layer]
            xa = x[la]
            xb = x[lb]
            d = xa - xb
            active = d != 0.0
            tally.waits += int(np.count_nonzero(~active))
            half_sq = 0.5 * (d * d)
            tally.descent_violations += int(np.count_nonzero(active & ~(half_sq > 0.0)))
            tally.ledger += float(np.sum(half_sq[active]))
            m = 0.5 * (xa + xb)
            ni = _neighbour_sums(x, table, la)
            nj = _neighbour_sums(x, table, lb)
            delta = _seam_sum_change(d, xa, xb, ni, nj, degree[la], degree[lb])
            raising = active & (delta > 0.0)
            if np.any(raising):
                if tally.phi_raising_moves == 0 and tally.first_raise is None:
                    k = int(np.flatnonzero(raising)[np.argmin(layer[raising])])
                    tally.first_raise = {
                        "index": int(start + layer[k]),
                        "seam": int(ids[layer[k]]),
                        "x_i": float(xa[k]),
                        "x_j": float(xb[k]),
                        "N_i": float(ni[k]),
                        "N_j": float(nj[k]),
                        "deg_i": float(degree[la[k]]),
                        "deg_j": float(degree[lb[k]]),
                        "delta_phi": float(delta[k]),
                    }
                tally.phi_raising_moves += int(np.count_nonzero(raising))
            x[la[active]] = m[active]
            x[lb[active]] = m[active]


def mean_moves_native(x: np.ndarray, graph: "PortGraph", seq: np.ndarray, tally: MeanTally) -> None:
    lib = native_kernel()
    if lib is None:
        raise RuntimeError("native kernel unavailable")
    ledger = np.zeros(1, dtype=np.float64)
    counts = np.array([0, 0, 0, -1, -1], dtype=np.int64)
    witness = np.zeros(7, dtype=np.float64)
    seq = np.ascontiguousarray(seq, dtype=np.int64)
    i64p = ctypes.POINTER(ctypes.c_int64)
    f64p = ctypes.POINTER(ctypes.c_double)
    lib.oph_fed_mean(
        x.ctypes.data_as(f64p),
        graph.a.ctypes.data_as(i64p),
        graph.b.ctypes.data_as(i64p),
        graph.indptr.ctypes.data_as(i64p),
        graph.indices.ctypes.data_as(i64p),
        seq.ctypes.data_as(i64p),
        ctypes.c_int64(seq.size),
        ledger.ctypes.data_as(f64p),
        counts.ctypes.data_as(i64p),
        witness.ctypes.data_as(f64p),
    )
    if counts[2] > 0 and tally.phi_raising_moves == 0 and tally.first_raise is None:
        tally.first_raise = {
            "index": int(counts[3]),
            "seam": int(counts[4]),
            "x_i": float(witness[0]),
            "x_j": float(witness[1]),
            "N_i": float(witness[2]),
            "N_j": float(witness[3]),
            "deg_i": float(witness[4]),
            "deg_j": float(witness[5]),
            "delta_phi": float(witness[6]),
        }
    tally.waits += int(counts[0])
    tally.descent_violations += int(counts[1])
    tally.phi_raising_moves += int(counts[2])
    tally.ledger += float(ledger[0])


MEAN_ENGINES = {"python": mean_moves_python, "numpy": mean_moves_numpy, "native": mean_moves_native}


# --------------------------------------------------------------------------
# Move engines: integer law
# --------------------------------------------------------------------------


@dataclass
class IntegerTally:
    descents: int = 0
    swaps: int = 0
    waits: int = 0
    unit_transfers: int = 0
    decrement_violations: int = 0


def integer_moves_python(
    x: np.ndarray, graph: "PortGraph", seq: np.ndarray, coin: np.ndarray, v: int, v_min: int, tally: IntegerTally
) -> tuple[int, int]:
    """Reference engine.  Returns ``(v, first_index)``."""

    xl = x.tolist()
    al = graph.a.tolist()
    bl = graph.b.tolist()
    first = -1
    for t, (s, c) in enumerate(zip(seq.tolist(), coin.tolist())):
        i = al[s]
        j = bl[s]
        xi = xl[i]
        xj = xl[j]
        ni, nj = carrier.integer_nearest_agreement(xi, xj, ceiling_to_first=bool(c))
        d = xi - xj
        if d == 0:
            tally.waits += 1
            continue
        if d == 1 or d == -1:
            if ni == xi:
                tally.waits += 1
                continue
            tally.swaps += 1
        else:
            ad = abs(d)
            change = ni * ni + nj * nj - xi * xi - xj * xj
            if change != -((ad * ad - (ad & 1)) // 2):
                tally.decrement_violations += 1
            tally.unit_transfers += ad >> 1
            tally.descents += 1
            v += change
        xl[i] = ni
        xl[j] = nj
        if first < 0 and v == v_min:
            first = t
    x[:] = xl
    return int(v), first


def integer_moves_numpy(
    x: np.ndarray, graph: "PortGraph", seq: np.ndarray, coin: np.ndarray, v: int, v_min: int, tally: IntegerTally, chunk: int | None = None
) -> tuple[int, int]:
    first = -1
    ports = x.size
    a = graph.a
    b = graph.b
    chunk = _chunk_size(ports, chunk)
    for start in range(0, seq.size, chunk):
        ids = seq[start : start + chunk]
        ea = a[ids]
        eb = b[ids]
        ec = coin[start : start + chunk].astype(bool)
        delta = np.zeros(ids.size, dtype=np.int64)
        for layer in _dependency_layers(ea, eb, ports):
            la = ea[layer]
            lb = eb[layer]
            xa = x[la]
            xb = x[lb]
            tot = xa + xb
            lo = tot >> 1
            hi = tot - lo
            na = np.where(ec[layer], hi, lo)
            nb = tot - na
            d = xa - xb
            ad = np.abs(d)
            unit = ad == 1
            wait = (d == 0) | (unit & (na == xa))
            swap = unit & ~wait
            desc = ~unit & (d != 0)
            change = na * na + nb * nb - xa * xa - xb * xb
            tally.waits += int(np.count_nonzero(wait))
            tally.swaps += int(np.count_nonzero(swap))
            tally.descents += int(np.count_nonzero(desc))
            tally.unit_transfers += int(np.sum((ad >> 1)[desc]))
            tally.decrement_violations += int(np.count_nonzero(desc & (change != -((ad * ad - (ad & 1)) // 2))))
            move = ~wait
            x[la[move]] = na[move]
            x[lb[move]] = nb[move]
            delta[layer[desc]] = change[desc]
        if first < 0:
            running = v + np.cumsum(delta)
            hit = np.flatnonzero(running == v_min)
            if hit.size:
                first = start + int(hit[0])
        v += int(delta.sum())
    return int(v), first


def integer_moves_native(
    x: np.ndarray, graph: "PortGraph", seq: np.ndarray, coin: np.ndarray, v: int, v_min: int, tally: IntegerTally
) -> tuple[int, int]:
    lib = native_kernel()
    if lib is None:
        raise RuntimeError("native kernel unavailable")
    counts = np.zeros(5, dtype=np.int64)
    vv = np.array([v], dtype=np.int64)
    seq = np.ascontiguousarray(seq, dtype=np.int64)
    coin = np.ascontiguousarray(coin, dtype=np.int64)
    i64p = ctypes.POINTER(ctypes.c_int64)
    first = lib.oph_fed_integer(
        x.ctypes.data_as(i64p),
        graph.a.ctypes.data_as(i64p),
        graph.b.ctypes.data_as(i64p),
        seq.ctypes.data_as(i64p),
        coin.ctypes.data_as(i64p),
        ctypes.c_int64(seq.size),
        counts.ctypes.data_as(i64p),
        vv.ctypes.data_as(i64p),
        ctypes.c_int64(v_min),
    )
    tally.descents += int(counts[0])
    tally.swaps += int(counts[1])
    tally.waits += int(counts[2])
    tally.unit_transfers += int(counts[3])
    tally.decrement_violations += int(counts[4])
    return int(vv[0]), int(first)


INTEGER_ENGINES = {"python": integer_moves_python, "numpy": integer_moves_numpy, "native": integer_moves_native}


# --------------------------------------------------------------------------
# Schedules
# --------------------------------------------------------------------------


def _trace_sweep(sweep: int) -> bool:
    return sweep == 0 or (sweep & (sweep - 1)) == 0


def run_mean_law(
    fed: ExactFederation,
    loads: np.ndarray,
    seed: int,
    *,
    engine: str = "auto",
    threshold: float = FLOAT_PHI_THRESHOLD,
    max_sweeps: int = MAX_SWEEPS_DEFAULT,
) -> dict[str, Any]:
    """One asynchronous float schedule of the mean law to ``Phi < threshold``."""

    engine = resolve_engine(engine, fed.ports)
    move = MEAN_ENGINES[engine]
    graph = fed.graph
    x = np.asarray(loads, dtype=np.float64).copy()
    rng = np.random.default_rng(seed)
    seams = fed.seams
    phi0 = fed.mismatch_potential(x)
    v0 = fed.descent_functional(x)
    v_min = fed.mean_minimum_descent(loads)
    trace = [[0, _sig(phi0, PHI_TRACE_ROUND)]]
    tally = MeanTally()
    sweep = 0
    phi = phi0
    terminated = phi < threshold
    first_raise_sweep = None
    while not terminated and sweep < max_sweeps:
        seq = rng.integers(0, seams, size=seams, dtype=np.int64)
        before = tally.first_raise is None
        move(x, graph, seq, tally)
        if before and tally.first_raise is not None:
            first_raise_sweep = sweep
        sweep += 1
        phi = fed.mismatch_potential(x)
        if _trace_sweep(sweep):
            trace.append([sweep, _sig(phi, PHI_TRACE_ROUND)])
        terminated = phi < threshold
    if trace[-1][0] != sweep:
        trace.append([sweep, _sig(phi, PHI_TRACE_ROUND)])
    v_end = fed.descent_functional(x)
    ledger_error = abs((v0 - v_end) - tally.ledger) / max(v0, 1.0)
    means = fed.component_means(loads)
    residual_vector = x - means[fed.component_of_port]
    deviation = float(np.max(np.abs(residual_vector)))
    centered_terminal = float(np.dot(residual_vector, residual_vector))
    quotient_hash, residual = fed.terminal_quotient_hash(x)
    # The snap is a certificate only when it is unambiguous; an unsettled state gets no hash.
    snap_ok = bool(residual < LATTICE_SNAP_MARGIN)
    counterexample = None
    if tally.first_raise is not None:
        fr = tally.first_raise
        counterexample = {
            "attempt": int(first_raise_sweep * seams + fr["index"] + 1),
            "seam": int(fr["seam"]),
            "ports": [int(fed.seam_a[fr["seam"]]), int(fed.seam_b[fr["seam"]])],
            "x_i": _round(fr["x_i"]),
            "x_j": _round(fr["x_j"]),
            "neighbour_sum_i": _round(fr["N_i"]),
            "neighbour_sum_j": _round(fr["N_j"]),
            "deg_i": int(fr["deg_i"]),
            "deg_j": int(fr["deg_j"]),
            "delta_phi": _round(fr["delta_phi"]),
        }
    attempts = sweep * seams
    return {
        "seed": int(seed),
        "terminated": bool(terminated),
        "sweeps": int(sweep),
        "attempts": int(attempts),
        "waits": int(tally.waits),
        "non_wait_moves": int(attempts - tally.waits),
        "phi_initial": _sig(phi0, PHI_TRACE_ROUND),
        "phi_terminal": _sig(phi, PHI_TRACE_ROUND),
        "phi_trace": trace,
        "descent_functional_initial": _round(v0),
        "descent_functional_terminal": _round(v_end),
        "descent_functional_minimum": _round(v_min),
        "centered_norm_initial": _round(v0 - v_min),
        "centered_norm_terminal": _sig(centered_terminal, 6),
        "descent_ledger_relative_error_below_1e-9": bool(ledger_error < 1e-9),
        "strict_descent_violations": int(tally.descent_violations),
        "phi_raising_moves": int(tally.phi_raising_moves),
        "phi_raising_fraction_of_non_wait_moves": _sig(tally.phi_raising_moves / max(1, attempts - tally.waits), 4),
        "phi_raising_counterexample": counterexample,
        "max_abs_deviation_from_component_mean": _sig(deviation, 3),
        "component_mean_within_1e-9": bool(deviation < 1e-9),
        "lattice_residual_max": _sig(residual, 3),
        "lattice_snap_unambiguous": snap_ok,
        "terminal_quotient_hash": quotient_hash if snap_ok else None,
        "state": x,
    }


def run_integer_law(
    fed: ExactFederation,
    loads: np.ndarray,
    seed: int,
    *,
    engine: str = "auto",
    max_sweeps: int = MAX_SWEEPS_DEFAULT,
) -> dict[str, Any]:
    """One asynchronous integer-law schedule to the balanced class."""

    engine = resolve_engine(engine, fed.ports)
    move = INTEGER_ENGINES[engine]
    graph = fed.graph
    x = np.asarray(loads, dtype=np.int64).copy()
    if np.any(x < 0):
        raise ValueError("integer loads must be nonnegative")
    rng = np.random.default_rng(seed)
    seams = fed.seams
    v = int(np.dot(x, x))
    v0 = v
    v_min = fed.integer_minimum_descent(loads)
    tally = IntegerTally()
    sweep = 0
    first_attempt = 0 if v == v_min else -1
    while first_attempt < 0 and sweep < max_sweeps:
        seq = rng.integers(0, seams, size=seams, dtype=np.int64)
        coin = rng.integers(0, 2, size=seams, dtype=np.int64)
        v, first = move(x, graph, seq, coin, v, v_min, tally)
        if first >= 0:
            first_attempt = sweep * seams + first + 1
        sweep += 1
    v_state = int(np.dot(x, x))
    d = np.abs(x[fed.seam_a] - x[fed.seam_b])
    return {
        "seed": int(seed),
        "terminated": bool(first_attempt >= 0),
        "sweeps": int(sweep),
        "attempts_to_balanced_class": int(first_attempt),
        "descents": int(tally.descents),
        "swaps": int(tally.swaps),
        "waits": int(tally.waits),
        "unit_transfers": int(tally.unit_transfers),
        "unit_transfer_decrement_identity_violations": int(tally.decrement_violations),
        "descent_functional_initial": int(v0),
        "descent_functional_terminal": int(v_state),
        "descent_functional_minimum": int(v_min),
        "centered_norm_initial": int(v0 - v_min),
        "descent_ledger_exact": bool(v == v_state),
        "strict_descent_violations": int(tally.decrement_violations),
        "odd_tie_seams_at_termination": int(np.count_nonzero(d == 1)),
        "max_seam_difference_at_termination": int(d.max()) if d.size else 0,
        "quotient_hash": fed.integer_quotient_hash(x),
        "state": x,
    }


def run_mean_law_exact(
    fed: ExactFederation,
    loads: np.ndarray,
    seed: int,
    *,
    engine: str = "auto",
    max_sweeps: int = MAX_SWEEPS_DEFAULT,
    descent_sample_every: int = 1,
) -> dict[str, Any]:
    """One schedule in exact dyadic arithmetic ``x_p = n_p / 2^K``.

    ``K`` grows on demand (the state is invariant under ``(n, K) -> (n << r,
    K + r)``).  The float companion runs the same moves in IEEE arithmetic and
    triggers the exact ``Phi`` evaluation once it is below
    ``EXACT_TRIGGER_FACTOR`` times the threshold; the termination decision
    ``Phi <= 10^-18`` is exact.  The per-move descent identity
    ``V(x) - V(E_e x) = (x_i - x_j)^2 / 2`` is checked exactly on every
    ``descent_sample_every``-th non-wait move (every move when the sample
    step is one, in which case the global identity is checked as well).
    """

    engine = resolve_engine(engine, fed.ports)
    move = MEAN_ENGINES[engine]
    graph = fed.graph
    companion = MeanTally()
    al = graph.a.tolist()
    bl = graph.b.tolist()
    seams = fed.seams
    rescale_bits = 1024
    K = rescale_bits
    rescales = 0
    n = [int(v) << K for v in np.asarray(loads, dtype=np.int64).tolist()]
    xf = np.asarray(loads, dtype=np.float64).copy()
    rng = np.random.default_rng(seed)
    total0 = sum(int(v) for v in np.asarray(loads).tolist())
    v0 = sum(v * v for v in n)  # units 2^(-2K) at the initial K
    v0_bits = 2 * K
    ledger = 0  # sum of (n_i - n_j)^2 in units 2^(-2K) at the current K (full ledger only)
    full_ledger = descent_sample_every == 1
    step = max(1, int(descent_sample_every))
    waits = 0
    non_wait = 0
    identity_checked = 0
    identity_ok = 0
    sweep = 0
    exact_checks = 0
    threshold_ok = False
    phi_exact_num = None
    ratio_min = None
    ratio_max = None
    trigger = EXACT_TRIGGER_FACTOR * FLOAT_PHI_THRESHOLD

    def exact_phi_numerator() -> int:
        return sum((n[i] - n[j]) ** 2 for i, j in zip(al, bl))

    def evaluate_exact(phi_float: float) -> bool:
        nonlocal phi_exact_num, exact_checks, ratio_min, ratio_max
        phi_exact_num = exact_phi_numerator()
        exact_checks += 1
        if phi_exact_num > 0:
            ratio = phi_float / float(Fraction(phi_exact_num, 1 << (2 * K)))
            ratio_min = ratio if ratio_min is None else min(ratio_min, ratio)
            ratio_max = ratio if ratio_max is None else max(ratio_max, ratio)
        return phi_exact_num * EXACT_PHI_THRESHOLD_DENOMINATOR <= (1 << (2 * K))

    phi_f = fed.mismatch_potential(xf)
    if phi_f <= trigger:
        threshold_ok = evaluate_exact(phi_f)
    float_first_sweep = 0 if phi_f < FLOAT_PHI_THRESHOLD else -1
    while not threshold_ok and sweep < max_sweeps:
        seq = rng.integers(0, seams, size=seams, dtype=np.int64)
        for s in seq.tolist():
            i = al[s]
            j = bl[s]
            ni = n[i]
            nj = n[j]
            if ni == nj:
                waits += 1
                continue
            tot = ni + nj
            if tot & 1:
                n = [v << rescale_bits for v in n]
                ledger <<= 2 * rescale_bits
                K += rescale_bits
                rescales += 1
                ni = n[i]
                nj = n[j]
                tot = ni + nj
            m = tot >> 1
            non_wait += 1
            if full_ledger:
                dd = ni - nj
                ledger += dd * dd
            elif non_wait % step == 0:
                dd = ni - nj
                identity_checked += 1
                identity_ok += int(2 * (ni * ni + nj * nj - 2 * m * m) == dd * dd)
            n[i] = m
            n[j] = m
        move(xf, graph, seq, companion)
        sweep += 1
        phi_f = fed.mismatch_potential(xf)
        if float_first_sweep < 0 and phi_f < FLOAT_PHI_THRESHOLD:
            float_first_sweep = sweep
        if phi_f <= trigger:
            threshold_ok = evaluate_exact(phi_f)
    if phi_exact_num is None:
        evaluate_exact(phi_f)
    scale_sq = 1 << (2 * K)
    v_end = sum(v * v for v in n)
    global_identity = ((v0 << (2 * K - v0_bits)) - v_end) * 2 == ledger if full_ledger else None
    if full_ledger:
        identity_checked = non_wait
        identity_ok = non_wait if global_identity else 0
    conservation = sum(n) == total0 << K
    comp = fed.component_of_port.tolist()
    comp_totals = [0] * fed.components
    for p, c in enumerate(comp):
        comp_totals[c] += n[p]
    load_totals = np.bincount(fed.component_of_port, weights=np.asarray(loads, dtype=float), minlength=fed.components)
    component_conservation = all(comp_totals[c] == int(round(float(load_totals[c]))) << K for c in range(fed.components))
    sizes = fed.component_sizes.tolist()
    deviation = Fraction(0)
    residual_ok = True
    snapped = []
    for p, c in enumerate(comp):
        m_c = sizes[c]
        s_c = int(round(float(load_totals[c])))
        dev = Fraction(abs(m_c * n[p] - (s_c << K)), m_c << K)
        if dev > deviation:
            deviation = dev
        scaled = Fraction(m_c * n[p], 1 << K)
        q = round(scaled)
        if abs(scaled - q) >= Fraction(1, 2):
            residual_ok = False
        snapped.append(int(q))
    quotient_hash = sha256_of(
        {
            "canonicalizer": "component_lattice_snap",
            "component_sizes": sizes,
            "component_of_port": comp,
            "q": snapped,
        }
    )
    denominator = 1 << K
    x_exact_float = np.array([v / denominator for v in n], dtype=np.float64)
    float_vs_exact = float(np.max(np.abs(xf - x_exact_float)))
    phi_exact = float(Fraction(phi_exact_num, scale_sq))
    return {
        "seed": int(seed),
        "precision_bits": K,
        "precision_rescales": rescales,
        "terminated": bool(threshold_ok),
        "sweeps": int(sweep),
        "attempts": int(sweep * seams),
        "waits": int(waits),
        "non_wait_moves": int(non_wait),
        "exact_phi_evaluations": int(exact_checks),
        "float_over_exact_phi_ratio_range": [
            _sig(ratio_min, 6) if ratio_min is not None else None,
            _sig(ratio_max, 6) if ratio_max is not None else None,
        ],
        "phi_terminal": _sig(phi_exact, 6),
        "phi_terminal_at_most_1e-18_exact": bool(threshold_ok),
        "descent_identity_sample_step": step,
        "descent_identity_moves_checked": int(identity_checked),
        "descent_identity_exact": bool(identity_ok == identity_checked),
        "global_descent_identity_exact": global_identity,
        "total_conservation_exact": bool(conservation),
        "component_conservation_exact": bool(component_conservation),
        "strict_descent_violations": 0,
        "max_abs_deviation_from_component_mean": _sig(float(deviation), 3),
        "lattice_snap_unambiguous": bool(residual_ok),
        "terminal_quotient_hash": quotient_hash,
        "float_companion_waits": int(companion.waits),
        "float_companion_first_sweep_below_threshold": int(float_first_sweep),
        "float_companion_phi_at_termination": _sig(phi_f, 6),
        "float_vs_exact_max_abs_deviation": _sig(float_vs_exact, 3),
    }


# --------------------------------------------------------------------------
# Synchronous operator and response kernels
# --------------------------------------------------------------------------


def synchronous_operator_receipt(fed: ExactFederation) -> dict[str, Any]:
    lap = fed.laplacian
    D = fed.mean_denominator
    block_ok = None
    if fed.gluing == "isolated":
        expected = carrier.repair_mean_exact()
        block_ok = True
        for c in (0, fed.carriers // 2, fed.carriers - 1):
            block = lap[c * PORTS : (c + 1) * PORTS, c * PORTS : (c + 1) * PORTS].toarray()
            for r in range(PORTS):
                for s in range(PORTS):
                    entry = Fraction(int(r == s)) - Fraction(int(round(block[r, s]))) / D
                    if entry != expected[r][s]:
                        block_ok = False
    ports = fed.ports
    if fed.gluing == "isolated":
        eig = np.linalg.eigvalsh(carrier.laplacian().astype(float))
        lam2 = float(eig[1])
        lam_max = float(eig[-1])
        method = "single_carrier_block_eigvalsh"
    elif ports <= 4096:
        eig = np.linalg.eigvalsh(lap.toarray())
        lam2 = float(eig[fed.components]) if fed.components < ports else 0.0
        lam_max = float(eig[-1])
        method = "dense_eigvalsh"
    else:
        lam2 = float(
            np.sort(sparse_linalg.eigsh(lap, k=fed.components + 1, sigma=-1e-3, which="LM", return_eigenvectors=False))[-1]
        )
        lam_max = float(sparse_linalg.eigsh(lap, k=1, which="LA", return_eigenvectors=False)[0])
        method = "sparse_shift_invert"
    return {
        "definition": "T_fed = I - L_fed/D with D = 2|S|/N (one attempt per carrier per tick); per attempt the expectation operator is I - L_fed/(2|S|)",
        "D": str(D),
        "per_attempt_denominator": int(2 * fed.seams),
        "isolated_block_equals_carrier_repair_mean_exact": block_ok,
        "laplacian_lambda_2": _round(lam2, 10),
        "laplacian_lambda_max": _round(lam_max, 10),
        "spectrum_method": method,
        "t_fed_second_eigenvalue": _round(1.0 - lam2 / float(D), 12),
        "t_fed_smallest_eigenvalue": _round(1.0 - lam_max / float(D), 12),
        "attempts_per_e_fold_of_slowest_mode": _sig(2 * fed.seams / lam2, 6) if lam2 > 0 else None,
    }


def sample_cells(fed: ExactFederation, count: int = KERNEL_SAMPLE) -> list[int]:
    pent = np.flatnonzero(fed.pentagonal_cell)
    interior = np.flatnonzero(~fed.pentagonal_cell)
    chosen: list[int] = []
    half = max(4, count // 3)
    chosen.extend(int(v) for v in pent[:half])
    chosen.extend(int(v) for v in interior[:half])
    rng = np.random.default_rng(LOAD_SEED_BASE + 17 * fed.level)
    pool = [c for c in range(fed.carriers) if c not in chosen]
    while len(chosen) < count and pool:
        pick = int(rng.integers(0, len(pool)))
        chosen.append(pool.pop(pick))
    return sorted(set(chosen))


def response_kernels(fed: ExactFederation, cell: int, steps: Sequence[int] = KERNEL_STEPS) -> dict[int, np.ndarray]:
    """``K_n = 12 C_n / tr C_n`` for one carrier from centered impulse propagation.

    The twelve centered impulses ``e_{c,p} - u_c/12`` are propagated by
    ``T_fed`` (``2n`` applications, common rescaling each step); reading back
    carrier ``c`` gives ``R_n Q`` and ``C_n = Q R_n Q``.
    """

    ports = fed.ports
    y = np.zeros((ports, PORTS))
    block = slice(cell * PORTS, (cell + 1) * PORTS)
    y[block, :] = np.eye(PORTS) - np.ones((PORTS, PORTS)) / PORTS
    q = np.eye(PORTS) - np.ones((PORTS, PORTS)) / PORTS
    out: dict[int, np.ndarray] = {}
    target = {2 * n: n for n in steps}
    last = max(target)
    for step in range(1, last + 1):
        y = fed.apply_synchronous(y)
        scale = np.max(np.abs(y))
        if scale > 0:
            y /= scale
        if step in target:
            c = q @ y[block, :]
            c = 0.5 * (c + c.T)
            out[target[step]] = PORTS * c / np.trace(c)
    return out


def green_pattern_slow_share_reference() -> float:
    """Slow-band share of an isotropic boundary Green's pattern on one carrier.

    For inflow through glued ports with covariance proportional to the identity
    the late-time centered kernel is spanned by columns of the carrier Green's
    function ``G = sum_b P_b / mu_b`` over the nonconstant bands, and the share
    is ``(3/mu_s^2) / (3/mu_s^2 + 5/mu_m^2 + 3/mu_f^2)`` by vertex transitivity.
    """

    weights = {}
    for key, (mu, mult) in zip(("slow", "middle", "fast"), carrier.LAPLACIAN_BANDS[1:]):
        weights[key] = mult / carrier.q5_float(mu) ** 2
    return weights["slow"] / sum(weights.values())


def kernel_receipt(fed: ExactFederation, cells: Sequence[int] | None = None) -> dict[str, Any]:
    cells = list(cells) if cells is not None else sample_cells(fed)
    p_slow = carrier.slow_band_projector()
    gram = carrier.intrinsic_gram()
    per_cell = []
    isolated_dev = 0.0
    gram_dev_300 = 0.0
    shares: dict[int, list[float]] = {n: [] for n in KERNEL_STEPS}
    tops: dict[int, list[list[float]]] = {n: [] for n in KERNEL_STEPS}
    for cell in cells:
        kernels = response_kernels(fed, cell)
        rows = []
        for n in KERNEL_STEPS:
            k = kernels[n]
            eig = np.sort(np.linalg.eigvalsh(k))[::-1]
            share = float(np.trace(p_slow @ k @ p_slow) / np.trace(k))
            rank = int(np.linalg.matrix_rank(k, tol=1e-9))
            if fed.gluing == "isolated":
                isolated_dev = max(isolated_dev, float(np.max(np.abs(k - carrier.normalized_response_kernel(n)))))
            if n == 300:
                gram_dev_300 = max(gram_dev_300, float(np.max(np.abs(k - gram))))
            shares[n].append(share)
            tops[n].append([float(v) for v in eig[:4]])
            rows.append(
                {
                    "n": n,
                    "top_eigenvalues": [_round(v, 9) for v in eig[:4]],
                    "eigenvalues": [_round(v, 9) for v in eig],
                    "slow_band_share": _round(share, 9),
                    "rank_1e-9": rank,
                }
            )
        per_cell.append({"cell": int(cell), "pentagonal": bool(fed.pentagonal_cell[cell]), "per_step": rows})
    summary = {}
    for n in KERNEL_STEPS:
        s = np.asarray(shares[n])
        t = np.asarray(tops[n])
        summary[str(n)] = {
            "slow_band_share": {"min": _round(s.min(), 9), "median": _round(float(np.median(s)), 9), "max": _round(s.max(), 9)},
            "top_four_eigenvalues_min": [_round(v, 9) for v in t.min(axis=0)],
            "top_four_eigenvalues_median": [_round(v, 9) for v in np.median(t, axis=0)],
            "top_four_eigenvalues_max": [_round(v, 9) for v in t.max(axis=0)],
        }
    return {
        "definition": (
            "probe port p of carrier c with a unit impulse, apply T_fed 2n times, read carrier c: column p of R_n; "
            "C_n = Q R_n Q with Q = I - J/12; K_n = 12 C_n / tr C_n; computed by propagating the centered impulses "
            "R_n Q directly with a common rescaling per step"
        ),
        "steps": list(KERNEL_STEPS),
        "sampled_cells": [int(c) for c in cells],
        "pentagonal_cells_in_sample": int(sum(bool(fed.pentagonal_cell[c]) for c in cells)),
        "slow_band_projector": "carrier.slow_band_projector() (5 - sqrt5 band, trace 3)",
        "slow_band_share_definition": "tr(P_slow K_n P_slow) / tr(K_n)",
        "green_pattern_slow_share_reference": _round(green_pattern_slow_share_reference(), 9),
        "green_pattern_reading": (
            "reference for the late-time glued kernel: with every centered impulse draining through the glued ports, "
            "the readback converges to a pattern spanned by the carrier Green's columns of those ports; for an isotropic "
            "inflow covariance its slow-band share is the reference value and its rank is at most the number of glued ports"
        ),
        "isolated_max_abs_deviation_from_carrier_kernel": _sig(isolated_dev, 3) if fed.gluing == "isolated" else None,
        "isolated_kernel_matches_carrier_within_1e-12": bool(isolated_dev < 1e-12) if fed.gluing == "isolated" else None,
        "max_abs_deviation_from_intrinsic_gram_at_n_300": _sig(gram_dev_300, 3),
        "per_cell": per_cell,
        "summary": summary,
    }


# --------------------------------------------------------------------------
# Exact local confluence data (small rungs)
# --------------------------------------------------------------------------


def _apply_mean_exact(x: list[Fraction], i: int, j: int) -> list[Fraction]:
    y = list(x)
    m = (y[i] + y[j]) / 2
    y[i] = m
    y[j] = m
    return y


def _component_projection_exact(fed: ExactFederation, x: list[Fraction]) -> list[Fraction]:
    comp = fed.component_of_port.tolist()
    totals = [Fraction(0)] * fed.components
    for p, c in enumerate(comp):
        totals[c] += x[p]
    sizes = fed.component_sizes.tolist()
    return [totals[c] / sizes[c] for c in comp]


def local_confluence_receipt(fed: ExactFederation, loads: np.ndarray, *, pairs: int = 64) -> dict[str, Any]:
    rng = np.random.default_rng(LOAD_SEED_BASE + 31 * fed.level + GLUINGS.index(fed.gluing))
    x = [Fraction(int(v)) for v in np.asarray(loads).tolist()]
    a = fed.seam_a.tolist()
    b = fed.seam_b.tolist()
    seams = fed.seams
    disjoint_checked = 0
    disjoint_commute = 0
    conflict_checked = 0
    conflict_one_step_diamond = 0
    conflict_aggregate_join = 0
    conflict_terminal_join = 0
    port_seams: dict[int, list[int]] = {}
    for e in range(seams):
        port_seams.setdefault(a[e], []).append(e)
        port_seams.setdefault(b[e], []).append(e)
    while disjoint_checked < pairs:
        e = int(rng.integers(0, seams))
        f = int(rng.integers(0, seams))
        if len({a[e], b[e], a[f], b[f]}) < 4:
            continue
        ef = _apply_mean_exact(_apply_mean_exact(x, a[e], b[e]), a[f], b[f])
        fe = _apply_mean_exact(_apply_mean_exact(x, a[f], b[f]), a[e], b[e])
        disjoint_checked += 1
        disjoint_commute += int(ef == fe)
    pi_x = _component_projection_exact(fed, x)
    while conflict_checked < pairs:
        e = int(rng.integers(0, seams))
        shared = a[e] if rng.integers(0, 2) else b[e]
        candidates = [f for f in port_seams[shared] if f != e]
        f = candidates[int(rng.integers(0, len(candidates)))]
        union = sorted({a[e], b[e], a[f], b[f]})
        if len(union) != 3:
            continue
        xe = _apply_mean_exact(x, a[e], b[e])
        xf = _apply_mean_exact(x, a[f], b[f])
        ef = _apply_mean_exact(xe, a[f], b[f])
        fe = _apply_mean_exact(xf, a[e], b[e])
        conflict_checked += 1
        conflict_one_step_diamond += int(ef == fe)

        def aggregate(y: list[Fraction]) -> list[Fraction]:
            z = list(y)
            m = sum(z[p] for p in union) / 3
            for p in union:
                z[p] = m
            return z

        conflict_aggregate_join += int(aggregate(xe) == aggregate(xf) == aggregate(x))
        conflict_terminal_join += int(_component_projection_exact(fed, ef) == _component_projection_exact(fed, fe) == pi_x)
    # Terminal-quotient invariance of every single seam move: Pi E_e x = Pi x.
    invariant = 0
    for e in range(seams):
        if _component_projection_exact(fed, _apply_mean_exact(x, a[e], b[e])) == pi_x:
            invariant += 1
    # Exact per-move descent identity on a random walk of moves.
    walk_moves = 500
    y = list(x)
    identity_ok = 0
    for _ in range(walk_moves):
        e = int(rng.integers(0, seams))
        i, j = a[e], b[e]
        d = y[i] - y[j]
        z = _apply_mean_exact(y, i, j)
        drop = sum(v * v for v in y) - sum(v * v for v in z)
        identity_ok += int(drop == d * d / 2)
        y = z
    return {
        "arithmetic": "fractions.Fraction",
        "disjoint_pairs_checked": disjoint_checked,
        "disjoint_pairs_commute_exactly": disjoint_commute,
        "conflicting_pairs_checked": conflict_checked,
        "conflicting_pairs_one_step_diamond": conflict_one_step_diamond,
        "conflicting_pairs_joined_by_aggregate_component_mean": conflict_aggregate_join,
        "conflicting_pairs_same_terminal_quotient": conflict_terminal_join,
        "single_moves_fixing_terminal_quotient": invariant,
        "seams": seams,
        "descent_identity_moves_checked": walk_moves,
        "descent_identity_exact": identity_ok,
        "reading": (
            "disjoint moves commute exactly; conflicting single moves have no one-step diamond in general, "
            "their conflict component is joined exactly by the aggregate three-port mean (the canonical "
            "aggregate payload of the transactional diamond), and every single move fixes the terminal "
            "quotient Pi x (component means), which is the exact reason every schedule shares one limit"
        ),
    }


def phi_nonmonotone_witness() -> dict[str, Any]:
    """An exact single-carrier witness that ``Phi`` rises along one seam move."""

    seams = carrier.seams()
    adj = {p: set() for p in range(PORTS)}
    for i, j in seams:
        adj[i].add(j)
        adj[j].add(i)
    a, b = seams[0]
    common = adj[a] & adj[b]
    # Phi(E_e x) - Phi(x) = -2 d^2 + d (N_a - N_b) with d = x_a - x_b and N_p the
    # neighbour sum of p, so private neighbours of a above x_a make Phi rise.
    x = [Fraction(1)] * PORTS
    x[a] = Fraction(2)
    x[b] = Fraction(0)
    for p in adj[a] - {b} - common:
        x[p] = Fraction(5)
    for p in adj[b] - {a} - common:
        x[p] = Fraction(0)
    y = _apply_mean_exact(x, a, b)

    def phi(z: list[Fraction]) -> Fraction:
        return sum((z[i] - z[j]) ** 2 for i, j in seams)

    def v(z: list[Fraction]) -> Fraction:
        return sum(t * t for t in z)

    return {
        "seam": [a, b],
        "state": [str(t) for t in x],
        "phi_before": str(phi(x)),
        "phi_after": str(phi(y)),
        "phi_rises": bool(phi(y) > phi(x)),
        "descent_functional_before": str(v(x)),
        "descent_functional_after": str(v(y)),
        "descent_drop_equals_half_square_difference": bool(v(x) - v(y) == (x[a] - x[b]) ** 2 / 2),
        "reading": "Phi is the termination potential, V = sum x_p^2 is the single-move descent functional",
    }


# --------------------------------------------------------------------------
# Receipt assembly
# --------------------------------------------------------------------------


def _strip(entry: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in entry.items() if k != "state"}


def exact_descent_sample_step(fed: ExactFederation) -> int:
    """Full ledger (every non-wait move) on L0 and on isolated rungs; every 4096th move on glued L1."""

    return 1 if (fed.gluing == "isolated" or fed.level == 0) else 4096


def _schedule_task(task: dict[str, Any]) -> dict[str, Any]:
    fed = build_federation(task["level"], task["gluing"])
    loads = initial_loads(task["level"], fed.ports)
    seed = task["seed"]
    law = task["law"]
    started = time.perf_counter()
    if law == "float":
        result = _strip(run_mean_law(fed, loads, seed, engine=task["engine"], max_sweeps=task["max_sweeps"]))
    elif law == "integer":
        result = _strip(run_integer_law(fed, loads, seed, engine=task["engine"], max_sweeps=task["max_sweeps"]))
    elif law == "exact":
        result = run_mean_law_exact(
            fed, loads, seed, engine=task["engine"], max_sweeps=task["max_sweeps"], descent_sample_every=exact_descent_sample_step(fed)
        )
    else:
        raise ValueError(law)
    return {"task": task, "result": result, "seconds": time.perf_counter() - started}


def _run_tasks(tasks: list[dict[str, Any]], workers: int) -> list[dict[str, Any]]:
    if workers <= 1 or len(tasks) <= 1:
        return [_schedule_task(t) for t in tasks]
    import concurrent.futures as futures
    import multiprocessing as mp

    native_kernel()  # compile once before the pool starts
    context = mp.get_context("spawn")
    with futures.ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
        return list(pool.map(_schedule_task, tasks, chunksize=1))


def _aggregate_schedules(entries: list[dict[str, Any]], hash_key: str, expected: str) -> dict[str, Any]:
    hashes = sorted({e[hash_key] for e in entries if e[hash_key] is not None})
    ambiguous = sum(1 for e in entries if e[hash_key] is None)
    return {
        "schedules": len(entries),
        "all_terminated": bool(all(e["terminated"] for e in entries)),
        "unique_terminal_hash_count": len(hashes),
        "ambiguous_terminal_count": int(ambiguous),
        "terminal_hashes_identical_across_schedules": bool(len(hashes) == 1 and ambiguous == 0),
        "terminal_hash_equals_expected": bool(len(hashes) == 1 and ambiguous == 0 and hashes[0] == expected),
        "expected_terminal_hash": expected,
    }


def _stats(values: Sequence[float]) -> dict[str, float]:
    arr = np.asarray(values, dtype=float)
    return {"min": _sig(arr.min(), 6), "median": _sig(float(np.median(arr)), 6), "max": _sig(arr.max(), 6)}


def build_rung_block(
    level: int,
    gluing: str,
    *,
    schedules: int = SCHEDULES,
    engine: str = "auto",
    workers: int = 1,
    laws: Sequence[str] = ("float", "integer", "exact"),
    with_kernel: bool = True,
    max_sweeps: int = MAX_SWEEPS_DEFAULT,
    timings: dict[str, float] | None = None,
) -> dict[str, Any]:
    engine = resolve_engine(engine)
    fed = build_federation(level, gluing)
    loads = initial_loads(level, fed.ports)
    block: dict[str, Any] = {
        "level": level,
        "gluing": gluing,
        "federation": {
            "carriers": fed.carriers,
            "ports": fed.ports,
            "seams": fed.seams,
            "intra_seams": fed.intra_count,
            "inter_seams": fed.inter_count,
            "components": fed.components,
            "component_size_min": int(fed.component_sizes.min()),
            "component_size_max": int(fed.component_sizes.max()),
            "pentagonal_cells": int(np.count_nonzero(fed.pentagonal_cell)),
            "gluing_receipt": fed.gluing_receipt,
            "seam_order": "intra seams carrier-major in carrier.seams() order, then inter seams in dual-edge order",
        },
        "loads": {
            "rule": f"numpy.random.default_rng({LOAD_SEED_BASE} + level).integers(0, {LOAD_MAX + 1}, size=ports)",
            "seed": LOAD_SEED_BASE + level,
            "sha256": sha256_of(loads.tolist()),
            "total": int(loads.sum()),
            "component_totals_sha256": sha256_of(
                np.bincount(fed.component_of_port, weights=loads.astype(float), minlength=fed.components).astype(np.int64).tolist()
            ),
        },
    }
    t0 = time.perf_counter()
    block["synchronous_operator"] = synchronous_operator_receipt(fed)
    if timings is not None:
        timings[f"L{level}/{gluing}/synchronous"] = time.perf_counter() - t0
    seeds = [schedule_seed(level, gluing, s) for s in range(schedules)]
    tasks = []
    for law in laws:
        if law == "exact" and level not in EXACT_LEVELS:
            continue
        for seed in seeds:
            tasks.append({"level": level, "gluing": gluing, "law": law, "seed": seed, "engine": engine, "max_sweeps": max_sweeps})
    t0 = time.perf_counter()
    outcomes = _run_tasks(tasks, workers)
    if timings is not None:
        timings[f"L{level}/{gluing}/schedules"] = time.perf_counter() - t0
        for o in outcomes:
            timings[f"L{level}/{gluing}/{o['task']['law']}/{o['task']['seed']}"] = o["seconds"]
    by_law: dict[str, list[dict[str, Any]]] = {}
    for o in outcomes:
        by_law.setdefault(o["task"]["law"], []).append(o["result"])
    if "float" in by_law:
        entries = sorted(by_law["float"], key=lambda e: e["seed"])
        agg = _aggregate_schedules(entries, "terminal_quotient_hash", fed.expected_terminal_quotient_hash(loads))
        agg.update(
            {
                "threshold": FLOAT_PHI_THRESHOLD,
                "termination": "first sweep end with Phi < threshold; one sweep = |S| attempts",
                "attempts": _stats([e["attempts"] for e in entries]),
                "waits": _stats([e["waits"] for e in entries]),
                "strict_descent_violations_total": int(sum(e["strict_descent_violations"] for e in entries)),
                "descent_ledger_ok_all": bool(all(e["descent_ledger_relative_error_below_1e-9"] for e in entries)),
                "phi_raising_moves": _stats([e["phi_raising_moves"] for e in entries]),
                "phi_raising_fraction_of_non_wait_moves": _stats([e["phi_raising_fraction_of_non_wait_moves"] for e in entries]),
                "phi_raising_counterexample": entries[0]["phi_raising_counterexample"],
                "centered_norm_initial": _round(entries[0]["centered_norm_initial"]),
                "centered_norm_terminal_max": _sig(max(e["centered_norm_terminal"] for e in entries), 6),
                "component_mean_within_1e-9_all": bool(all(e["component_mean_within_1e-9"] for e in entries)),
                "max_abs_deviation_from_component_mean_max": _sig(max(e["max_abs_deviation_from_component_mean"] for e in entries), 3),
                "lattice_residual_max": _sig(max(e["lattice_residual_max"] for e in entries), 3),
                "lattice_snap_unambiguous_all": bool(all(e["lattice_snap_unambiguous"] for e in entries)),
                "entries": entries,
            }
        )
        block["mean_law_float"] = agg
    if "exact" in by_law:
        entries = sorted(by_law["exact"], key=lambda e: e["seed"])
        agg = _aggregate_schedules(entries, "terminal_quotient_hash", fed.expected_terminal_quotient_hash(loads))
        agg.update(
            {
                "threshold": "Phi <= 1/10^18 decided in exact dyadic arithmetic at a sweep end",
                "attempts": _stats([e["attempts"] for e in entries]),
                "descent_identity_exact_all": bool(all(e["descent_identity_exact"] for e in entries)),
                "descent_identity_sample_step": int(entries[0]["descent_identity_sample_step"]),
                "descent_identity_moves_checked_total": int(sum(e["descent_identity_moves_checked"] for e in entries)),
                "global_descent_identity_exact_all": (
                    bool(all(e["global_descent_identity_exact"] for e in entries))
                    if all(e["global_descent_identity_exact"] is not None for e in entries)
                    else None
                ),
                "float_over_exact_phi_ratio_range": [
                    _sig(min(e["float_over_exact_phi_ratio_range"][0] for e in entries), 6),
                    _sig(max(e["float_over_exact_phi_ratio_range"][1] for e in entries), 6),
                ],
                "conservation_exact_all": bool(all(e["total_conservation_exact"] and e["component_conservation_exact"] for e in entries)),
                "lattice_snap_unambiguous_all": bool(all(e["lattice_snap_unambiguous"] for e in entries)),
                "max_abs_deviation_from_component_mean_max": _sig(max(e["max_abs_deviation_from_component_mean"] for e in entries), 3),
                "float_vs_exact_max_abs_deviation_max": _sig(max(e["float_vs_exact_max_abs_deviation"] for e in entries), 3),
                "entries": entries,
            }
        )
        block["mean_law_exact"] = agg
    if "integer" in by_law:
        entries = sorted(by_law["integer"], key=lambda e: e["seed"])
        agg = _aggregate_schedules(entries, "quotient_hash", fed.expected_integer_quotient_hash(loads))
        agg.update(
            {
                "termination": "first attempt after which V equals the balanced-class minimum (exact attempt index)",
                "attempts_to_balanced_class": _stats([e["attempts_to_balanced_class"] for e in entries]),
                "descents": _stats([e["descents"] for e in entries]),
                "swaps": _stats([e["swaps"] for e in entries]),
                "waits": _stats([e["waits"] for e in entries]),
                "odd_tie_seams_at_termination": _stats([e["odd_tie_seams_at_termination"] for e in entries]),
                "max_seam_difference_at_termination": int(max(e["max_seam_difference_at_termination"] for e in entries)),
                "unit_transfers": _stats([e["unit_transfers"] for e in entries]),
                "unit_transfer_decrement_identity_violations_total": int(sum(e["unit_transfer_decrement_identity_violations"] for e in entries)),
                "descent_ledger_exact_all": bool(all(e["descent_ledger_exact"] for e in entries)),
                "strict_descent_violations_total": int(sum(e["strict_descent_violations"] for e in entries)),
                "entries": entries,
            }
        )
        block["integer_law"] = agg
    if level in EXACT_LEVELS:
        t0 = time.perf_counter()
        block["local_confluence"] = local_confluence_receipt(fed, loads)
        if timings is not None:
            timings[f"L{level}/{gluing}/local_confluence"] = time.perf_counter() - t0
    if with_kernel:
        t0 = time.perf_counter()
        block["response_kernel"] = kernel_receipt(fed)
        if timings is not None:
            timings[f"L{level}/{gluing}/kernel"] = time.perf_counter() - t0
    return block


def production_law_note() -> dict[str, Any]:
    return {
        "report_path": PRODUCTION_LAW_REPORT,
        "law": "bw_array overwrite kernel: port_left <- g_ij * port_right or port_right <- inverse(g_ij) * port_left, branch by an independent Bernoulli(0.5) per selected edge, plus a sector-link mutation with probability 0.08",
        "witness": (
            "exact_endpoint_branch_structural_confluence_v1: at a shared node of degree 12 an unchanged incident edge "
            "fixes the node frame, so repairing the left endpoint and repairing the right endpoint leave two distinct "
            "local-frame quotient orbits; unique_terminal_quotient_hash_count = 2; FINITE_CONSENSUS_THEOREM_RECEIPT = false"
        ),
        "reading": (
            "that witness is a property of the overwrite law, which violates clause (ii) of the transactional local "
            "diamond (one canonical aggregate payload per conflict component); the canonical seam-mean law of this "
            "receipt is a different law, symmetric between the two sides of a seam, conservative and idempotent, and its "
            "terminal quotient is the component mean on every schedule"
        ),
    }


def build_receipt(
    *,
    levels: Sequence[int] = RECEIPT_LEVELS,
    gluings: Sequence[str] = GLUINGS,
    schedules: int = SCHEDULES,
    engine: str = "auto",
    workers: int = 1,
    timings: dict[str, float] | None = None,
) -> dict[str, Any]:
    engine = resolve_engine(engine)
    carrier_checks = carrier.check_carrier()
    rungs: dict[str, Any] = {}
    for level in levels:
        for gluing in gluings:
            rungs[f"L{level}/{gluing}"] = build_rung_block(
                level, gluing, schedules=schedules, engine=engine, workers=workers, timings=timings
            )
    summary_rows = []
    for key, block in rungs.items():
        row: dict[str, Any] = {
            "rung": key,
            "carriers": block["federation"]["carriers"],
            "seams": block["federation"]["seams"],
            "components": block["federation"]["components"],
            "laplacian_lambda_2": block["synchronous_operator"]["laplacian_lambda_2"],
        }
        if "mean_law_float" in block:
            f = block["mean_law_float"]
            row.update(
                {
                    "float_all_terminated": f["all_terminated"],
                    "float_attempts_median": f["attempts"]["median"],
                    "float_unique_terminal_hash_count": f["unique_terminal_hash_count"],
                    "float_terminal_hash_equals_expected": f["terminal_hash_equals_expected"],
                }
            )
        if "mean_law_exact" in block:
            e = block["mean_law_exact"]
            row.update(
                {
                    "exact_all_terminated": e["all_terminated"],
                    "exact_attempts_median": e["attempts"]["median"],
                    "exact_unique_terminal_hash_count": e["unique_terminal_hash_count"],
                    "exact_descent_and_conservation": bool(e["descent_identity_exact_all"] and e["conservation_exact_all"]),
                }
            )
        if "integer_law" in block:
            i = block["integer_law"]
            row.update(
                {
                    "integer_attempts_median": i["attempts_to_balanced_class"]["median"],
                    "integer_unique_quotient_hash_count": i["unique_terminal_hash_count"],
                    "integer_odd_tie_seams_median": i["odd_tie_seams_at_termination"]["median"],
                }
            )
        if "response_kernel" in block:
            k = block["response_kernel"]
            row["slow_band_share_n300"] = k["summary"]["300"]["slow_band_share"]
            row["top_four_eigenvalues_median_n300"] = k["summary"]["300"]["top_four_eigenvalues_median"]
            row["isolated_kernel_matches_carrier_within_1e-12"] = k["isolated_kernel_matches_carrier_within_1e-12"]
        summary_rows.append(row)
    all_float_hash = all(r.get("float_terminal_hash_equals_expected", True) for r in summary_rows)
    all_float_term = all(r.get("float_all_terminated", True) for r in summary_rows)
    all_exact = all(r.get("exact_descent_and_conservation", True) and r.get("exact_all_terminated", True) for r in summary_rows)
    all_integer = all(r.get("integer_unique_quotient_hash_count", 1) == 1 for r in summary_rows)
    all_isolated_kernel = all(
        r.get("isolated_kernel_matches_carrier_within_1e-12", True) in (True, None) for r in summary_rows
    )
    receipt = {
        "schema": SCHEMA,
        "lane": "L1 exact federation: confluence and the slow band",
        "claim_boundary": (
            "Finite receipt for the canonical seam-mean law on N = 20 * 4^L twelve-port carriers with a declared "
            "inter-carrier gluing. Supplied: the carrier incidence, the production port-pair gluing convention, seeded "
            "integer loads, the uniform seam schedule and the per-carrier clock normalization of the synchronous "
            "operator. Established: exact single-move descent of V, exact conservation, exact invariance of the "
            "terminal quotient under every seam move, schedule-independent terminal hashes, exact integer-law "
            "termination modulo the odd-tie orbit, and the per-carrier normalized response kernel from the run. "
            "The gluing is a declared convention (the open (M1) item), so the slow-band numbers under gluing are "
            "measurements of that convention. Physical position, physical length, cofinal gluing, the refinement "
            "limit and any field attachment are not claimed."
        ),
        "scope": {
            "supplied": {
                "carrier_incidence": True,
                "gluing_convention": True,
                "initial_loads": True,
                "schedule_law": True,
                "clock_normalization": True,
            },
            "claimed": {
                "physical_position": False,
                "physical_length": False,
                "cofinal_gluing": False,
                "refinement_limit": False,
                "continuum_limit": False,
                "field_attachment": False,
                "source_derived_gluing": False,
            },
        },
        "conventions": {
            "carrier": {
                "ports": PORTS,
                "seams": carrier.SEAM_COUNT,
                "faces": carrier.FACE_COUNT,
                "rotations": carrier.ROTATION_COUNT,
                "antipode": list(carrier.antipode()),
                "self_checks": carrier_checks,
            },
            "tower": "geodesic icosahedral tower, cell basis, N = 20 * 4^L carriers (oph_fpe.core.icosahedral)",
            "mean_law": "E_e = I - b_e b_e^T / 2 on seam e (both endpoints replaced by their mean)",
            "descent_functional": (
                "V(x) = sum_p x_p^2 (flagship descent functional, tex lines 889-895); the seam mean is an orthogonal "
                "projection, so V drops by exactly (x_i - x_j)^2 / 2 per non-wait move; the centered norm V - V_min = "
                "sum_p (x_p - mean_c)^2 drops by the same amount since the total is conserved"
            ),
            "integer_descent": (
                "a nearest-agreement move with mismatch |d| >= 2 is the composition of floor(|d|/2) conservative unit "
                "transfers with mismatches |d|, |d|-2, ...; each lowers V by exactly 2(d_k - 1), so the move lowers V by "
                "(d^2 - (d mod 2)) / 2; the identity is checked on every descent move"
            ),
            "termination_potential": "Phi(x) = sum_seams (x_i - x_j)^2; Phi = 0 is consensus",
            "phi_diagnostic": (
                "Phi is not monotone along single moves: Phi(E_e x) - Phi(x) = -d ((deg_i x_i - N_i) - (deg_j x_j - N_j)) "
                "+ d^2 (deg_i + deg_j + 2) / 4 with N the neighbour sums; the count of accepted moves raising Phi and the "
                "first such move of each schedule are reported as diagnostics"
            ),
            "phi_single_move_witness": phi_nonmonotone_witness(),
            "schedule": (
                "uniform over all seams: per sweep numpy.random.default_rng(seed).integers(0, |S|, size=|S|); "
                "the integer law draws integers(0, 2, size=|S|) tie coins after each sweep's seam draw"
            ),
            "schedule_seeds": f"{SCHEDULE_SEED_BASE} + 1000*level + 100*gluing_index + schedule_index, gluing_index in {dict(zip(GLUINGS, range(len(GLUINGS))))}",
            "integer_law": (
                "carrier.integer_nearest_agreement: (floor(s/2), ceil(s/2)) with the ceiling placed by the coin (1/2 each); "
                "a unit-difference seam is a wait in one placement and a swap in the other"
            ),
            "integer_terminal_class": "balanced class: every reading in {q, q+1} per component; one orbit of the unit swaps",
            "terminal_canonicalizer_mean_law": (
                "component-lattice snap: q_p = round(m_c x_p) with m_c the port count of the component of p; the hash "
                "covers component sizes, component labels and q; the maximal residual |m_c x_p - q_p| is reported"
            ),
            "terminal_canonicalizer_integer_law": "multiset of readings per component, components ordered by lowest port",
            "float_rounding_in_receipt": f"{FLOAT_ROUND} decimals for kernel data, {PHI_TRACE_ROUND} significant digits for Phi",
            "exact_arithmetic": "dyadic integers x_p = n_p / 2^K with the parity of every endpoint total asserted",
            "engines": "python reference, layered numpy and native C apply identical IEEE operations; states are bit-identical",
        },
        "rungs": rungs,
        "summary": summary_rows,
        "verdicts": {
            "float_terminal_hash_identical_and_expected_all_rungs": bool(all_float_hash),
            "float_all_schedules_terminated_all_rungs": bool(all_float_term),
            "exact_rungs_descent_conservation_termination": bool(all_exact),
            "integer_quotient_unique_all_rungs": bool(all_integer),
            "isolated_kernel_equals_carrier_all_rungs": bool(all_isolated_kernel),
            "carrier_self_checks": bool(all(carrier_checks.values())),
        },
        "production_law_note": production_law_note(),
        "pins": file_pins(),
        "open_items": [
            "source-derived gluing (M1) is work in progress; the port_pair convention is declared",
            "the refinement limit of the glued slow band across levels is work in progress",
        ],
    }
    return receipt


def file_pins() -> dict[str, str | None]:
    pins: dict[str, str | None] = {}
    for name, path in (
        ("oph_exact/federation.py", PRODUCER_PATH),
        ("oph_exact/verify_federation_independent.py", VERIFIER_PATH),
        ("tests/test_exact_federation.py", TEST_PATH),
        ("oph_exact/carrier.py", CARRIER_PATH),
    ):
        pins[name] = sha256_file(path) if Path(path).exists() else None
    return pins


def write_receipt(receipt: dict[str, Any], path: Path = RECEIPT_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(canonical_json(receipt), encoding="ascii")


def load_receipt(path: Path = RECEIPT_PATH) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="ascii"))


def check_receipt(path: Path = RECEIPT_PATH, **kwargs: Any) -> None:
    stored = load_receipt(path)
    rebuilt = build_receipt(**kwargs)
    if canonical_json(stored) != canonical_json(rebuilt):
        raise RuntimeError(f"stale receipt: {path} differs from an in-memory rebuild")


# --------------------------------------------------------------------------
# Scale demonstration (runs/ summaries)
# --------------------------------------------------------------------------


def scale_run(
    level: int,
    gluing: str,
    schedules: int,
    out: Path,
    *,
    engine: str = "auto",
    workers: int = 1,
    max_sweeps: int = 256,
    kernel_cells: int = 8,
) -> dict[str, Any]:
    """Run the federation at scale with a declared sweep budget and write a summary."""

    engine = resolve_engine(engine)
    timings: dict[str, float] = {}
    t0 = time.perf_counter()
    fed = build_federation(level, gluing)
    timings["build_seconds"] = time.perf_counter() - t0
    loads = initial_loads(level, fed.ports)
    t0 = time.perf_counter()
    sync = synchronous_operator_receipt(fed)
    timings["spectrum_seconds"] = time.perf_counter() - t0
    tasks = []
    for s in range(schedules):
        seed = schedule_seed(level, gluing, s)
        tasks.append({"level": level, "gluing": gluing, "law": "integer", "seed": seed, "engine": engine, "max_sweeps": MAX_SWEEPS_DEFAULT})
        tasks.append({"level": level, "gluing": gluing, "law": "float", "seed": seed, "engine": engine, "max_sweeps": max_sweeps})
    t0 = time.perf_counter()
    outcomes = _run_tasks(tasks, workers)
    timings["schedules_seconds"] = time.perf_counter() - t0
    floats = sorted([o["result"] for o in outcomes if o["task"]["law"] == "float"], key=lambda e: e["seed"])
    ints = sorted([o["result"] for o in outcomes if o["task"]["law"] == "integer"], key=lambda e: e["seed"])
    for entry, o in zip(floats, [o for o in outcomes if o["task"]["law"] == "float"]):
        entry["seconds"] = round(o["seconds"], 3)
        entry["attempts_per_second"] = _sig(entry["attempts"] / max(o["seconds"], 1e-9), 4)
    for entry, o in zip(ints, [o for o in outcomes if o["task"]["law"] == "integer"]):
        entry["seconds"] = round(o["seconds"], 3)
    t0 = time.perf_counter()
    cells = sample_cells(fed, kernel_cells)
    kernel = kernel_receipt(fed, cells)
    timings["kernel_seconds"] = time.perf_counter() - t0
    lam2 = sync["laplacian_lambda_2"]
    extrapolated = None
    if floats and lam2 and lam2 > 0:
        phi_end = floats[0]["phi_terminal"]
        if phi_end > FLOAT_PHI_THRESHOLD:
            e_folds = math.log(phi_end / FLOAT_PHI_THRESHOLD)
            extrapolated = _sig(floats[0]["attempts"] + e_folds * fed.seams / lam2, 4)
    summary = {
        "schema": "oph.exact.federation-scale-run.v1",
        "level": level,
        "gluing": gluing,
        "engine": engine,
        "native_kernel_error": _NATIVE["error"],
        "workers": workers,
        "carriers": fed.carriers,
        "ports": fed.ports,
        "seams": fed.seams,
        "inter_seams": fed.inter_count,
        "components": fed.components,
        "gluing_receipt": {k: v for k, v in fed.gluing_receipt.items() if k != "port_pairs"},
        "synchronous_operator": sync,
        "float_sweep_budget": max_sweeps,
        "mean_law_float": floats,
        "float_unique_terminal_hash_count": len({h for e in floats if (h := e["terminal_quotient_hash"]) is not None}),
        "float_ambiguous_terminal_count": sum(1 for e in floats if e["terminal_quotient_hash"] is None),
        "float_extrapolated_attempts_to_threshold": extrapolated,
        "integer_law": ints,
        "integer_unique_quotient_hash_count": len({e["quotient_hash"] for e in ints}),
        "integer_expected_quotient_hash": fed.expected_integer_quotient_hash(loads),
        "response_kernel": kernel,
        "timings_seconds": {k: round(v, 3) for k, v in timings.items()},
        "wall_clock_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(canonical_json(summary), encoding="ascii")
    return summary


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--write", action="store_true", help="build and write the receipt")
    parser.add_argument("--check", action="store_true", help="rebuild and compare with the stored receipt (stale raises)")
    parser.add_argument("--level", type=int, default=None, help="scale run: tower level")
    parser.add_argument("--gluing", default="port_pair", choices=GLUINGS)
    parser.add_argument("--schedules", type=int, default=4)
    parser.add_argument("--max-sweeps", type=int, default=256, help="scale run: float sweep budget")
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--engine", default="auto", choices=("auto", "native", "numpy", "python"))
    parser.add_argument("--workers", type=int, default=max(1, min(8, (os.cpu_count() or 2) - 2)))
    parser.add_argument("--receipt", type=Path, default=RECEIPT_PATH)
    args = parser.parse_args(argv)
    if args.write or args.check:
        timings: dict[str, float] = {}
        started = time.perf_counter()
        if args.check:
            check_receipt(args.receipt, engine=args.engine, workers=args.workers, timings=timings)
            print(f"receipt fresh: {args.receipt} ({time.perf_counter() - started:.1f} s)")
        else:
            receipt = build_receipt(engine=args.engine, workers=args.workers, timings=timings)
            write_receipt(receipt, args.receipt)
            total = time.perf_counter() - started
            summary_dir = REPO_ROOT / "runs" / "exact_federation_receipt"
            summary_dir.mkdir(parents=True, exist_ok=True)
            (summary_dir / "timing.json").write_text(
                canonical_json(
                    {
                        "engine": resolve_engine(args.engine),
                        "workers": args.workers,
                        "total_seconds": round(total, 2),
                        "timings_seconds": {k: round(v, 3) for k, v in timings.items()},
                        "wall_clock_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    }
                ),
                encoding="ascii",
            )
            print(f"wrote {args.receipt} in {total:.1f} s; verdicts: {receipt['verdicts']}")
        return 0
    if args.level is not None:
        out = args.out or (REPO_ROOT / "runs" / f"exact_federation_L{args.level}")
        started = time.perf_counter()
        summary = scale_run(
            args.level,
            args.gluing,
            args.schedules,
            out,
            engine=args.engine,
            workers=args.workers,
            max_sweeps=args.max_sweeps,
        )
        print(
            f"level {args.level} {args.gluing}: {summary['carriers']} carriers, {summary['seams']} seams, "
            f"engine {summary['engine']}, {time.perf_counter() - started:.1f} s; summary at {out / 'summary.json'}"
        )
        return 0
    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
