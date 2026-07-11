#!/usr/bin/env python3
"""Machine receipts for the conditional Lorentzian event manifold (GitHub #525).

Implements a synthetic end-to-end run of the compact paper's subsection
`subsec:event-manifold`:

* ground truth: events (tau, X) with X on the unit-radius hyperboloid H^3
  embedded in R^{1,3}; located records carry calibrated cap responses
  h_j = eta(X, n_j) (Proposition `prop:h3-general-cap-frame-condition`), one
  modular clock read, and stage-indexed noise/mesh moduli shrinking under
  refinement (receipts E1-E3);
* localization + germ formation: per-stage frame reconstruction (the
  closed-form left inverse of the rank-4 response frame), certified boxes,
  and coincidence classes across stages (Definition `def:record-coincidence`);
* verdicts: SEPARATED (certified disjoint boxes at some stage), COINCIDENT
  (boxes intersect cofinally with shrinking radii), AMBIGUOUS (neither
  certificate: the non-Hausdorff detection of
  Proposition `prop:event-countermodels`(iii));
* chart receipts (Theorem `thm:event-4-manifold`): bi-Lipschitz distortion of
  the reconstructed chart against ground truth, local dimension 4 by PCA
  rank;
* causal/signature receipts (Theorem `thm:event-causal-structure`): from the
  intrinsic precedence labels alone (record ancestry = cone order in the
  synthetic model) plus reconstructed chart coordinates, fit a quadratic form
  separating causal from spacelike displacements and verify Lorentz signature
  (-+++), one time direction, three spatial directions;
* invariance receipts (Lemma `lem:event-quotient-invariance`): global Lorentz
  boost of all record data, gauge relabeling, and schedule permutation leave
  the event classes and all pairwise eta-intervals unchanged;
* countermodels (Proposition `prop:event-countermodels`): population deficit
  drops the PCA dimension; label inflation doubles the event count with
  identical kinematics; a seeded reconciliation branching yields AMBIGUOUS.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


def eta_dot(u: np.ndarray, v: np.ndarray) -> float:
    return float(u @ ETA @ v)


def embed_h3(y: np.ndarray) -> np.ndarray:
    """Map y in R^3 to the unit hyperboloid H^3 in R^{1,3}."""
    return np.concatenate(([np.sqrt(1.0 + y @ y)], y))


# ---------------------------------------------------------------------------
# response frame (rank-4: three-plus clock) and closed-form reconstruction
# ---------------------------------------------------------------------------

def icosahedral_cap_normals() -> np.ndarray:
    """Twelve unit spacelike cap normals from the centered tight frame of
    Definition def:centered-tight-cap-frame."""
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    dirs = []
    for s1 in (-1.0, 1.0):
        for s2 in (-1.0, 1.0):
            dirs += [
                (0.0, s1, s2 * phi),
                (s1, s2 * phi, 0.0),
                (s2 * phi, 0.0, s1),
            ]
    u = np.array(dirs)
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    normals = np.stack(
        [np.concatenate(([1.0 / np.sqrt(2.0)], np.sqrt(1.5) * d)) for d in u]
    )
    return normals


def response_matrix(normals: np.ndarray) -> np.ndarray:
    """Rows (-n^0, n^1, n^2, n^3): h(X) = B X_E with X_E Euclideanized."""
    b = normals.copy()
    b[:, 0] *= -1.0
    return b


def frame_reconstruct(responses: np.ndarray, b_matrix: np.ndarray) -> np.ndarray:
    """Closed-form left inverse of Proposition prop:h3-general-cap-frame-condition."""
    x_e, *_ = np.linalg.lstsq(b_matrix, responses, rcond=None)
    return x_e


# ---------------------------------------------------------------------------
# synthetic record layer
# ---------------------------------------------------------------------------

@dataclass
class LocatedRecord:
    token: int
    stage: int
    tau: float
    x_e: np.ndarray       # reconstructed Euclideanized 4-vector (t, x, y, z)
    radius: float         # certified box radius (space + clock combined)
    certified: bool       # residual-gap certificate available at this stage


def stage_noise(stage: int, base: float = 0.2) -> float:
    return base * (0.5 ** stage)


def generate_records(events: np.ndarray, taus: np.ndarray, stages: int,
                     seed: int = 5, certify: bool = True) -> list[LocatedRecord]:
    """Per stage, one located record per event with shrinking noise (E1)."""
    rng = np.random.default_rng(seed)
    normals = icosahedral_cap_normals()
    b = response_matrix(normals)
    records = []
    for stage in range(stages):
        sigma = stage_noise(stage)
        for token, (x, tau) in enumerate(zip(events, taus)):
            h = np.array([eta_dot(x, n) for n in normals])
            h_noisy = h + rng.normal(scale=sigma / np.sqrt(len(normals)), size=len(h))
            x_rec = frame_reconstruct(h_noisy, b)
            tau_rec = tau + rng.normal(scale=sigma)
            records.append(
                LocatedRecord(
                    token=token,
                    stage=stage,
                    tau=tau_rec,
                    x_e=x_rec,
                    radius=3.0 * sigma + 1e-6,
                    certified=certify,
                )
            )
    return records


# ---------------------------------------------------------------------------
# germs, coincidence, verdicts
# ---------------------------------------------------------------------------

def chart_point(rec: LocatedRecord) -> np.ndarray:
    return np.concatenate(([rec.tau], rec.x_e[1:]))


def germs_from_records(records: list[LocatedRecord]) -> dict[int, list[LocatedRecord]]:
    germs: dict[int, list[LocatedRecord]] = {}
    for rec in records:
        germs.setdefault(rec.token, []).append(rec)
    for chain in germs.values():
        chain.sort(key=lambda r: r.stage)
    return germs


def pair_verdict(g1: list[LocatedRecord], g2: list[LocatedRecord]) -> str:
    """SEPARATED / COINCIDENT / AMBIGUOUS from certified boxes only."""
    stages = range(min(len(g1), len(g2)))
    separated_certified = False
    intersect_tail = True
    for s in stages:
        r1, r2 = g1[s], g2[s]
        gap = np.linalg.norm(chart_point(r1) - chart_point(r2))
        boxes_disjoint = gap > (r1.radius + r2.radius)
        if boxes_disjoint:
            intersect_tail = False
            if r1.certified and r2.certified:
                separated_certified = True
    if separated_certified:
        return "SEPARATED"
    if intersect_tail:
        return "COINCIDENT"
    return "AMBIGUOUS"


def event_classes(records: list[LocatedRecord]) -> tuple[list[list[int]], list[tuple]]:
    """Union-find over COINCIDENT germs; returns classes and ambiguous pairs."""
    germs = germs_from_records(records)
    tokens = sorted(germs)
    parent = {t: t for t in tokens}

    def find(t: int) -> int:
        while parent[t] != t:
            parent[t] = parent[parent[t]]
            t = parent[t]
        return t

    ambiguous = []
    for i, t1 in enumerate(tokens):
        for t2 in tokens[i + 1:]:
            verdict = pair_verdict(germs[t1], germs[t2])
            if verdict == "COINCIDENT":
                parent[find(t1)] = find(t2)
            elif verdict == "AMBIGUOUS":
                ambiguous.append((t1, t2))
    classes: dict[int, list[int]] = {}
    for t in tokens:
        classes.setdefault(find(t), []).append(t)
    return list(classes.values()), ambiguous


def class_positions(records: list[LocatedRecord], classes: list[list[int]]) -> np.ndarray:
    """Final-stage chart position per event class."""
    germs = germs_from_records(records)
    out = []
    for cls in classes:
        pts = [chart_point(germs[t][-1]) for t in cls]
        out.append(np.mean(pts, axis=0))
    return np.array(out)


# ---------------------------------------------------------------------------
# chart, dimension, and signature receipts
# ---------------------------------------------------------------------------

def bilipschitz_distortion(true_pts: np.ndarray, rec_pts: np.ndarray) -> float:
    """max over pairs of ratio between reconstructed and true chart distance."""
    ratios = []
    n = len(true_pts)
    for i in range(n):
        for j in range(i + 1, n):
            dt = np.linalg.norm(true_pts[i] - true_pts[j])
            dr = np.linalg.norm(rec_pts[i] - rec_pts[j])
            if dt > 1e-9:
                ratios.append(dr / dt)
    ratios = np.array(ratios)
    return float(max(ratios.max(), 1.0 / ratios.min()))


def pca_dimension(points: np.ndarray, rel_threshold: float = 1e-2) -> int:
    centered = points - points.mean(axis=0)
    svals = np.linalg.svd(centered, compute_uv=False)
    if svals[0] < 1e-12:
        return 0
    return int(np.sum(svals / svals[0] > rel_threshold))


def causal_labels(chart_pts: np.ndarray) -> list[tuple[int, int, bool]]:
    """Intrinsic precedence labels of the synthetic transactional ancestry:
    (i, j, causal?) for displacement chart_pts[j] - chart_pts[i]."""
    labels = []
    n = len(chart_pts)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            v = chart_pts[j] - chart_pts[i]
            labels.append((i, j, eta_dot(v, v) < 0.0 and v[0] > 0.0))
    return labels


def fit_cone_quadric(chart_pts: np.ndarray,
                     labels: list[tuple[int, int, bool]]) -> np.ndarray:
    """Least-squares quadratic form q with q(v) ~ -1 on causal displacements
    and q(v) ~ +1 on spacelike ones, from labels + chart coordinates only."""
    rows, targets = [], []
    for i, j, causal in labels:
        v = chart_pts[j] - chart_pts[i]
        v = v / (np.linalg.norm(v) + 1e-12)
        outer = np.outer(v, v)
        rows.append(outer[np.triu_indices(4)] * (2.0 - np.eye(4))[np.triu_indices(4)])
        targets.append(-1.0 if causal else 1.0)
    coef, *_ = np.linalg.lstsq(np.array(rows), np.array(targets), rcond=None)
    q = np.zeros((4, 4))
    q[np.triu_indices(4)] = coef
    q = (q + q.T) / 2.0
    return q


def signature_of(q: np.ndarray, tol: float = 1e-6) -> tuple[int, int]:
    evals = np.linalg.eigvalsh(q)
    return int(np.sum(evals < -tol)), int(np.sum(evals > tol))


# ---------------------------------------------------------------------------
# invariance drivers
# ---------------------------------------------------------------------------

def boost_matrix(rapidity: float, axis: int = 0) -> np.ndarray:
    lam = np.eye(4)
    c, s = np.cosh(rapidity), np.sinh(rapidity)
    j = axis + 1
    lam[0, 0] = c
    lam[0, j] = s
    lam[j, 0] = s
    lam[j, j] = c
    return lam


def pairwise_intervals(chart_pts: np.ndarray) -> np.ndarray:
    n = len(chart_pts)
    out = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            v = chart_pts[j] - chart_pts[i]
            out[i, j] = eta_dot(v, v)
    return out


# ---------------------------------------------------------------------------
# countermodel generators (Proposition prop:event-countermodels)
# ---------------------------------------------------------------------------

def worldline_events(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Population deficit: events along a single timelike worldline."""
    taus = np.linspace(0.0, 3.0, n)
    xs = np.stack([embed_h3(np.zeros(3)) for _ in taus])
    return xs, taus


def generic_events(n: int, seed: int = 9) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    ys = rng.uniform(-1.0, 1.0, size=(n, 3))
    xs = np.stack([embed_h3(y) for y in ys])
    taus = rng.uniform(0.0, 3.0, size=n)
    return xs, taus


def label_inflated_records(events: np.ndarray, taus: np.ndarray, stages: int,
                           seed: int = 5) -> list[LocatedRecord]:
    """Doubled record family with a persistent internal label surviving all
    separation tests: every event appears as two never-coincident germs
    displaced by a label offset that does NOT shrink under refinement."""
    base = generate_records(events, taus, stages, seed=seed)
    # (tau, x_e) block offset: the label shifts a chart-visible coordinate
    offset = np.array([0.0, 0.0, 10.0, 0.0, 0.0])
    doubled = []
    n = len(events)
    for rec in base:
        doubled.append(rec)
        shifted = LocatedRecord(
            token=rec.token + n,
            stage=rec.stage,
            tau=rec.tau + offset[0],
            x_e=rec.x_e + offset[1:],
            radius=rec.radius,
            certified=rec.certified,
        )
        doubled.append(shifted)
    return doubled


def branching_records(stages: int, split_stage: int = 2,
                      seed: int = 13) -> list[LocatedRecord]:
    """Reconciliation branching: two germ families share every certified box
    up to split_stage, then separate with permanently uncertified gaps
    (Delta_loc <= 0): the finite-stage witness of non-Hausdorff formation."""
    x = embed_h3(np.array([0.2, -0.1, 0.4]))
    tau = 1.0
    records = generate_records(np.stack([x]), np.array([tau]), stages, seed=seed)
    twin = []
    for rec in records:
        drift = 0.0 if rec.stage < split_stage else 0.5
        twin.append(
            LocatedRecord(
                token=1,
                stage=rec.stage,
                tau=rec.tau + drift,
                x_e=rec.x_e.copy(),
                radius=rec.radius,
                certified=rec.stage < split_stage,  # gap never certified after split
            )
        )
    # the original germ also loses certification after the split
    original = [
        LocatedRecord(rec.token, rec.stage, rec.tau, rec.x_e, rec.radius,
                      certified=rec.stage < split_stage)
        for rec in records
    ]
    return original + twin
