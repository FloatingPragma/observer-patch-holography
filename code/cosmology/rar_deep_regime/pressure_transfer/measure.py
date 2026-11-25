"""Brouwer Figure 8 type contrasts in units of the held SPARC envelope."""

import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile

import numpy as np
from scipy.optimize import lsq_linear, minimize_scalar

HERE = Path(__file__).resolve().parent
FROZEN = HERE.parent / "receipts/lensing_profile.json"


def extract(path, frozen):
    archive_sha = hashlib.sha256(path.read_bytes()).hexdigest()
    if archive_sha != frozen["archive_sha256"]:
        raise ValueError("Unrecognized release checksum")
    edges = np.geomspace(1e-15, 5e-12, 16)
    keep = (edges[:-1] >= 1e-13) & (edges[1:] <= 1.2e-11)
    result = dict(
        archive_sha256=archive_sha, bin_edges_m_s2=edges[9:].tolist(), splits={}
    )
    np.testing.assert_array_equal(np.flatnonzero(keep), np.arange(9, 15))
    with tarfile.open(path) as archive:
        for split, threshold in (("Color", 2.5), ("Sersic", 2.0)):
            prefix = f"Fig-8_RAR-KiDS-isolated_{split}"
            tables = {}
            hashes = {}
            for suffix in ("bin_1", "bin_2", "bins_covmatrix"):
                name = f"{prefix}{suffix}.txt"
                with archive.extractfile(name) as member:
                    data = member.read()
                hashes[name] = hashlib.sha256(data).hexdigest()
                tables[suffix] = np.loadtxt(io.BytesIO(data))
            raw = tables["bins_covmatrix"]
            row_type, ti = np.unique(raw[:, 0], return_inverse=True)
            col_type, tj = np.unique(raw[:, 1], return_inverse=True)
            b, bi = np.unique(raw[:, 2], return_inverse=True)
            other_b, bj = np.unique(raw[:, 3], return_inverse=True)
            np.testing.assert_array_equal(row_type, [0.0, threshold])
            np.testing.assert_array_equal(row_type, col_type)
            np.testing.assert_array_equal(b, other_b)
            np.testing.assert_allclose(
                b, (edges[:-1] + edges[1:]) / 2, rtol=2e-10, atol=0
            )
            i, j = 15 * ti + bi, 15 * tj + bj
            np.testing.assert_array_equal(np.sort(30 * i + j), np.arange(900))
            cov = np.empty((30, 30))
            cov[i, j] = raw[:, 4] / raw[:, 6]
            np.testing.assert_allclose(cov, cov.T, rtol=5e-10, atol=1e-11)
            stacks = np.stack((tables["bin_1"], tables["bin_2"]))
            for stack in stacks:
                np.testing.assert_allclose(stack[:, 0], b, rtol=5e-5, atol=0)
            errors = (stacks[:, :, 3] / stacks[:, :, 4]).ravel()
            np.testing.assert_allclose(np.sqrt(cov.diagonal()), errors, rtol=2e-4)
            both = np.tile(keep, 2)
            result["splits"][split] = dict(
                member_sha256=hashes,
                bin_centers_m_s2=b[keep].tolist(),
                tangential_esd=(stacks[:, keep, 1] / stacks[:, keep, 4]).tolist(),
                cross_esd=(stacks[:, keep, 2] / stacks[:, keep, 4]).tolist(),
                covariance=cov[np.ix_(both, both)].tolist(),
            )
    return result


def weights(template, cov):
    precision_template = np.linalg.solve(cov, template)
    return precision_template / (template @ precision_template)


def chisq(residual, cov):
    return float(residual @ np.linalg.solve(cov, residual))


def joint_edge_distance(y, cov, edges, a0, four_g):
    """Project both stacks onto all possible within-bin means of a common law."""
    lo, hi = np.tile(edges[:-1], 2), np.tile(edges[1:], 2)
    core_lo, core_hi = (4 / np.pi) * lo / four_g, (4 / np.pi) * hi / four_g
    env_lo, env_hi = np.sqrt(a0 * lo) / four_g, np.sqrt(a0 * hi) / four_g
    white = np.linalg.solve(np.linalg.cholesky(cov), np.eye(len(y)))
    white_y = white @ y

    def fit(q):
        result = lsq_linear(
            white,
            white_y,
            bounds=(core_lo + q * env_lo, core_hi + q * env_hi),
            method="bvls",
        )
        if not result.success:
            raise RuntimeError(result.message)
        return result

    def certify(q, common_q_is_free):
        fitted = fit(q)
        upper = chisq(y - fitted.x, cov)
        v = np.linalg.solve(cov, y - fitted.x)
        if upper <= 1e-12:
            return dict(
                envelope_relative_to_sparc=float(q),
                distance_upper=upper,
                distance_lower=0.0,
            )
        if common_q_is_free:
            # Make the support of v over every q >= 0 finite. Correct tiny
            # numerical error in the stationary slope of the supporting plane.
            slope = v @ np.where(v > 0, env_hi, env_lo)
            j = int(np.argmin(v))
            if v[j] >= 0:
                raise ValueError(
                    "No negative component to certify the supporting plane"
                )
            v[j] -= max(0.0, slope + 1e-12) / env_lo[j]
            if v @ np.where(v > 0, env_hi, env_lo) > 0:
                raise ValueError("Unbounded dual support")
            support = v @ np.where(v > 0, core_hi, core_lo)
        else:
            support = v @ np.where(v > 0, core_hi + q * env_hi, core_lo + q * env_lo)
        # Every allowed mean m obeys v.m <= support. Complete the square in
        # (y-m)' C^-1 (y-m) to certify a global lower bound on the distance.
        lower = 2 * (v @ y - support) - v @ cov @ v
        if not -1e-9 <= upper - lower <= 1e-6:
            raise ValueError("Joint global distance certificate failed")
        return dict(
            envelope_relative_to_sparc=float(q),
            distance_upper=upper,
            distance_lower=float(lower),
        )

    limit = 1.0
    while fit(limit).cost > fit(2 * limit).cost:
        limit *= 2
    result = minimize_scalar(
        lambda q: fit(q).cost,
        bounds=(0, 2 * limit),
        method="bounded",
        options={"xatol": 1e-12},
    )
    if not result.success:
        raise RuntimeError(result.message)
    # The certificate also covers q outside the numerical search bracket.
    best_q = 0.0 if fit(0.0).cost <= result.fun else result.x
    return dict(
        frozen=certify(1.0, False), any_nonnegative_envelope=certify(best_q, True)
    )


def read_split(data, edges, a0, four_g):
    b = np.array(data["bin_centers_m_s2"])
    core = (4 / np.pi) * b / four_g
    env = np.sqrt(a0 * b) / four_g
    y = np.array(data["tangential_esd"]).ravel()
    cross = np.array(data["cross_esd"]).ravel()
    joint = np.array(data["covariance"])
    np.linalg.cholesky(joint)
    ident = np.eye(len(b))
    diff = np.hstack((-ident, ident))
    cov = diff @ joint @ diff.T
    w = weights(env, cov)
    d = diff @ y
    q = float(w @ d)
    se = float(np.sqrt(w @ cov @ w))
    cross_q = float(w @ (diff @ cross))
    null = chisq(d, cov)
    residual = chisq(d - q * env, cov)
    np.testing.assert_allclose(null - residual, (q / se) ** 2)

    # All physical means over pairs of a monotone law lie between its edge values.
    # For any common q >= 0, sqrt(upper) = k sqrt(lower) cancels q exactly.
    low, high = edges[:-1], edges[1:]
    k = float(np.sqrt(high[0] / low[0]))
    np.testing.assert_allclose(np.sqrt(high / low), k, rtol=1e-12)
    edge_diff = np.hstack((-k * ident, ident))
    edge_cov = edge_diff @ joint @ edge_diff.T
    edge_weights = weights(env, edge_cov)
    if np.min(edge_weights) < 0:
        raise ValueError("A signed weight cannot certify the stacked bin-edge bound")
    upper = (4 / np.pi) * (high - k * low) / four_g
    edge_excess = float(edge_weights @ (edge_diff @ y - upper))
    edge_se = float(np.sqrt(edge_weights @ edge_cov @ edge_weights))

    at_centers = {}
    for label, sl in (("low", slice(0, len(b))), ("high", slice(len(b), None))):
        type_cov = joint[sl, sl]
        type_w = weights(env, type_cov)
        type_q = float(type_w @ (y[sl] - core))
        at_centers[label] = dict(
            envelope_relative_to_sparc=type_q,
            statistical_se=float(np.sqrt(type_w @ type_cov @ type_w)),
            residual_chisq=chisq(y[sl] - core - type_q * env, type_cov),
        )
    return dict(
        difference_at_centers=dict(
            envelope_high_minus_low=q,
            statistical_se=se,
            cross_high_minus_low=cross_q,
            chisq_identical_profiles=null,
            identical_profiles_residual_dimension=len(b),
            chisq_after_one_envelope_difference=residual,
            envelope_difference_residual_dimension=len(b) - 1,
        ),
        independent_within_bin_weighting=dict(
            sqrt_edge_ratio=k,
            high_minus_k_low_minus_upper_bound=edge_excess,
            statistical_se=edge_se,
            excess_in_statistical_se=edge_excess / edge_se,
            weights=edge_weights.tolist(),
        ),
        minimum_joint_within_bin_chisq=joint_edge_distance(y, joint, edges, a0, four_g),
        fixed_core_single_stack_diagnostics=at_centers,
    )


def analyse(profiles, frozen):
    edges = np.array(profiles["bin_edges_m_s2"])
    return {
        name: read_split(
            data,
            edges,
            frozen["fitted_sparc_a0_m_s2"],
            frozen["four_G_m_s2_per_Msun_pc2"],
        )
        for name, data in profiles["splits"].items()
    }


def check_values(a, b):
    if isinstance(a, dict):
        if set(a) != set(b):
            raise AssertionError("Result keys differ")
        for key in a:
            check_values(a[key], b[key])
    else:
        np.testing.assert_allclose(a, b, rtol=1e-7, atol=1e-11)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--archive", type=Path, help="Verify extraction against the pinned public tar"
    )
    p.add_argument(
        "--write", action="store_true", help="Write receipts from the pinned public tar"
    )
    args = p.parse_args()
    if args.write and not args.archive:
        p.error("--write requires --archive")
    frozen = json.loads(FROZEN.read_text())
    path = HERE / "profiles.json"
    profiles = (
        extract(args.archive, frozen) if args.archive else json.loads(path.read_text())
    )
    results = analyse(profiles, frozen)
    for file, obj in ((path, profiles), (HERE / "results.json", results)):
        if args.write:
            file.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")
        elif file == path:
            if obj != json.loads(file.read_text()):
                raise AssertionError("Published profiles differ from stored extraction")
        else:
            check_values(obj, json.loads(file.read_text()))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
