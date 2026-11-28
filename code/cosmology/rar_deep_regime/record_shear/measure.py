"""Stationary mark visits, the sign of the projected curvature, and a shear test."""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import tarfile

import numpy as np
from scipy.fft import dstn, idstn
from scipy.special import ndtr, ndtri

from matched_einstein_limit.collar_identifiability.refinement import (
    boundary_sources as boundary,
)
from .. import lensing_profile

HERE = Path(__file__).parent
RECEIPT = HERE / "results.json"
SHIFTS = ((11, 12, 15), (15, 16, 18), (17, 13, 11), (15, 12, 13))
WIDTH = 31


def laplacian(x):
    padded = np.pad(x, 1)
    core = tuple(slice(1, -1) for _ in x.shape)
    total = np.zeros_like(x)
    for ax in range(x.ndim):
        for region in (slice(None, -2), slice(2, None)):
            ind = list(core)
            ind[ax] = region
            total += padded[tuple(ind)] - x
    return total


def finite_walk(survival):
    # Descent from a triangular face fixes the four native source cells.
    board = boundary.federation(1)
    sites = list(range(4))
    incidence = np.stack(
        [boundary.boundary_coefficients(board, sites, v) for v in np.eye(4)], axis=1
    )
    incidence = incidence[np.any(incidence, axis=1)]
    # The chart locations and absorbing box are declared comparison choices.
    xyz = np.arange(1, WIDTH + 1)
    c = np.cos(np.pi * xyz / (WIDTH + 1))
    eigen = (
        1 - survival * (1 + c[:, None, None] + c[None, :, None] + c[None, None, :]) / 4
    )
    injections = np.zeros((4, WIDTH, WIDTH, WIDTH))
    for u, where in zip(injections, SHIFTS, strict=True):
        u[where] = 1
    visits = np.stack(
        [
            idstn(dstn(v, type=1, norm="ortho") / eigen, type=1, norm="ortho")
            for v in injections
        ]
    )
    source = np.einsum("ei,ixyz->exyz", incidence, injections)
    fields = np.einsum("ei,ixyz->exyz", incidence, visits)
    squares = (fields**2).sum(axis=0)
    gradient = np.zeros_like(squares)
    for field in fields:
        padded = np.pad(field, 1)
        for ax in range(3):
            for region in (slice(None, -2), slice(2, None)):
                cut = [slice(1, -1)] * 3
                cut[ax] = region
                gradient += (padded[tuple(cut)] - field) ** 2
    rhs = (
        gradient
        + 16 * (1 / survival - 1) * squares
        - 16 / survival * (source * fields).sum(axis=0)
    )
    error = float(abs(laplacian(squares) - rhs).max())
    if error > 1e-11:
        raise ArithmeticError("Squared-visit identity failed")
    projected = squares.sum(axis=2)
    curvature = laplacian(projected)
    # Summing the z-Laplacian over the interior leaves the two boundary faces.
    check = rhs.sum(axis=2) + squares[:, :, 0] + squares[:, :, -1]
    np.testing.assert_allclose(curvature, check, rtol=0, atol=1e-11)
    injecting = np.any(injections, axis=(0, 3))
    return {
        "survival": survival,
        "native_cells": board.cells,
        "native_incidence": incidence.tolist(),
        "cube_side": WIDTH,
        "injection_sites": [list(x) for x in SHIFTS],
        "off_source_sky_columns": int(np.sum(~injecting)),
        "maximum_balance_error": error,
        "minimum_projected_laplacian_off_source": float(curvature[~injecting].min()),
        "minimum_projected_laplacian_on_sources": float(curvature[injecting].min()),
        "minimum_after_reversing_modular_sign_off_source": float(
            (-curvature)[~injecting].min()
        ),
    }


def contrast(signal, cov, i, j, k, weight):
    # Excess of endpoint chord over the middle datum; concavity permits <= 0.
    row = np.zeros(len(signal))
    row[[i, j, k]] = [weight, -1, 1 - weight]
    effect = float(row @ signal)
    sigma = float(np.sqrt(row @ cov @ row))
    return {"chord_excess": effect, "standard_error": sigma, "z": effect / sigma}


def shape_test(accel, lo, hi, signal, cov):
    triples = []
    for i, j, k in combinations(range(len(signal)), 3):

        def weight(a, b, c):
            return (c - b) / (c - a)

        lam = weight(accel[i], accel[j], accel[k])
        min_w = weight(lo[i], hi[j], lo[k])
        max_w = weight(hi[i], lo[j], hi[k])
        np.testing.assert_allclose(
            np.clip([min_w, max_w], 0, 1), [min_w, max_w], atol=1e-13
        )
        center = contrast(signal, cov, i, j, k, lam)
        uncertain = [contrast(signal, cov, i, j, k, x) for x in (min_w, max_w)]
        triples.append(
            {
                "bins": [i, j, k],
                "centre_weight": lam,
                "edge_weights": [min_w, max_w],
                "centre": center,
                "edge_contrasts": uncertain,
                "z_min_at_edge_weights": min(v["z"] for v in uncertain),
            }
        )
    n = len(triples)

    def summarize(zs):
        index = int(np.argmax(zs))
        z = zs[index]
        return {
            "triple": triples[index]["bins"],
            "largest_chord_z": z,
            "one_sided_marginal_p": float(ndtr(-z)),
            "bonferroni_p_upper": float(min(1, n * ndtr(-z))),
        }

    return {
        "simultaneous_95pct_z_threshold": float(ndtri(1 - 0.05 / n)),
        "triples": triples,
        "bin_centres": summarize([v["centre"]["z"] for v in triples]),
        "arbitrary_within_bin_weights": summarize(
            [v["z_min_at_edge_weights"] for v in triples]
        ),
    }


def read_data(archive):
    if sha256(archive.read_bytes()).hexdigest() != lensing_profile.ARCHIVE_SHA256:
        raise ValueError("Unexpected KiDS release")
    data = {}
    with tarfile.open(archive) as f:
        for key in ("KiDS-bright", "GAMA"):
            panel = lensing_profile.load_profile(f, lensing_profile.PROFILES[key])
            keep = np.array(panel["lower_bin_edge_m_s2"]) >= 1e-13
            for name in (
                "g_bar_m_s2",
                "lower_bin_edge_m_s2",
                "tangential_esd",
                "cross_esd",
            ):
                panel[name] = np.asarray(panel[name])[keep].tolist()
            panel["covariance"] = np.array(panel["covariance"])[
                np.ix_(keep, keep)
            ].tolist()
            if np.sum(keep) != 6:
                raise ValueError("Different set of enrolled shear bins")
            data[key] = panel
    return data


def run(data):
    survey = {}
    for key, panel in data.items():
        x = np.array(panel["g_bar_m_s2"])
        lo = np.array(panel["lower_bin_edge_m_s2"])
        hi = 2 * x - lo
        cov = np.array(panel["covariance"])
        survey[key] = {
            name: shape_test(x, lo, hi, np.array(panel[name]), cov)
            for name in ("tangential_esd", "cross_esd")
        }
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "source_url": lensing_profile.ARCHIVE,
        "source_sha256": lensing_profile.ARCHIVE_SHA256,
        "source_data": data,
        "cubic_walk": [finite_walk(s) for s in (1.0, 7 / 8)],
        "measurements": survey,
    }


def verify(saved, actual):
    # The residual of an exact identity can change by several ulps across FFTs.
    # Keep the equation's absolute tolerance and compare all other measurements.
    def checked_values(result):
        walks = []
        for row in result["cubic_walk"]:
            values = row.copy()
            if not 0 <= values.pop("maximum_balance_error") < 1e-11:
                raise ArithmeticError("Squared-visit residual exceeds tolerance")
            walks.append(values)
        return {**result, "cubic_walk": walks}

    lensing_profile.check(checked_values(saved), checked_values(actual))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--archive", type=Path, help="authors' original Brouwer et al. archive"
    )
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    if args.write and not args.archive:
        p.error("The original archive is required to create a receipt")
    data = (
        read_data(args.archive)
        if args.archive
        else json.loads(RECEIPT.read_text())["source_data"]
    )
    out = run(data)
    if args.write:
        with RECEIPT.open("x") as f:
            json.dump(out, f, allow_nan=False, indent=2)
            f.write("\n")
    else:
        verify(json.loads(RECEIPT.read_text()), out)
    print(
        json.dumps(
            {
                k: {
                    ch: (s["bin_centres"], s["arbitrary_within_bin_weights"])
                    for ch, s in r.items()
                }
                for k, r in out["measurements"].items()
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
