"""Test a family-independent spatial floor as a correction to charged guide rest lines."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

import mpmath as mp
import networkx as nx
import numpy as np

from model.executable_model.carrier import federation, seam_digest
from ..dirac_read.measure import register_effect
from .search import blocks, certify_positive

HERE = Path(__file__).resolve().parent


def run():
    comparison = json.loads(
        (HERE.parent / "arrival_spin/anchor_spectrum.json").read_text()
    )
    reference = comparison["scale_audit"]
    # The frozen comparison supplies the target and uncertainty only at the readout stage.
    geometric = json.loads((HERE / "results.json").read_text())
    with mp.workdps(60):
        phase = mp.mpf(2) / 9
        masses = sorted(
            (1 + mp.sqrt(2) * mp.cos(phase + 2 * j * mp.pi / 3)) ** 2 / 27
            for j in range(3)
        )
        np.testing.assert_allclose(
            [float(x) for x in masses], register_effect()[1], rtol=1e-12
        )
        light, middle, heavy = masses
        ratio = mp.mpf(str(reference["codata_eq_187_ratio"]))
        sigma = mp.mpf(str(reference["ratio_standard_uncertainty"]))

        def required_floor(r):
            return (middle**2 - r**2 * light**2) / (r**2 - 1)

        floor = required_floor(ratio)
        lo, hi = required_floor(ratio + 5 * sigma), required_floor(ratio - 5 * sigma)
        assert 0 < lo < floor < hi

        def ratio_at(x, mass):
            return mp.sqrt((x + mass**2) / (x + light**2))

        assert abs(ratio_at(floor, middle) - ratio) < mp.mpf("1e-50")
        target = {
            "zero_floor_middle_ratio": float(middle / light),
            "frozen_middle_ratio": float(ratio),
            "frozen_standard_uncertainty": float(sigma),
            "lambda_required_at_central_ratio": float(floor),
            "lambda_range_at_five_quoted_uncertainties": [float(lo), float(hi)],
            "heavy_light_ratio_conditional_on_fitted_floor": float(
                ratio_at(floor, heavy)
            ),
            "lambda_units": "squared exposure energy, declared register coupling one",
        }
        by_board = []
        for level, saved in enumerate(geometric["federations"]):
            board = federation(level)
            assert seam_digest(board) == saved["seam_sha256"]
            high = nx.Graph()
            high.add_nodes_from(np.flatnonzero(board.degree == 6))
            high.add_edges_from(
                (u, v)
                for u, v in board.seams
                if board.degree[u] == board.degree[v] == 6
            )
            assert nx.is_forest(high) and max(dict(high.degree).values()) <= 2
            limit = max(len(c) for c in nx.connected_components(high))
            lower = Fraction(
                saved["pinned_connected_subgraphs"]["certified_strict_lower_bound"]
            )
            forms = list(blocks(limit, "branched"))
            for _, (stiffness, norm) in forms:
                certify_positive(stiffness, norm, lower)
            assert (
                limit == saved["largest_consecutive_unpinned_run_bound"] and lower < 8
            )
            stiff = mp.mpf(lower.numerator) / lower.denominator
            # A connected all-degree-six scalar guide is a path. Every nonzero
            # standing mode lies above its longest possible path's first mode.
            first_mode = 16 * mp.sin(mp.pi / (2 * limit)) ** 2
            assert hi < stiff < first_mode
            by_board.append(
                {
                    "cells": board.cells,
                    "longest_degree_six_path": limit,
                    "all_degree_six_ground": 0,
                    "least_possible_positive_path_mode": float(first_mode),
                    "pinned_guide_strict_floor": str(lower),
                    "branch_block_certificates": len(forms),
                    "all_positive_floors_over_five_sigma_upper_bound": float(
                        stiff / hi
                    ),
                    "minimum_mass_coupling_if_fitted_on_pinned_guides": float(
                        mp.sqrt(stiff / hi)
                    ),
                    "upper_middle_ratio_on_pinned_guides": float(
                        ratio_at(stiff, middle)
                    ),
                }
            )
        return {
            "frozen_ratio_source": reference["ratio_source"],
            "comparison": target,
            "boards": by_board,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    path = HERE / "floor_calibration.json"
    result = run()
    if parser.parse_args().write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    elif result != json.loads(path.read_text()):
        raise ValueError(
            "Measured geometric floors differ from the recorded comparison"
        )
    print(json.dumps(result["comparison"], indent=2))


if __name__ == "__main__":
    main()
