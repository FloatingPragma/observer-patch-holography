"""Project the round midpoint source onto the two released tangential shear profiles."""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from math import sqrt
from pathlib import Path

import numpy as np
from scipy.stats import norm

from .. import lensing_profile
from ..record_shear import measure as shear

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "MIDPOINT_SHEAR.md"
INPUT = shear.HERE / "results.json"


def esd(radius, edge):
    """Excess surface density in units of twice the 1/r coefficient."""
    if radius >= edge:
        return (edge / radius) ** 2
    return 1 / (1 + sqrt(1 - (radius / edge) ** 2))


def bounds(squared_ratio):
    return squared_ratio, 1 + sqrt(1 - squared_ratio)


def contrast(panel, channel, far, near, weights):
    y = np.asarray(panel[channel])[[far, near]]
    cov = np.asarray(panel["covariance"])[np.ix_([far, near], [far, near])]
    w = np.array(weights)
    value = float(w @ y)
    error = sqrt(w @ cov @ w)
    return {"value": value, "error": error, "z": value / error}


def compare(panel, channel):
    centres = panel["g_bar_m_s2"]
    low = panel["lower_bin_edge_m_s2"]
    pairs = []
    for far, near in combinations(range(len(centres)), 2):
        hi_near = 2 * centres[near] - low[near]
        readings = {}
        for tag, ratio in (
            ("centre", centres[far] / centres[near]),
            ("full_bins", low[far] / hi_near),
        ):
            lower, upper = bounds(ratio)
            readings[tag] = {
                "squared_radius_ratio": ratio,
                "lower": contrast(panel, channel, far, near, [-1, lower]),
                "upper": contrast(panel, channel, far, near, [1, -upper]),
            }
        pairs.append({"bins": [far, near], **readings})
    maxima = {}
    for tag in ("centre", "full_bins"):
        z, pair, side = max(
            (p[tag][side]["z"], p["bins"], side)
            for p in pairs
            for side in ("lower", "upper")
        )
        maxima[tag] = {
            "z": z,
            "bins": pair,
            "bound": side,
            "family_tail_upper": min(1.0, 2 * len(pairs) * float(norm.sf(z))),
        }
    return {"pairs": pairs, "maxima": maxima}


def run(archive=None):
    saved = json.loads(INPUT.read_bytes())
    data = saved["source_data"]
    if archive is not None:
        assert shear.read_data(archive) == data
    return {
        "source_url": saved["source_url"],
        "source_sha256": saved["source_sha256"],
        "binned_input_sha256": sha256(INPUT.read_bytes()).hexdigest(),
        "protocol_sha256": sha256(PROTOCOL.read_bytes()).hexdigest(),
        "one_sided_family_threshold": float(norm.isf(0.05 / 30)),
        "surveys": {
            name: {
                channel: compare(panel, channel)
                for channel in ("tangential_esd", "cross_esd")
            }
            for name, panel in data.items()
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--archive", type=Path)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run(args.archive)
    receipt = HERE / "midpoint_shear.json"
    if args.write:
        receipt.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    else:
        lensing_profile.check(json.loads(receipt.read_bytes()), result)
    print({s: p["tangential_esd"]["maxima"] for s, p in result["surveys"].items()})
