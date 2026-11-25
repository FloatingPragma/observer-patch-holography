"""Retrospective acoustic-feature comparison at Planck PR3 bin centres."""

import argparse
import hashlib
import json
from pathlib import Path

import camb
import numpy as np

from likelihood_slice import EDGE_NS


BASELINE_NS = 0.9649
BACKGROUND = dict(H0=67.36, ombh2=0.02237, omch2=0.12, mnu=0.06, omk=0.0, tau=0.0544)
AS = 2.1e-9
LMAX = 2600
DATA = {
    "TT": ("R3.01", "ed2b8b51b033e668e661f6cf4eb7705f15bcaf358189dfc4a4481cc7020130cd"),
    "TE": ("R3.02", "191a588f9000885e931739ddb4f88717c29f79a45abd337837dac8032aed0137"),
    "EE": ("R3.02", "bd30fc1121c0d1c2e1a6ac997b30ac14ea1fdd42bc91ecc40d594c00463ffdb2"),
}
ARCHIVE = (
    "https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/"
)


def transfer(ns):
    p = camb.CAMBparams()
    helium = np.asarray(camb.bbn.get_predictor().Y_He(BACKGROUND["ombh2"], 0.0)).item()
    p.set_cosmology(**BACKGROUND, YHe=helium)
    p.InitPower.set_params(As=AS, ns=ns)
    p.set_for_lmax(LMAX, lens_potential_accuracy=1)
    return camb.get_results(p).get_total_cls(lmax=LMAX, CMB_unit="muK")


def feature_rows(ell, observed, errors, baseline, branch, sign):
    # Define the search windows using the smooth control alone. Positive TE
    # crests and negative TE troughs are distinct series with signed heights.
    reference = sign * baseline
    is_peak = (reference[1:-1] > reference[:-2]) & (reference[1:-1] > reference[2:])
    centers = (np.flatnonzero(is_peak & (reference[1:-1] > 0)) + 1)[:5]
    if len(centers) != 5:
        raise ValueError("Five resolved acoustic features required in this channel")
    rows = []
    for center in centers:
        window = np.flatnonzero(np.abs(ell - ell[center]) <= 60)
        row = {}
        for label, curve in (
            ("measured", observed),
            ("control", baseline),
            ("branch", branch),
        ):
            index = window[np.argmax(sign * curve[window])]
            row[label] = [round(float(ell[index]), 5), round(float(curve[index]), 6)]
            if label == "measured":
                row["marginal_error_minus_plus"] = errors[:, index].tolist()
        rows.append(row)
    return rows


def compare(data_dir):
    measured = {}
    sources = {}
    for channel, (version, expected) in DATA.items():
        name = f"COM_PowerSpect_CMB-{channel}-binned_{version}.txt"
        path = data_dir / name
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Unexpected Planck file content for {name}")
        measured[channel] = np.loadtxt(path, usecols=(0, 1, 2, 3), unpack=True)
        sources[channel] = dict(url=ARCHIVE + name, sha256=expected)
    control, branch = transfer(BASELINE_NS), transfer(EDGE_NS)
    features = {}
    for channel, column in (("TT", 0), ("TE", 3), ("EE", 1)):
        ell, observed = measured[channel][:2]
        errors = measured[channel][2:]
        curves = [
            np.interp(ell, np.arange(len(curve)), curve[:, column])
            for curve in (control, branch)
        ]
        features[channel + "_max"] = feature_rows(ell, observed, errors, *curves, 1)
        if channel == "TE":
            features["TE_min"] = feature_rows(ell, observed, errors, *curves, -1)
    return dict(
        camb=camb.__version__,
        background=BACKGROUND,
        As=AS,
        lmax=LMAX,
        control_ns=BASELINE_NS,
        branch_ns=EDGE_NS,
        planck=sources,
        columns=["effective_multipole", "D_ell_microkelvin_squared"],
        feature_window_halfwidth=60,
        features=features,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "data_dir", type=Path, help="Directory with the three PR3 binned text files"
    )
    args = parser.parse_args()
    print(json.dumps(compare(args.data_dir), indent=2, allow_nan=False))
