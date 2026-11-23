"""Decode ring measurements from the pinned public vector figures for AGC 114905."""

import argparse
import hashlib
import io
import json
import re
import tarfile
from pathlib import Path

import numpy as np
from pypdf import PdfReader

HERE = Path(__file__).resolve().parent
SOURCE_HASH = "c921b5f148d5484811d5d1dbd6f5ae149bcf33feb6e0360c27775afd2720a337"
KPC = 3.085677581491367e19
# Axis calibrations are the origins and annotated tick positions in PDF points.
PANELS = {
    "figure7": (
        "plot_mm_noDM.pdf",
        "f903f6adf6a8fa505971712608cd9c0fe6fc716afaeff6d0963611c58115aada",
        (59.904, 131.956364, 2),
        (47.3472, 104.038892, 10),
    ),
    "figure8": (
        "plot_mond.pdf",
        "87f3734ea0d74fa891e207b5978199c63b852cec7fd092d2c8af416efd1ad132",
        (58.5, 132.103636, 2),
        (48.6, 115.513043, 20),
    ),
}


def read_panel(pdf, x_axis, y_axis):
    source = (
        PdfReader(io.BytesIO(pdf)).pages[0].get_contents().get_data().decode("latin1")
    )

    def calibrate(values, axis):
        zero, other, units = axis
        return (values - zero) * units / (other - zero)

    # Each scientific data line is clipped to the plotting rectangle. Legends
    # and axes have their own paths outside these clipped groups.
    sections = re.split(r"Q?\s*q\s+[\d.]+ [\d.]+ [\d.]+ [\d.]+ re W n ", source)[1:]
    grey = []
    red = []
    for section in sections:
        start, _, _ = section.partition(" cs\n")
        if len(start) > 200:
            continue
        pairs = re.findall(r"(?m)^([\d.]+) ([\d.]+) ([ml])$", section.split("\nQ")[0])
        points = np.array([[float(x), float(y)] for x, y, _ in pairs])
        if len(points) == 10 and "1.5 w 0.5019607843 G" in start:
            grey.append(points)
        if len(points) > 10 and "1 0 0 RG" in start:
            red.append((points, "[" in start))
    centres, errors = grey[:2]
    if not (
        np.allclose(centres[::2], centres[1::2])
        and np.allclose(centres[::2], (errors[::2] + errors[1::2]) / 2, atol=1e-5)
    ):
        raise ValueError("Unexpected grey marker or uncertainty coordinates")
    x, y = centres[::2].T
    half = calibrate(errors[1::2, 1], y_axis) - calibrate(y, y_axis)
    radius = calibrate(x, x_axis)
    velocities = [calibrate(np.interp(x, p[:, 0], p[:, 1]), y_axis) for p, _ in red]
    return radius, calibrate(y, y_axis), half, velocities, [dashed for _, dashed in red]


def baryons_from_published_interpolation(modeled, radius, scale):
    # Eq. 11 of the observational paper has Vmodel^4 = Vbar²(Vmodel² + a0 R).
    return modeled**2 / np.sqrt(modeled**2 + scale * radius * KPC / 1e6)


def extract(source_tar):
    if hashlib.sha256(source_tar.read_bytes()).hexdigest() != SOURCE_HASH:
        raise ValueError("Unrecognized author source archive")
    panels = {}
    with tarfile.open(source_tar) as archive:
        for label, (filename, digest, xaxis, yaxis) in PANELS.items():
            pdf = archive.extractfile(filename).read()
            if hashlib.sha256(pdf).hexdigest() != digest:
                raise ValueError(f"{filename} differs from registered source")
            r, v, e, model, dashed = read_panel(pdf, xaxis, yaxis)
            if label == "figure8":
                assert dashed == [True, False]
                alt = baryons_from_published_interpolation(model[0], r, 2e-12)
                vbar = baryons_from_published_interpolation(model[1], r, 1.2e-10)
            else:
                assert dashed == [False]
                vbar = model[0]
                alt = None
            rows = []
            for i in range(len(r)):
                row = dict(
                    radius_kpc=float(r[i]),
                    circular_kms=float(v[i]),
                    plotted_error_kms=float(e[i]),
                    baryonic_kms=float(vbar[i]),
                )
                if alt is not None:
                    row["baryonic_from_second_curve_kms"] = float(alt[i])
                rows.append(row)
            panels[label] = dict(
                source_filename=filename, source_sha256=digest, rings=rows
            )
    return dict(source_sha256=SOURCE_HASH, panels=panels)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_archive", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = extract(args.source_archive)
    path = HERE / "rings.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    elif data != json.loads(path.read_text()):
        raise SystemExit("Vector figure extraction differs from saved rings")
    print("Reconstructed five paired rings in each panel.")
