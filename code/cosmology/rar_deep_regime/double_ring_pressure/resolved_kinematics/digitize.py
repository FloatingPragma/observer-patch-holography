"""Read the black observed second-moment markers from the original Fig. 12 bitmap."""

import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile

import numpy as np
from PIL import Image
from scipy.ndimage import find_objects, label, minimum_filter

HERE = Path(__file__).resolve().parent
ARCHIVE_URL = "https://arxiv.org/src/2401.08771v1"
IMAGE_NAME = "no_dm_no_bh.png"
IMAGE_SHA = "8a0a2911bc1442d29a880169783c2f161da6b95ec303dac2a5accc8936d709b0"


def decode(raw):
    assert hashlib.sha256(raw).hexdigest() == IMAGE_SHA
    image = np.asarray(Image.open(io.BytesIO(raw)).convert("RGB"))
    assert image.shape == (498, 654, 3)
    dark = image.max(axis=2) < 64
    xticks = np.where(image[439, 85:639].max(axis=-1) < 180)[0] + 85
    yticks = np.where(dark[7:440, 83])[0] + 7
    assert len(xticks) == 6 and len(yticks) == 9
    # First/last axis ticks set the scales; intermediate ticks check one-pixel rounding.
    sx = (xticks[-1] - xticks[0]) / 2.5
    sy = (yticks[-1] - yticks[0]) / 80
    assert np.max(np.abs(xticks - xticks[0] - sx * np.arange(6) / 2)) <= 1
    assert np.max(np.abs(yticks - yticks[0] - sy * np.arange(9) * 10)) <= 1

    def radius(x):
        return round(float((x - xticks[0]) / sx), 5)

    def velocity(y):
        return round(float(290 - (y - yticks[0]) / sy), 4)

    # A 3 x 3 black marker core survives; a one-pixel horizontal/vertical error bar does not.
    cores = minimum_filter(dark, size=3)
    cores[:75] = cores[410:] = False
    cores[:, :102] = cores[:, 540:] = False
    objects = find_objects(label(cores)[0])
    pixels = sorted(
        (
            round((s[1].start + s[1].stop - 1) / 2),
            round((s[0].start + s[0].stop - 1) / 2),
        )
        for s in objects
    )
    centre = min(pixels, key=lambda pos: radius(pos[0]))
    held = min(pixels, key=lambda pos: abs(radius(pos[0]) - 1))
    assert centre == (106, 114) and held == (319, 291)
    # Checked cap tips on the registered bin nearest 1 arcsec and the central vertical bar.
    bars = [(centre, [106, 106], [90, 138]), (held, [274, 362], [262, 319])]
    selected = []
    for (x, y), (left, right), (top, bottom) in bars:
        assert dark[top, x] and dark[bottom, x] and dark[y, left] and dark[y, right]
        selected.append(
            {
                "marker_pixel": [x, y],
                "horizontal_tip_pixels": [left, right],
                "vertical_tip_pixels": [top, bottom],
                "radial_bar_arcsec": [radius(left), radius(right)],
                "moment_bar_km_s": [velocity(bottom), velocity(top)],
            }
        )
    return {
        "source_archive": ARCHIVE_URL,
        "figure_number": 12,
        "preprint_page": 11,
        "archive_member": IMAGE_NAME,
        "figure_sha256": IMAGE_SHA,
        "tick_x_pixels": xticks.tolist(),
        "tick_y_pixels": yticks.tolist(),
        "radius_and_vrms": [[radius(x), velocity(y)] for x, y in pixels],
        "central_and_nearest_to_one_arcsec": selected,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    with tarfile.open(args.archive) as tar:
        record = decode(tar.extractfile(IMAGE_NAME).read())
    path = HERE / "muse.json"
    if args.write:
        path.write_text(json.dumps(record, indent=2) + "\n")
    else:
        assert record == json.loads(path.read_text())
        print("Reproduced the independent second-moment markers from Figure 12.")
