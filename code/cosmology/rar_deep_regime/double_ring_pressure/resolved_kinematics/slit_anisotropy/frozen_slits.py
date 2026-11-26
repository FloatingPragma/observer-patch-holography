"""Zero-adjustment predictions at all Keck apertures for the two frozen stellar fields."""

import json
from math import hypot, log, sqrt
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from .. import measure as m
from .strong_radial import smooth


def predictions(refined=False):
    data = json.loads(Path(__file__).with_name("keck.json").read_text())
    seeing = data["seeing_fwhm_arcsec"] / sqrt(8 * log(2))
    sky = np.linspace(0, 6, 2401 if refined else 1801)
    offsets = np.linspace(0, 1.5, 901 if refined else 601)
    light = smooth(sky, m.intensity(sky), offsets, seeing)
    pixels, weights = roots_legendre(20 if refined else 16)
    weights = np.outer(weights, weights).ravel() / 4
    transverse = (pixels[:, None] * 0.5 + np.zeros_like(pixels)).ravel()
    longitudinal = (np.zeros_like(pixels[:, None]) + pixels * 0.105).ravel()
    rows = data["measurements"]
    scale = m.lens.geometry()[0]
    a0 = m.lens.acceleration()[0]
    result = {}
    for name, logs in m.lens.LOG_MASSES.items():
        moment = smooth(
            sky,
            m.los_moment(sky, logs, scale, a0, None, 180 if refined else 128),
            offsets,
            seeing,
        )
        panel = []
        for offset, v, dv, dispersion, ds in rows:
            o = abs(offset)
            radii = np.hypot(transverse, longitudinal + o)
            nominal = sqrt(
                (np.interp(radii, offsets, moment) @ weights)
                / (np.interp(radii, offsets, light) @ weights)
            )
            # Each positive rectangular aperture integral decreases in either absolute centre
            # coordinate. Bound its numerator/denominator by opposite corners per centre box.
            dx = np.linspace(0, 0.03, 9)
            dy = np.linspace(max(0, o - 0.03), o + 0.03, 33)
            rr = np.hypot(transverse + dx[:, None, None], longitudinal + dy[:, None])
            j = np.interp(rr, offsets, moment) @ weights
            i = np.interp(rr, offsets, light) @ weights
            lower = np.sqrt(np.min(j[1:, 1:] / i[:-1, :-1]))
            upper = np.sqrt(np.max(j[:-1, :-1] / i[1:, 1:]))
            lo = hypot(max(abs(v) - dv, 0), dispersion - ds)
            hi = hypot(abs(v) + dv, dispersion + ds)
            panel.append(
                dict(
                    offset_arcsec=offset,
                    predicted_km_s=round(nominal, 2),
                    center_box_km_s=[
                        round(float(lower) - 0.01, 2),
                        round(float(upper) + 0.01, 2),
                    ],
                    measured_component_box_km_s=[round(lo, 2), round(hi, 2)],
                    inside=bool(lo <= nominal <= hi),
                )
            )
        result[name] = dict(inside_count=sum(r["inside"] for r in panel), panel=panel)
    return result
