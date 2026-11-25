"""Verify the original acquisition exactly, then recompute its physical readout law."""

from hashlib import sha256
import json
from math import log, sqrt

import numpy as np

from . import measure

RESULT_SHA256 = "0c87fca5662b37a7e968b2f5e4bedec0bcc2c4743c7f71db26a6b3b0f9c1a112"
POPULATIONS_SHA256 = "e9ebeb92255098d41f66d8f75e18858665bcee3a7f0b4aaf5c485c230001ff0b"


def archived(name, digest):
    data = (measure.HERE / name).read_bytes().replace(b"\r\n", b"\n")
    if sha256(data).hexdigest() != digest:
        raise ValueError(f"Changed impulse acquisition {name}")
    return json.loads(data)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def resolution(p, shots, log_family):
    """Bernstein radius for a multinomial bin's Bernoulli marginal."""
    return np.sqrt(2 * p * (1 - p) * log_family / shots) + 2 * log_family / (3 * shots)


def close(a, b):
    require(abs(a - b) <= 3e-12, "Inconsistent derived impulse readout")


def shot_record(arm, acquisition, log_family):
    seed, admission = acquisition["seed"], acquisition["admission"]
    probability = np.array(acquisition["probabilities"])
    values = np.array(arm["values"])
    require(
        len(probability) == len(values) + 1 and np.all(probability >= 0),
        "Incomplete acquisition bins",
    )
    close(probability.sum(), 1)
    require(arm["enrolled"] == measure.SHOTS, "Changed trial enrollment")
    rng = np.random.default_rng(seed)
    selected = int(rng.binomial(measure.SHOTS, admission))
    counts = rng.multinomial(selected, probability).tolist()
    require(
        arm["failed_preparation"] == measure.SHOTS - selected
        and [*arm["counts"], arm["failed_band_filter"]] == counts,
        "Seeded impulse counts differ",
    )
    close(arm["mean"], np.dot(counts[:-1], values) / selected)
    close(arm["exact_mean"], np.dot(probability[:-1], values))
    close(
        arm["halfwidth_99_simultaneous"],
        sqrt(2 * log(2 * 24 / 0.01) / selected) + 1e-12,
    )
    # Count custody and the independent population comparison are separate checks.
    return selected, resolution(probability, selected, log_family)


def verify():
    saved = archived("results.json", RESULT_SHA256)
    scheme = archived("acquisition.json", POPULATIONS_SHA256)
    bank = []
    fresh = measure.run(bank)
    require(len(scheme) == len(bank) == 24, "Incomplete source probability table")
    # Count every admission, terminal address and failed band filter in the
    # comparison family. This uses a source-resolution test, not
    # a claim that two deterministic eigensolvers are random samplings.
    log_family = log(2 * sum(1 + len(s["probabilities"]) for s in scheme) / 0.01)
    require(
        all(
            fresh[k] == saved[k]
            for k in ("theta", "native_energy_units", "spatial_units")
        ),
        "Changed instrument units",
    )
    index = 0
    for src, ref in zip(fresh["preparations"], saved["preparations"], strict=True):
        require(
            (src["owner"], src["axis"]) == (ref["owner"], ref["axis"]),
            "Changed preparation axis",
        )
        for k in ("native", "contrast_removed"):
            # Rest projection of a sharp source has amplitude sqrt(admission).
            # Projector residuals bound its two normalized preparations; a
            # projector-event expectation moves by <=2 times their vector distance.
            scale = (
                4
                * sum(v["rest_projector_error"] for v in (fresh, saved))
                / sqrt(min(src["rest_admission"], ref["rest_admission"]))
            )
            require(
                abs(src["controls"][k] - ref["controls"][k]) <= scale,
                "Changed rest-packet control",
            )
        close(
            src["controls"]["uniform_common_response"],
            ref["controls"]["uniform_common_response"],
        )
        for row in (src, ref):
            require(
                row["controls"]["native"] > 0.99
                and row["controls"]["contrast_removed"] < 0.32
                and row["controls"]["uniform_common_response"] > 0.01,
                "Failed recoil control",
            )
        for number, (label, sign) in enumerate((("common", -1), ("signed_load", 1))):
            here, previous = src["cases"][label], ref["cases"][label]
            for tag, orientation in (("all_ports", sign), ("within_cluster", 1)):
                a, b = (
                    here["derivative_intervals"][tag],
                    previous["derivative_intervals"][tag],
                )
                require(
                    max(a[0], b[0]) <= min(a[1], b[1])
                    and min(orientation * v for v in a + b) > 0,
                    "Incompatible signed spectral enclosures",
                )
            for label_time, t, orientation in (
                ("free", 1 / 16, sign),
                ("band_selected", 32768, 1),
            ):
                panel, old = here["finite"][label_time], previous["finite"][label_time]
                require(
                    panel["time"] == old["time"] == t, "Changed impulse appointment"
                )
                for i, (now, then) in enumerate(
                    zip(panel["records"], old["records"], strict=True)
                ):
                    actual, frozen = bank[index], scheme[index]
                    require(
                        actual["seed"]
                        == frozen["seed"]
                        == [174920, src["owner"], number, int(t), i],
                        "Changed source/receiver registration",
                    )
                    require(
                        actual["admission"] == src["rest_admission"]
                        and frozen["admission"] == ref["rest_admission"],
                        "Changed admission registration",
                    )
                    count, bounds = shot_record(then, frozen, log_family)
                    shot_record(now, actual, log_family)
                    require(
                        count > 0
                        and abs(actual["admission"] - frozen["admission"])
                        <= resolution(frozen["admission"], measure.SHOTS, log_family),
                        "Rest admission differs above counting resolution",
                    )
                    pa, pb = (
                        np.array(actual["probabilities"]),
                        np.array(frozen["probabilities"]),
                    )
                    require(
                        now["values"] == then["values"]
                        and pa.shape == pb.shape
                        and np.all(np.abs(pa - pb) <= bounds),
                        "Recomputed terminal probability differs above counting resolution",
                    )
                    index += 1
                for p in (panel, old):
                    low, high = p["records"]
                    unit = 2 * measure.THETA * t
                    close(
                        p["signed_position_response"],
                        (high["mean"] - low["mean"]) / unit,
                    )
                    close(
                        p["nominal_response"],
                        (high["exact_mean"] - low["exact_mean"]) / unit,
                    )
                    close(
                        p["halfwidth_99_simultaneous"],
                        (
                            high["halfwidth_99_simultaneous"]
                            + low["halfwidth_99_simultaneous"]
                        )
                        / unit,
                    )
                    require(
                        abs(p["signed_position_response"] - p["nominal_response"])
                        < p["halfwidth_99_simultaneous"],
                        "Nominal recoil misses its acquisition interval",
                    )
                    require(
                        orientation * p["signed_position_response"]
                        > p["halfwidth_99_simultaneous"],
                        "Signed recoil gate fails",
                    )
                require(
                    abs(panel["nominal_response"] - old["signed_position_response"])
                    < old["halfwidth_99_simultaneous"]
                    and abs(old["nominal_response"] - panel["signed_position_response"])
                    < panel["halfwidth_99_simultaneous"],
                    "Independent recoil forecast misses measured interval",
                )
    require(index == len(bank), "Unchecked acquisition records")
    return fresh
