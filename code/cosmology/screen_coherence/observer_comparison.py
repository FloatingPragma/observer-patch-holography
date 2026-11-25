"""Read equilibrium, radar and oscillator measurements from archived detector words."""

import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
COSMOLOGY = HERE.parent


def read_json(path):
    return json.loads(path.read_text())


def thermal_readings():
    directory = COSMOLOGY / "guarded_slab/energy_exchange"
    banks = read_json(directory / "precision.json")
    starts = read_json(directory / "precision_protocol.json")["launches_per_bank"]
    # The first bit in each probe/reservoir pair is the probe's energy reading.
    probes = (np.arange(64)[:, None] >> [0, 2, 4]) & 1
    all_agree = probes.min(axis=1) == probes.max(axis=1)
    result = {}
    for name, bank in banks.items():
        words = np.asarray(bank["age_by_terminal_word"]).sum(axis=0)
        delivered = int(words.sum())
        refused = sum(map(sum, bank["first_filter_leaks"]))
        assert delivered == bank["delivered"]
        assert delivered + refused + bank["timeout"] == starts
        result[name] = dict(
            starts=starts,
            delivered=delivered,
            refused=refused,
            timeout=bank["timeout"],
            upper_given_delivery=(words @ probes / delivered).tolist(),
            all_agree_given_delivery=float(words[all_agree].sum() / delivered),
        )
    return result


def aperture_readings():
    with np.load(HERE / "survey_counts.npz") as archive:
        counts = archive["L2_counts"]
        blank = archive["L2_blank"]
    # The two settled schedules give duplicate observations for each preparation.
    np.testing.assert_array_equal(counts[:, 0, -1], counts[:, 1, -1])
    assert np.all((counts[:, 0, -1] > 0) & (counts[:, 0, -1] <= blank))
    contrast = np.log(counts[:, 0, -1] / blank) / 3
    # Center across preparations for each individual observer and aperture.
    centered = contrast - contrast.mean(axis=0)
    variance = contrast.var(axis=0, ddof=1).mean(axis=1)
    four, six = centered.transpose(1, 0, 2)
    correlation = np.corrcoef(four.ravel(), six.ravel())[0, 1]
    return dict(
        preparations=len(contrast),
        observers=contrast.shape[-1],
        variance_4_and_6=variance.tolist(),
        six_to_four=float(variance[1] / variance[0]),
        paired_aperture_correlation=float(correlation),
    )


def clock_readings():
    directory = COSMOLOGY / "release_horizon/transverse_clock"
    clock = read_json(directory / "records.json")["oscillator"]
    assert clock["marker"] == [1, 0, 1, 0]
    result = {}
    eps = clock["signed_halfwidth"]
    for mode in clock["modes"]:
        if mode["mode"] not in ("port", "oblique"):
            continue
        signed = {}
        for row in mode["reads"]:
            counts = np.asarray(row["counts"])
            assert counts.shape == (2, 2, 3) and (counts >= 0).all()
            assert (counts.sum(axis=-1) == clock["shots"]).all()
            signed[row["input"], row["t"]] = (
                counts[:, :, 0] - counts[:, :, 1]
            ) / clock["shots"]
        cosine = signed["vertex", 1][0, 0]
        errors = {}
        for prep in ("vertex", "edge_plus", "edge_minus"):
            # Two marker readings at receiver 0 predict the third at receiver 1.
            prediction = 2 * cosine * signed[prep, 2][0] - signed[prep, 1][0]
            errors[prep] = (signed[prep, 3][1] - prediction).tolist()
        # A nonzero determinant excludes even a signed rescaling of one time profile.
        a, b = (signed["vertex", t][1, 0] for t in (1, 2))
        c, d = (signed["edge_plus", t][1, 0] for t in (1, 2))
        result[mode["mode"]] = dict(
            calibration_cosine=cosine,
            third_mark_residuals=errors,
            max_abs_residual=float(np.abs(list(errors.values())).max()),
            recurrence_sampling_allowance=6 * eps,
            two_preparation_determinant=a * d - b * c,
            determinant_sampling_allowance=4 * eps,
        )
    return result


if __name__ == "__main__":
    out = dict(
        thermometers=thermal_readings(),
        apertures=aperture_readings(),
        clocks=clock_readings(),
    )
    print(json.dumps(out, indent=2, allow_nan=False))
