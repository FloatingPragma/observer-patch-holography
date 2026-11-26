"""Centroid impulse, current and body-centred material probes on the row-14 sources."""

import gzip
import json
from pathlib import Path

import numpy as np

from .. import measure as stock
from ..scattering.measure import detector, gain
from ...row_16_binding.signed_pair.material import vertices

HERE = Path(__file__).resolve().parent
STARTS = 2**44
SEED = 1401624
BUDGET = 4096
KICKS = (0.0, -0.25, 0.25, -0.5, 0.5)


def radius(scores, n):
    return float(np.ptp(scores) * np.sqrt(np.log(2 * BUDGET / 0.01) / (2 * n)))


def sample(rng, probabilities, acceptance):
    p = np.asarray(probabilities, float)
    if np.min(p) < -2e-12 or abs(np.sum(p) - 1) > 1e-10:
        raise ValueError("Incomplete detector outcomes")
    p = np.clip(p, 0, 1)
    accepted = rng.binomial(STARTS, acceptance)
    draw = rng.multinomial(accepted, p / p.sum())
    return dict(counts=draw.tolist(), refused=STARTS - accepted)


def read(bank, scores):
    counts = np.asarray(bank["counts"])
    return float(counts @ scores / counts.sum()), radius(scores, counts.sum())


def inertia(rng, rho, swaps, coord, acceptance):
    """Random seam parity and a separate imaginary quadrature of the same pair."""
    idx = np.arange(len(rho))
    diagonal = np.diag(rho)
    # Complete both orientations of every unordered pair. For fixed words the
    # current vanishes, and the swap parity is deterministically even.
    coh = rho[idx, swaps]
    pop = np.broadcast_to(diagonal, swaps.shape)
    step = coord[swaps] - coord
    magnitudes, lab = np.unique(np.round(abs(step), 14), return_inverse=True)
    lab = lab.reshape(step.shape)
    pm = np.zeros((len(magnitudes), 2))
    scores = np.stack([15 * magnitudes, -15 * magnitudes], axis=1).ravel()
    banks = {}
    for u in KICKS:
        odd = (1 - (coh * np.cos(u * step)).sum() / 30) / 2
        e = sample(rng, [1 - odd, odd], acceptance)
        y = coh * np.sign(step) * np.sin(u * step)
        for sign, side in ((1, 0), (-1, 1)):
            pm[:, side] = np.bincount(
                lab.ravel(),
                weights=((pop + sign * y) / 60).ravel(),
                minlength=len(magnitudes),
            )
        flow = sample(rng, pm.ravel(), acceptance)
        banks[str(u)] = {"energy": e, "current": flow}
    # For every density matrix, 2|rho_ij| <= rho_ii + rho_jj. Sum over
    # both orientations of the involution to replace the coherence by p_i.
    fourth_ceiling = float(np.max((step**4).sum(axis=0)))
    return {
        "current_levels": scores.tolist(),
        "fourth_ceiling": fourth_ceiling,
        "banks": banks,
    }


def recoil_decode(axes, u):
    estimates = []
    for axis in axes:
        banks = axis["banks"]
        eo, de = read(banks["0.0"]["energy"], [0, 30])
        ep, se = read(banks[str(u)]["energy"], [0, 30])
        em, te = read(banks[str(-u)]["energy"], [0, 30])
        vp, sv = read(banks[str(u)]["current"], axis["current_levels"])
        vm, tv = read(banks[str(-u)]["current"], axis["current_levels"])
        b = axis["fourth_ceiling"] * u**2
        estimates.append(
            [
                (ep + em - 2 * eo) / u**2,
                (se + te + 2 * de) / u**2 + b / 24,
                (vp - vm) / (2 * u),
                (sv + tv) / (2 * u) + b / 12,
            ]
        )
    return np.mean(estimates, axis=0).tolist()


def radial(rng, words, p, xyz, flags, acceptance, g):
    occ = words != 0
    centroid = occ @ xyz / occ.sum(axis=1)[:, None]
    d2 = np.sum((xyz[None, :, :] - centroid[:, None, :]) ** 2, axis=2)
    tied = np.round(d2, 12)
    scales = {"second": float(d2.max()), "nearest": 1.0, "farthest": 1.0}
    weights = {
        "second": d2 / scales["second"],
        "nearest": occ & (tied == np.where(occ, tied, np.inf).min(axis=1)[:, None]),
        "farthest": occ & (tied == np.where(occ, tied, -np.inf).max(axis=1)[:, None]),
    }
    output = {}
    for aperture, weight in weights.items():
        output[aperture] = {}
        for tag, polarity in (("forward", 1), ("inverted", -1), ("host", 0)):
            blocks = []
            for flag in (False, True):
                subset = flags == flag
                weighted = p[subset, None] * weight[subset] / 12
                through = weighted.sum()
                mark = words[subset] * polarity
                mean = (weighted * mark).sum() / through
                power = (weighted * mark**2).sum() / through
                blocks.extend(through * detector(mean, power, g, True))
                blocks.append(p[subset].sum() - through)
            output[aperture][tag] = sample(rng, blocks, acceptance)
    return dict(scales=scales, banks=output)


def radial_read(bank, scale, calibration, connected):
    """Divide the raw balanced fringe by the independently measured unit gain."""
    counts = np.asarray(bank["counts"]).reshape(2, 4)
    if connected is None:
        counts = counts.sum(axis=0)
    else:
        counts = counts[int(connected)]
    t, e = read({"counts": counts}, [1, -1, 0, 0])
    g, err = read(calibration, [1, -1, 0])
    if g <= err:
        return None
    value = 12 * scale * t / g
    return [value, (12 * scale * e + abs(value) * err) / (g - err)]


def run():
    rng = np.random.default_rng(SEED)
    mask, _ = stock.native_geometry()
    xyz = vertices()
    result = {"starts": STARTS, "seed": SEED, "calibration": {}, "panels": {}}
    for tag, k, _ in stock.PANELS:
        result["calibration"][tag] = {
            str(q): sample(rng, detector(q, q * q, gain(k), True), 1)
            for q in (-1, 0, 1)
        }
    for name in ("neutral", "positive"):
        w, swaps = stock.sector(*stock.NAMES[name])
        _, _, flags = stock.graph_read(w)
        flags = np.asarray(flags)
        geom = (w != 0) @ xyz / (w[0] != 0).sum()
        kin, co, ref = stock.operators(w, swaps, mask)
        for tag, k, enable in stock.PANELS:
            rho, spec = stock.filtered(kin, co, ref, k, enable)
            acceptance = spec["acceptance"]
            axes = [inertia(rng, rho, swaps, x, acceptance) for x in geom.T]
            p = np.diag(rho)
            profile = radial(rng, w, p, xyz, flags, acceptance, gain(k))
            row = {
                "acceptance": acceptance,
                "source_clock": spec["target"],
                "axes": axes,
                "recoil": {str(u): recoil_decode(axes, u) for u in (0.25, 0.5)},
                "radial": profile,
            }
            if name == "neutral":
                # Following the unsigned support read, swap each positive mark
                # with one vacancy. This leaves the target centroid unchanged.
                row["conjugated"] = radial(rng, -w, p, xyz, flags, acceptance, gain(k))
            result["panels"][name + "/" + tag] = row
    return result


def summary(result):
    for key, row in result["panels"].items():
        print(key, "even work and odd velocity", row["recoil"])
        for subset in (None, True):
            prof = row["radial"]
            print(
                "all" if subset is None else "connected",
                {
                    name: radial_read(
                        bank["forward"],
                        prof["scales"][name],
                        result["calibration"][key.split("/")[1]]["1"],
                        subset,
                    )
                    for name, bank in prof["banks"].items()
                },
            )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    parsed = parser.parse_args()
    result = run()
    if parsed.record:
        with parsed.record.open("xb") as f:
            f.write(
                gzip.compress(
                    json.dumps(result, separators=(",", ":")).encode(), mtime=0
                )
            )
    else:
        from ..scattering.measure import verify

        saved = json.loads(gzip.decompress((HERE / "counts.json.gz").read_bytes()))
        verify(result, saved)
    summary(result)
