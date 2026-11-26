"""Independent integer-propagation enclosures for same-source phase slopes."""

import argparse
import json

from mpmath import iv
import numpy as np

from . import measure


def contraction(left, right, coin):
    """Tr(left coin right†), with real and imaginary integer limbs."""
    a, b = left
    c, d = coin
    u, v = right
    re, im = a @ c - b @ d, a @ d + b @ c
    return int(np.sum(re * u + im * v)), int(np.sum(im * u - re * v))


def compute():
    # Only the separate exact replay needs the code/ package path and SymPy.
    from matched_einstein_limit.source_implication.directional_stress.traction.certificate import (
        four_tick_columns,
    )
    from matched_einstein_limit.source_implication.full_geometry.station_balls.formed_optics import (
        certificate as formation,
    )

    iv.dps = 45
    native = formation.native
    rho = native.formed_state()
    tau = formation.product(rho)
    coins = [formation.quantize(stock) for stock in (tau, rho)]
    board, ports, emit = measure.board_setup()
    start = [4 * x + i for x in emit for i in range(4)]
    end = [4 * x + i for x in ports for i in range(4)]
    coefficients = native.prior.weights()[7]
    den = native.prior.FRAC**2 * formation.DEN
    err = iv.mpf(["-3e-18", "3e-18"])
    raw = json.loads((measure.HERE / "results.json").read_text())
    slopes = {}
    min_margin = 1
    from .readout import radius

    for charge in (-1, 1):
        sectors = {}
        for odd in (False, True):
            h = native.generator(board, ports, charge, odd)
            blocks = four_tick_columns(h, end, start, coefficients).reshape(
                2, 4, 4, 2, 4
            )
            for i, coin in enumerate(coins):
                z = []
                for j in range(4):
                    left, right = (blocks[:, j, :, s, :] for s in (0, 1))
                    dc = (
                        contraction(left, left, coin)[0]
                        + contraction(right, right, coin)[0]
                    )
                    re, im = contraction(right, left, coin)
                    z.append(
                        (
                            iv.mpf(dc) / (2 * den) + err,
                            iv.mpc(iv.mpf(re) / den + err, iv.mpf(im) / den + err),
                        )
                    )
                sectors[odd, i] = z
        prod, formed = sectors[False, 0], sectors[False, 1]
        slopes[str(charge)] = {}
        for endpoint, signs in (("product", (1, 0)), ("formed", (0, 1))):
            upper = (prod, formed)[signs[0]]
            lower = (prod, formed)[signs[1]]
            dp = [
                (1 if endpoint == "product" else -1) * (num[1] / deno[1]).imag
                for num, deno in zip(upper, lower, strict=True)
            ]
            slopes[str(charge)][endpoint] = [
                [float((v - dp[0]).a) - 1e-13, float((v - dp[0]).b) + 1e-13]
                for v in dp[1:]
            ]
        for stock in measure.STOCKS:
            base = prod if stock in ("product", "ordinary") else formed
            if stock == "half":
                panel = [
                    ((a + c) / 2, (b + d) / 2)
                    for (a, b), (c, d) in zip(prod, formed, strict=True)
                ]
            elif stock == "ordinary":
                panel = [
                    ((a + c) / 2, (b + d) / 2)
                    for (a, b), (c, d) in zip(prod, sectors[True, 0], strict=True)
                ]
            else:
                panel = [(a, 0 if stock == "scrambled" else b) for a, b in base]
            for p, counts in zip(
                measure.PHASES, raw["counts"][str(charge)][stock], strict=True
            ):
                exposure = iv.exp(iv.j * iv.pi * p / 4)
                probs = [dc + (f * exposure).real for dc, f in panel]
                probs.append(1 - sum(probs))
                seen = np.asarray(counts) / raw["shots"]
                for j, prob in enumerate(probs):
                    error = (
                        max(abs(float(prob.a) - seen[j]), abs(float(prob.b) - seen[j]))
                        + 1e-15
                    )
                    min_margin = min(min_margin, radius(seen[j], raw["shots"]) - error)
    assert min_margin > 1e-14
    return dict(
        phase_slopes_radians_per_pairing=slopes,
        minimum_probability_margin=float(min_margin),
        real_component_operator_and_source_error=3e-18,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = measure.HERE / "certificate.json"
    result = compute()
    if args.write:
        with path.open("x") as out:
            out.write(json.dumps(result, indent=2) + "\n")
    else:
        saved = json.loads(path.read_text())
        for q, panels in result["phase_slopes_radians_per_pairing"].items():
            for s, v in panels.items():
                np.testing.assert_allclose(
                    saved["phase_slopes_radians_per_pairing"][q][s],
                    v,
                    atol=1e-13,
                    rtol=0,
                )
        assert saved["minimum_probability_margin"] > 1e-14
        print(
            "Verified both species, four contacts, two pairing derivatives and all raw frequency intervals."
        )


if __name__ == "__main__":
    main()
