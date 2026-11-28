"""Source-to-exterior dispersive clock inference on the occupied finite chart."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json

import numpy as np

from . import clock_covariance as clock
from . import occupied_release as flight
from . import positive_resolvent as proof
from . import regional_tomography as meter

GAPS = (Fraction(1), Fraction(1, 4), Fraction(1, 16))
DIVISORS = (2, 3, 4, 5, 6)


def shell_spectra(twice, marks, rad):
    data = []
    for a, b in flight.SECTORS:
        e, v = np.linalg.eigh(twice / 2 + np.diag(a * marks[0] + b * marks[1]))
        # Uniformly chosen launch and return site, not a coherent shell mode.
        weights = np.array([(v[rad == r] ** 2).mean(0) for r in range(max(rad) + 1)])
        data.append((e, weights, v[(rad == 1) | (rad == 2)] ** 2))
    return data


def linewidth_records(data, omega, gap, law, delta, exact, rng):
    # Richardson coefficients for 1/(d^2+lambda^2), evaluated at zero;
    # the omitted kernel has the same positive sign on every spectral line.
    widths = [gap / k for k in DIVISORS]
    coeff = []
    for w in widths:
        v = Fraction(1)
        for other in widths:
            if other != w:
                v *= other**2 / (other**2 - w**2)
        coeff.append(v)
    tail = Fraction(1) / gap
    for w in widths:
        tail *= w * w / (gap * gap + w * w)
    measured, errors, records = [], [], []
    numerical = 5e-12
    for width in widths:
        eta = float(width)
        conditional = np.array(
            [p @ (eta * (e - omega) / ((e - omega) ** 2 + eta**2)) for e, p, _ in data]
        )
        probs = (1 + law @ conditional) / 2
        counts = rng.binomial(clock.N, probs)
        middle = 2 * (counts[1] - counts[0]) / clock.N
        margin = 4 * clock.EPS
        pred = flight.MIX @ flight.CONDITIONAL @ conditional
        forecast = np.sort(np.outer(pred, delta), axis=1) + np.array(
            [-numerical, numerical]
        )
        bands = middle[:, None] + np.array([-margin, margin])
        assert np.all(
            np.maximum(bands[:, 0], forecast[:, 0])
            <= np.minimum(bands[:, 1], forecast[:, 1])
        )
        records.append(
            dict(
                linewidth=str(width),
                counts=counts.tolist(),
                forecast=forecast.tolist(),
                signed_held=bands.tolist(),
            )
        )
        measured.append(middle / eta)
        errors.append(margin / eta)
    weights = np.array(coeff, float)
    estimate = weights @ measured
    margin = abs(weights) @ errors + delta[1] * float(tail) / 2
    bracket = estimate[:, None] + np.array([-margin, margin])
    assert np.all(
        np.maximum(bracket[:, 0], exact[:, 0]) <= np.minimum(bracket[:, 1], exact[:, 1])
    )
    return dict(
        receivers=records,
        extrapolation=[str(x) for x in coeff],
        one_spectrum_remainder=str(tail),
        held_narrow_line=bracket.tolist(),
    )


def finite_rail(spectra, rad, omega, gap, delta, work, exchange):
    # Scaled lowest-root equation: z = <x|(H-omega+g^2 z)^(-1)|x>.
    # Subtracting the two full eigenvalues near omega would lose precision.
    from_one = rad[(rad == 1) | (rad == 2)] == 1
    read = []
    for e, _, sites in spectra:
        lo = np.zeros(len(sites))
        hi = np.full(len(sites), 1 / float(gap))
        for _ in range(60):
            z = (lo + hi) / 2
            response = (
                sites / (e[None, :] - float(omega) + float(exchange) ** 2 * z[:, None])
            ).sum(1)
            lower = z < response
            lo = np.where(lower, z, lo)
            hi = np.where(lower, hi, z)
        scaled = (lo + hi) / 2
        read.append(scaled[from_one].mean() - scaled[~from_one].mean())
    mixed = flight.MIX @ flight.CONDITIONAL @ read
    g2 = float(exchange**2)
    # Each arm residual in z is between zero and g^2/gap^3; TV
    # per unit D is 1/2. Taking a difference of shells doubles the bound.
    padding = float(exchange**2 / gap**3)
    per_d = np.nextafter(
        work + np.array([-padding, padding]), np.array([-np.inf, np.inf])
    )
    assert 0 < per_d[0] <= mixed <= per_d[1]
    low, high = meter.product_interval(per_d, delta) * g2
    return dict(
        exchange=str(exchange),
        near_to_outer_work=[np.nextafter(low, -np.inf), np.nextafter(high, np.inf)],
        independent_scaled_root_per_D=mixed,
        primary_scaled_work_per_D=work.tolist(),
    )


def ordinary_change(before, after, marks, omega, gap, rad, p):
    # Native gain times the positive products of two resolvent columns.
    difference = (after - before).diagonal()
    sites = np.flatnonzero(difference)
    gain = Fraction(int(difference[sites[0]]), 2)
    np.testing.assert_array_equal(difference[sites], 2 * float(gain))
    sums = np.zeros((2, len(rad)), dtype=object)
    # Native head marginals/conditionals are dyadic multiples of 1/7168.
    nums = np.rint(7168 * p).astype(int)
    np.testing.assert_allclose(nums / 7168, p, rtol=0, atol=2e-15)
    for (a, b), weight in zip(flight.SECTORS, nums, strict=True):
        addition = np.diag(2 * (a * marks[0] + b * marks[1]))
        first = proof.columns(before + addition, omega, gap, sites)
        second = proof.columns(after + addition, omega, gap, sites)
        for side in (0, 1):
            sums[side] += int(weight) * (first[side] * second[side]).sum(1)
    assert sums.min() > 0
    answer = []
    for r in range(max(rad) + 1):
        take = rad == r
        raw = [
            float(-gain * Fraction(sum(v[take]), 7168 * proof.S**2 * int(take.sum())))
            for v in sums
        ]
        answer.append([np.nextafter(min(raw), -np.inf), np.nextafter(max(raw), np.inf)])
    return dict(gain=str(gain), signed_shell_response=answer)


def run():
    board, old, edges, rad, _ = flight.chart()
    paired = json.loads((flight.HERE / "clock_flight.json").read_bytes())
    clamp = paired["clock"]["ordinary_clamp"]
    r = rad[np.arange(board.n) != clamp]
    law = np.array(flight.sector_prior(old)) @ flight.CONDITIONAL
    delta = flight.lab.probe.log_factor(old["source"])
    rng = np.random.default_rng(6207031)
    result = {}
    for name, load in (("vacancy", 0), ("excess", 2)):
        before, marks = flight.apparatus(board, edges, load, clamp, 1)
        after, _ = flight.apparatus(board, edges, load, clamp, 2)
        edge = proof.floor_band(before)
        lower = min(edge, proof.floor_band(after))
        spectra = shell_spectra(before, marks, r)
        bath = {}
        cache = {}
        for category, band in (("inactive", edge), ("ordinary", lower)):
            selected = {}
            for gap in GAPS:
                omega = band - gap
                if omega not in cache:
                    mixed, bits = proof.mixed_columns(before, marks, omega, gap)
                    mixed_shell = proof.shell_intervals(mixed, r)
                    transfer = np.column_stack(
                        (mixed_shell[:, 0] * delta[0], mixed_shell[:, 1] * delta[1])
                    )
                    transfer = np.nextafter(transfer, np.array([-np.inf, np.inf]))
                    # Work per D at each adjacent pair, in inverse energy units.
                    work = np.array(
                        [
                            meter.difference_interval(a, b)
                            for a, b in zip(
                                mixed_shell[:-1], mixed_shell[1:], strict=True
                            )
                        ]
                    )
                    work = np.nextafter(work, np.array([-np.inf, np.inf]))
                    # At a detuned, weakly coupled additional rail, this
                    # bound on g^2 certifies the sign of the *finite* level difference.
                    limits = np.nextafter(np.maximum(0, work[:, 0]) * float(gap**3), 0)
                    records = linewidth_records(
                        spectra, float(omega), gap, law, delta, transfer, rng
                    )
                    cache[omega] = dict(
                        frequency=str(omega),
                        gap=str(gap),
                        residual_ticks=bits,
                        mixed_shell=mixed_shell.tolist(),
                        source_prediction=transfer.tolist(),
                        adjacent_work_per_source=work.tolist(),
                        positive_work_exchange_squared_below=limits.tolist(),
                        **records,
                    )
                row = dict(cache[omega])
                if category == "ordinary":
                    row["added_clamp"] = ordinary_change(
                        before, after, marks, omega, gap, r, law[0]
                    )
                selected[str(gap)] = row
            bath[category] = selected
        result[name] = bath
    # Choose one dyadic exchange certified for both probes and all primary
    # detunings; this is an instrument setting selected by the local proof.
    threshold = min(
        v["positive_work_exchange_squared_below"][1]
        for v0 in result.values()
        for v in v0["inactive"].values()
    )
    exchange = Fraction(1, 2)
    while float(exchange**2) >= threshold / 2:
        exchange /= 2
    for name, load in (("vacancy", 0), ("excess", 2)):
        before, marks = flight.apparatus(board, edges, load, clamp, 1)
        spectra = shell_spectra(before, marks, r)
        for row in result[name]["inactive"].values():
            row["finite_rail"] = finite_rail(
                spectra,
                r,
                Fraction(row["frequency"]),
                Fraction(row["gap"]),
                delta,
                np.array(row["adjacent_work_per_source"][1]),
                exchange,
            )
    return dict(
        sha256={
            f: sha256((flight.HERE / f).read_bytes()).hexdigest()
            for f in ("DISPERSIVE_CLOCK.md", "pair_contact.json", "clock_flight.json")
        },
        clamp=clamp,
        source_displacement=delta.tolist(),
        fraction_bits=proof.BITS,
        shell_counts=np.bincount(r).tolist(),
        trials=clock.N,
        results=result,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = run()
    path = flight.HERE / "dispersive_clock.json"
    if args.write:
        path.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    else:
        flight.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), output
        )
    for species, readings in output["results"].items():
        for gap, obs in readings["inactive"].items():
            print(
                species, gap, obs["source_prediction"][:3], obs["held_narrow_line"][0]
            )
