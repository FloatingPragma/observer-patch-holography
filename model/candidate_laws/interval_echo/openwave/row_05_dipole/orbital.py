"""Compare one vacancy's resolved angular recoil and odd magnetic splitting."""

import argparse
import json
import itertools
import numpy as np
from scipy.linalg import expm
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls import belt
from model.candidate_laws.interval_echo.openwave.row_29_maxwell.experiment import (
    geometry,
)
from model.executable_model import cell
from .experiment import HERE, CHARGE, portable

SHOTS = 2**25
FIELDS = (-1 / 32, -1 / 64, 0, 1 / 64, 1 / 32)
# 3 laws x at most 8 two-quadrature readings, 4 tour reads, 5 clock reads,
# plus the calibration controls. One bound covers all Bernoulli sample means.
RADIUS = np.sqrt(2 * np.log(2 * 400 / 0.01) / SHOTS)


def fixture():
    order = list(next(iter(belt.belts().values())))
    xyz = geometry()
    poles = sorted(set(range(12)) - set(order))
    axis = xyz[poles[0]].copy()
    if np.cross(xyz[order[0]], xyz[order[1]]) @ axis < 0:
        axis *= -1
    area = np.cross(xyz[:, None], xyz[None, :]) @ axis / 2
    r2 = 1 - (xyz[order[0]] @ axis) ** 2
    return order, area, r2


def generator(field, variant, twist=0):
    board = cell()
    order, area, _ = fixture()
    law = candidate.native_law
    strength = 0 if variant == "exchange" else candidate.KAPPA
    local = law.seam_term(strength)
    i, j = (law.SRC.index(w) for w in ((0, 1), (1, 0)))
    # The uniform host has zero seam cost on every equal occupied bond.
    h = np.diag(
        np.asarray(board.degree, float) * local[i, i] - (6 - board.degree) * strength
    ).astype(complex)
    k = np.diag(board.degree / 2).astype(complex)
    for a, b in board.seams:
        phase = CHARGE * field * area[a, b]
        if twist and a in order and b in order:
            step = (order.index(b) - order.index(a) + 5) % 10 - 5
            phase += 2 * np.pi * twist * step / 10
        hop = local[j, i] * np.exp(1j * phase)
        h[b, a] = k[b, a] = hop
        h[a, b] = k[a, b] = hop.conjugate()
    if variant == "two_hop":
        h += k @ k / 4
    return h[np.ix_(order, order)]


def read(z, rng):
    records = []
    for v in (z.real, z.imag):
        prob = (1 + v) / 2
        if not -1e-12 <= prob <= 1 + 1e-12:
            raise ArithmeticError("invalid fringe")
        n = int(rng.binomial(SHOTS, round(min(1, max(0, prob)), 12)))
        records.append([n, SHOTS - n])
    return records


def phase(counts):
    x, y = [(row[0] - row[1]) / sum(row) for row in counts]
    norm = np.hypot(x, y)
    return float(np.arctan2(y, x)), float(np.arcsin(min(1, np.sqrt(2) * RADIUS / norm)))


def modes():
    return np.exp(2j * np.pi * np.outer(np.arange(10), np.arange(10)) / 10) / np.sqrt(
        10
    )


def pair_counts(u, a, b, rng):
    v = u @ ((a + b) / np.sqrt(2))
    left, right = np.vdot(a, v), np.vdot(b, v)
    weight = abs(left) ** 2 + abs(right) ** 2
    score = 2 * left.conjugate() * right
    if abs(weight - 1) < 1e-12:
        return [[p, m, 0] for p, m in read(score, rng)]
    probabilities = [
        [(weight + v) / 2, (weight - v) / 2, 1 - weight]
        for v in (score.real, score.imag)
    ]
    return [rng.multinomial(SHOTS, row).tolist() for row in probabilities]


def calibration(field, rng, twist=0):
    order, area, r2 = fixture()
    law = candidate.native_law
    idx = [law.SRC.index(w) for w in ((0, 1), (1, 0))]
    base = law.seam_term(candidate.KAPPA)[np.ix_(idx, idx)]

    def tour(path):
        amp = 1 + 0j
        for a, b in zip(path, path[1:] + path[:1], strict=True):
            s = base.astype(complex)
            step = (order.index(b) - order.index(a) + 5) % 10 - 5
            s[1, 0] *= np.exp(
                1j * (CHARGE * field * area[a, b] + 2 * np.pi * twist * step / 10)
            )
            s[0, 1] = s[1, 0].conjugate()
            amp *= expm(-1j * np.pi * s)[1, 0]
        return amp

    count = read(tour(order) * tour(order[::-1]).conjugate(), rng)
    angle, err = phase(count)
    twice_area = sum(
        2 * area[a, b] for a, b in zip(order, order[1:] + order[:1], strict=True)
    )
    return dict(
        counts=count,
        f=angle * r2 / twice_area,
        error=err * r2 / twice_area,
        winding=angle / (4 * np.pi),
        winding_error=err / (4 * np.pi),
    )


def run():
    rng = np.random.default_rng(510247)
    order, _, r2 = fixture()
    f = modes()
    clocks = json.loads((HERE / "calibration.json").read_text())["clock_counts_z_y"]

    def clock(t):
        # The clock publishes consecutive short arcs; unwrap the local record
        # before taking the difference. Endpoint errors bound the interval.
        points = [phase(pair) for pair in clocks[: t + 1]]
        angles = np.unwrap([ang for ang, _ in points])
        # These existing clock counts used a different number of shots.
        n = sum(clocks[0][0])
        radius = np.sqrt(2 * np.log(2 * 400 / 0.01) / n)
        for pair in clocks[: t + 1]:
            v = np.array([(a - b) / n for a, b in pair])
            if np.linalg.norm(v) < 0.98:
                raise ArithmeticError("weak clock fringe")
        return abs(angles[-1] - angles[0]), 2 * np.arcsin(np.sqrt(2) * radius / 0.98)

    tours = {b: calibration(b, rng) for b in FIELDS}
    t0, dt0 = clock(1)
    t1, dt1 = clock(32)
    cases = []
    for variant in ("registered", "exchange", "two_hop"):
        u = expm(-1j * generator(0, variant))
        gaps = []
        # The short exposure resolves the long reading's phase branch.
        for m in (1, 2):
            short = pair_counts(
                expm(-0.25j * generator(0, variant)), f[:, m], f[:, 0], rng
            )
            long = pair_counts(u, f[:, m], f[:, 0], rng)
            a, _ = phase(short)
            z, e = phase(long)
            z += 2 * np.pi * round((4 * a - z) / (2 * np.pi))
            gaps.append(
                dict(
                    m=m,
                    short=short,
                    counts=long,
                    gap=z / t0,
                    lower=(z - e) / (t0 + dt0),
                    upper=(z + e) / (t0 - dt0),
                )
            )
        pulls = []
        for b in FIELDS:
            U = expm(-32j * generator(b, variant))
            c = pair_counts(U, f[:, -1], f[:, 1], rng)
            angle, er = phase(c)
            rec = dict(
                field=b,
                counts=c,
                splitting=angle / t1,
                zero_winding_return=float(abs(np.vdot(f[:, 0], U @ f[:, 0])) ** 2),
            )
            if b:
                norm = tours[b]
                pulls_f = [
                    abs(norm["f"]) - norm["error"],
                    abs(norm["f"]) + norm["error"],
                ]
                g = abs(angle) / (2 * t1 * abs(norm["f"]) * gaps[0]["gap"])
                values = [
                    v / (2 * t * bhat * gap)
                    for v, t, bhat, gap in itertools.product(
                        [max(0, abs(angle) - er), abs(angle) + er],
                        [t1 - dt1, t1 + dt1],
                        pulls_f,
                        [gaps[0]["lower"], gaps[0]["upper"]],
                    )
                ]
                rec.update(g=g, lower=min(values), upper=max(values))
            pulls.append(rec)
        cases.append(dict(variant=variant, gaps=gaps, flux=pulls))
    # These draws follow the original recoil/Zeeman runs in the same stream.
    curvature = []
    for twist in (-1 / 8, -1 / 16, 1 / 16, 1 / 8):
        cal = calibration(0, rng, twist)
        parked = expm(-32j * generator(0, "registered")) @ f[:, 0]
        thread = expm(-32j * generator(0, "registered", twist)) @ f[:, 0]
        counts = read(np.vdot(thread, parked), rng)
        a, err = phase(counts)
        nu, du = abs(cal["winding"]), cal["winding_error"]
        c = 2 * a / t1 / nu**2
        lo = 2 * (a - err) / (t1 + dt1) / (nu + du) ** 2
        hi = 2 * (a + err) / (t1 - dt1) / (nu - du) ** 2
        # The finite difference is a lower bound on the limiting curvature
        # for the registered two harmonic bands. The fourth-order alternating
        # cosine bound holds for the complete measured twist interval.
        delta = 2 * np.pi / 10
        bound = 1 - (2 * delta * (nu + du)) ** 2 / 12
        hi = hi / bound
        normalized = []
        for row in cases[0]["flux"]:
            b = row["field"]
            if not b:
                continue
            z, dz = phase(row["counts"])
            B, dB = abs(tours[b]["f"]), tours[b]["error"]
            normalized.append(
                dict(
                    field=b,
                    coefficient=abs(z) / (t1 * B * c),
                    lower=max(0, abs(z) - dz) / ((t1 + dt1) * (B + dB) * hi),
                    upper=(abs(z) + dz) / ((t1 - dt1) * (B - dB) * lo),
                )
            )
        curvature.append(
            dict(
                applied=twist,
                calibration=cal,
                counts=counts,
                curvature=c,
                limit_lower=lo,
                limit_upper=hi,
                normalized=normalized,
            )
        )
    sham = {}
    for b in (1 / 64, 1 / 32):
        left = expm(-32j * generator(-b, "registered")) @ f[:, 0]
        right = expm(-32j * generator(b, "registered")) @ f[:, 0]
        sham[str(b)] = read(np.vdot(left, right), rng)
    return dict(
        zero_winding_reversal=sham,
        curvature=curvature,
        seed=510247,
        shots=SHOTS,
        radius=float(RADIUS),
        order=order,
        radius_squared=r2,
        clock=[t0, t1],
        clock_errors=[dt0, dt1],
        tours={str(k): v for k, v in tours.items()},
        calibration_quarter_turn=read(1j, rng),
        cases=cases,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    out = portable(run())
    path = HERE / "orbital.json"
    if p.parse_args().write:
        path.write_text(json.dumps(out, indent=1) + "\n")
    else:
        if out != json.loads(path.read_text()):
            raise ValueError("Angular recoil or flux records differ")
        print("Verified single-vacancy angular recoil and signed flux splittings.")


if __name__ == "__main__":
    main()
