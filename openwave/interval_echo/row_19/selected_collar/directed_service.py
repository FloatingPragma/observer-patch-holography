"""Measure the transport of the parity-energy tag on the interrupted service flight."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy.sparse.linalg import expm_multiply

from . import distributed_source as dist

HERE, reg, iv, lab = dist.HERE, dist.reg, dist.iv, dist.lab
region = dist.regions
STEP = 1 / 256
SEED = 28112122
STENCIL = (-1, 9, -45, 0, 45, -9, 1)


def microsteps(generator, columns):
    return expm_multiply(
        generator,
        columns,
        start=0,
        stop=6 * STEP,
        num=7,
        traceA=generator.diagonal().sum(),
    ).reshape(7, -1, 4, 4)


def bond_set(H, ports):
    return sorted(
        {
            tuple(sorted((i, int(j))))
            for i in ports
            for j in H[4 * i : 4 * i + 4].nonzero()[1] // 4
            if i != j
        }
    )


def acquire():
    rng = np.random.default_rng(SEED)
    native = lab.six.survey.native
    board = native.iso.pairing.formation.federation(1)
    route = native.iso.pairing.formation.region(board)
    ports = list(dist.read("distributed_source")["ports"])
    knobs = {
        int(t): complex(*v)
        for t, v in dist.read("six_meter_pointer")["central_settings"].items()
    }
    _, sources, *_ = region.c.trials()
    histories = {}
    scales = {}
    H_bounds = {}
    for sign in (-1, 1):
        # The score is exactly zero in the even parity sector.
        H = native.generator(board, route, sign, True).tocsr() / 4
        H_bounds[str(sign)] = float(abs(H).sum(axis=1).max())
        generator = -1j * H
        bonds = bond_set(H, ports)
        if sign == -1:
            enrolled = bonds
        else:
            assert enrolled == bonds
        blocks = np.array(
            [H[4 * i : 4 * i + 4, 4 * j : 4 * j + 4].toarray() for i, j in bonds]
        )
        singular = np.linalg.svd(blocks, compute_uv=False)
        vmax = singular[:, 0]
        np.testing.assert_allclose(
            singular**2 * (singular**2 - vmax[:, None] ** 2), 0, atol=1e-12
        )
        # The native blocks have exact integer quarter entries and a single
        # nonzero singular value. Verify the squared projector identity exactly.
        for A, v in zip(blocks, vmax, strict=True):
            M = np.rint(4 * A).astype(np.int64)
            np.testing.assert_array_equal(4 * A, M)
            square = M @ M.T
            np.testing.assert_array_equal(square @ square, round((4 * v) ** 2) * square)
        vmax = np.round(vmax * 4) / 4
        scales[str(sign)] = vmax.tolist()
        left = np.array([i for i, j in bonds])
        right = np.array([j for i, j in bonds])

        def bond_means(flights, stock):
            gamma, click, miss = flights

            def local(amplitude):
                a, b = amplitude[left], amplitude[right]
                Ha = blocks.conj().transpose(0, 2, 1) @ a
                Hb = blocks @ b

                def inner(u, v):
                    return np.einsum("eij,jk,eik->e", u, stock, v.conj(), optimize=True)

                return inner(Hb, a), (inner(Ha, Ha) + inner(Hb, Hb)).real / vmax**2

            co, hit = local(miss + gamma * click)
            other, ohit = local(click)
            return co + (1 - abs(gamma) ** 2) * other, hit + (
                1 - abs(gamma) ** 2
            ) * ohit

        for home in ports:
            init = np.zeros((H.shape[0], 4))
            init[4 * home : 4 * home + 4] = np.eye(4)
            Uhome = microsteps(generator, init)
            tube = region.density.flight_columns(generator, init, max(knobs))
            for tau, gamma in knobs.items():
                free = microsteps(generator, tube[tau])
                click = Uhome @ tube[tau, 4 * home : 4 * home + 4]
                miss = free - click
                frame = {}
                for arm, g, stock in (
                    ("dense", 0, sources["reference"][1]),
                    ("ordinary", 0, sources["ordinary"][1]),
                    ("overlap", gamma, sources["reference"][1]),
                ):
                    co, hit = bond_means((g, click[3], miss[3]), stock)
                    readings = {}
                    for meter, score in (
                        ("current", 2 * co.imag),
                        ("work", 2 * co.real),
                    ):
                        readings[meter] = [
                            lab.counts(
                                np.array([(q + v / a) / 2, (q - v / a) / 2, 1 - q]), rng
                            )
                            for q, v, a in zip(hit, score, vmax, strict=True)
                        ]
                    selected = miss[:, ports] + g * click[:, ports]

                    def position(a):
                        return np.einsum(
                            "tsij,jk,tsik->ts", a, stock, a.conj(), optimize=True
                        ).real

                    probability = position(selected) + (1 - abs(g) ** 2) * position(
                        click[:, ports]
                    )
                    readings["positions"] = [
                        lab.counts(np.r_[row, 1 - row.sum()], rng)
                        for row in probability
                    ]
                    frame[arm] = readings
                histories[f"{tau}/{home}/{sign}"] = frame
    return {
        "seed": SEED,
        "starts": lab.SHOTS,
        "step": STEP,
        "ports": ports,
        "bonds": [list(e) for e in enrolled],
        "bond_scales": scales,
        "generator_row_norms": H_bounds,
        "histories": histories,
    }


def decode(z):
    iv.dps = 35
    dose = region.density.DOSE
    n = z["starts"]
    alpha = np.log(2 * reg.src.FAMILY / 0.01) / n

    def signed(tape, unit):
        plus, minus, zero = map(int, tape)
        assert min(tape) >= 0 and plus + minus + zero == n
        p = (plus + minus) / n
        width = unit * (
            np.sqrt(2 * alpha * (p + reg.src.probability_error(p, n)))
            + 4 * alpha / 3
            + 2e-13
        )
        mid = unit * ((plus - minus) / n)
        return reg.band(mid - width, mid + width)

    def positions(tape):
        assert len(tape) == 7 and sum(tape) == n and min(tape) >= 0
        p = np.asarray(tape[:-1]) / n
        e = reg.src.probability_error(p, n)
        return [reg.band(a - b, a + b) for a, b in zip(p, e, strict=True)]

    ports = z["ports"]
    edges = [tuple(v) for v in z["bonds"]]
    totals = {}
    archive = {}
    for h, arms in z["histories"].items():
        sign = h.rsplit("/", 1)[-1]
        scale = z["bond_scales"][sign]
        contrasts = {}
        for arm in ("ordinary", "overlap"):
            columns = {}
            rates = {}
            for meter in ("current", "work"):
                bond = {
                    e: dose * (signed(a, unit) - signed(b, unit))
                    for e, a, b, unit in zip(
                        edges,
                        arms[arm][meter],
                        arms["dense"][meter],
                        scale,
                        strict=True,
                    )
                }
                for i, j in edges:
                    if i in ports and j in ports:
                        columns[f"{meter}/{i}/{j}"] = bond[i, j]
                for site in ports:
                    outer = iv.mpf(0)
                    all_bonds = iv.mpf(0)
                    for e, v in bond.items():
                        if site not in e:
                            continue
                        oriented = -v if meter == "current" and site == e[1] else v
                        all_bonds += oriented
                        if not set(e) <= set(ports):
                            outer += oriented
                    columns[f"{meter}/{site}/outside"] = outer
                    if meter == "current":
                        rates[site] = all_bonds
            delta = []
            for high, low in zip(
                arms[arm]["positions"], arms["dense"]["positions"], strict=True
            ):
                delta.append(
                    [
                        dose * (a - b)
                        for a, b in zip(positions(high), positions(low), strict=True)
                    ]
                )
            # Seventh derivative of a probability on each separate arm is bounded
            # by (2||H||)^7; their difference has twice this bound.
            truncation = (
                2
                * dose
                * (2 * z["generator_row_norms"][sign]) ** 7
                * z["step"] ** 6
                / 140
            )
            slot = totals.setdefault(
                arm, {site: [iv.mpf(0) for _ in range(3)] for site in ports}
            )
            for i, site in enumerate(ports):
                snapshot = sum(
                    c * row[i] for c, row in zip(STENCIL, delta, strict=True)
                ) / (60 * z["step"])
                slot[site][0] += rates[site]
                slot[site][1] += snapshot + reg.band(-truncation, truncation)
                slot[site][2] += delta[0][i]
            contrasts[arm] = {k: lab.six.bounds(v) for k, v in columns.items()}
        archive[h] = contrasts
    continuity = {
        a: {
            str(site): dict(
                zip(
                    ("from_directed_bonds", "from_position_slope", "at_switch"),
                    [lab.six.bounds(v / len(archive)) for v in values],
                    strict=True,
                )
            )
            for site, values in panel.items()
        }
        for a, panel in totals.items()
    }
    return {"columns": archive, "continuity": continuity}


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = acquire()
    r.update(decode(r))
    r["inputs"] = {
        n: sha256((HERE / n).read_bytes()).hexdigest()
        for n in (
            "DIRECTED_SERVICE_PROTOCOL.md",
            "six_meter_pointer.json",
            "distributed_source.json",
        )
    }
    path = HERE / "directed_service.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    print(r["continuity"]["overlap"])
