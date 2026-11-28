"""Switch a real intermediate position interrogation while holding flight and energy fixed."""

import argparse
import json
from hashlib import sha256

import numpy as np
from scipy.sparse.linalg import expm_multiply

from . import position_source as pos

c, lab, iv, reg = pos.c, pos.lab, pos.iv, pos.reg
HERE = pos.HERE
WINDOWS = (1, 2, 4)
DOSE = 1 / 256
SEED = 28112103


def flight_columns(h, start, end):
    return expm_multiply(
        h, start, start=0, stop=end, num=end + 1, traceA=h.diagonal().sum()
    )


def history(
    input_source,
    ports,
    rng,
    pointers=None,
    *,
    signs=(-1, 1),
    read_ports=None,
    windows=WINDOWS,
    propagate=flight_columns,
):
    pointers = {} if pointers is None else pointers
    receivers = tuple(ports if read_ports is None else read_ports)
    native = lab.six.survey.native
    board = native.iso.pairing.formation.federation(1)
    region = native.iso.pairing.formation.region(board)
    local = {
        t: {
            k: np.zeros((len(receivers) + 1, *input_source.shape), dtype=complex)
            for k in ("sparse", "dense", "erased", *pointers)
        }
        for t in windows
    }
    yes = []
    for b in (0, 1):
        for sign in signs:
            h = native.generator(board, region, sign, bool(b)) * (-0.25j)
            for home in ports:
                start = np.zeros((h.shape[0], 4))
                start[4 * home : 4 * home + 4] = np.eye(4)
                end = 2 * max(windows)
                tube = propagate(h, start, end).reshape(end + 1, -1, 4, 4)

                def reduced(amplitude):
                    branches = amplitude.reshape(-1, board.n, 4, 4)
                    arrivals = np.einsum(
                        "hnij,jk,hnlk->nil",
                        branches,
                        input_source[b],
                        branches.conj(),
                        optimize=True,
                    )
                    enrolled = arrivals[list(receivers)]
                    return np.concatenate(
                        [enrolled, (arrivals.sum(0) - enrolled.sum(0))[None]]
                    ) / (len(signs) * len(ports))

                for tau in windows:
                    free = tube[2 * tau]
                    click = tube[tau] @ tube[tau, home]
                    miss = free - click
                    yes_prob = (
                        np.trace(
                            tube[tau, home] @ input_source[b] @ tube[tau, home].conj().T
                        ).real
                        / np.trace(input_source[b]).real
                    )
                    yes.append(
                        [b, sign, int(home), tau, int(rng.binomial(2**30, yes_prob))]
                    )
                    local[tau]["sparse"][:, b] += reduced(free)
                    local[tau]["dense"][:, b] += reduced(click) + reduced(miss)
                    # Trace the pointer explicitly with Kraus labels kept orthogonal.
                    local[tau]["erased"][:, b] += reduced(np.concatenate([click, miss]))
                    for name, gamma in pointers.items():
                        if abs(gamma) > 1:
                            raise ValueError(
                                "Pointer overlap must lie in the unit disk"
                            )
                        local[tau][name][:, b] += reduced(miss + gamma * click) + (
                            1 - abs(gamma) ** 2
                        ) * reduced(click)
    menus = {t: {k: z.sum(0) for k, z in m.items()} for t, m in local.items()}
    for menu in menus.values():
        for k, z in menu.items():
            np.testing.assert_allclose(
                np.trace(z, axis1=1, axis2=2),
                np.trace(input_source, axis1=1, axis2=2),
                rtol=0,
                atol=1e-12,
            )
            menu[k] = (1 - DOSE) * input_source + DOSE * z
        np.testing.assert_allclose(menu["erased"], menu["dense"], rtol=0, atol=1e-13)
    return menus, yes, local


def price_contrast(reference, high, low, which):
    mat, err = reference
    y, ey = high
    x, ex = low
    left, right = c.dose.heads.measure.marginals(mat)
    product = c.dose.heads.product_reference(mat)
    K = -c.station.log_blocks(mat)
    marg = -c.station.log_blocks(product)
    price = {"connected": K - marg, "full": K, "marginals": marg}[which]
    center = c.station.pair(price, y - x)
    width = 3e-12
    for b in (0, 1):
        referr = 2 * err[b]
        params = []
        if which in ("connected", "full"):
            params.append((np.linalg.eigvalsh(mat[b]).min(), referr))
        if which in ("connected", "marginals"):
            params.extend(
                [
                    (np.linalg.eigvalsh(left[b]).min(), 2 * referr),
                    (np.linalg.eigvalsh(right[b]).min(), 2 * referr),
                    (mat[b].trace().real, 2 * referr),
                ]
            )
        logerr = sum(e / (a - e) for a, e in params)
        coeff = abs(np.einsum("ijkl,lk->ij", c.dose.heads.PRODUCTS, price[b])).sum() / 4
        width += (
            coeff * (ex[b] + ey[b])
            + (np.linalg.norm(y[b] - x[b], ord="nuc") + 4 * (ex[b] + ey[b])) * logerr
        )
    return reg.band(center - width, center + width)


def run():
    iv.dps = 35
    rng = np.random.default_rng(SEED)
    _, inputs, *_ = c.trials()
    ref = inputs["reference"]
    ports, projectors = pos.calibrate()
    reader = pos.PositionReader(projectors)
    menus, yes, _ = history(ref, ports, rng)
    states = {"reference": ref}
    for tau, menu in menus.items():
        states.update({str(tau) + k: v for k, v in menu.items()})
    _, bins = lab.six.survey.return_bins()
    ids = [i for i, v in enumerate(lab.six.survey.ROUTES) if max(v) < 4]
    bins = bins[ids]
    raw_head = {}
    raw_radar = {}
    raw_positions = {}
    n = 2**62
    for k, x in states.items():
        raw_head[k] = lab.pairs(x, rng, False)
        raw_radar[k] = [
            lab.counts(np.append(p, 1 - p.sum()), rng) for p in c.station.pair(bins, x)
        ]
        raw_positions[k] = pos.emit(c.station.pair(projectors, x), n, rng)
    heads = {k: reg.src.observed(t, lab.SHOTS, False) for k, t in raw_head.items()}
    positions = {k: reader.decode(t, n) for k, t in raw_positions.items()}
    times, _ = reg.timing_slopes(raw_radar, lab.SHOTS, (("1sparse", "1dense", 1),))
    g = iv.mpf([23.1042361, 23.1042363])
    source = HERE / "link_source_contraction.json"
    multiplier = iv.mpf(json.loads(source.read_bytes())["ordinary_volume_multiplier"])
    findings = {}
    for tau in WINDOWS:
        dst = {}
        for hi, lo in (("sparse", "dense"), ("erased", "dense")):
            a, b = (str(tau) + k for k in (hi, lo))
            prices = {
                k: price_contrast(heads["reference"], heads[a], heads[b], k)
                for k in ("connected", "full", "marginals")
            }
            src_a, src_b = (c.dose.source_read(raw_head[k]) for k in (a, b))
            ra, rb = (np.array(times[k], dtype=object) for k in (a, b))
            dr = ra - rb
            box = np.array(
                [
                    iv.mpf([min(u.a, v.a), max(u.b, v.b)])
                    for u, v in zip(ra, rb, strict=True)
                ]
            )
            A, V = pos.geometry_jet(box)
            # One common direction in raw times for both components.
            corrected = (V.g * multiplier - A.g * g) @ dr
            px = positions["reference"]
            Kref = c.station.log_blocks(
                c.dose.heads.product_reference(px["matrix"])
            ) - c.station.log_blocks(px["matrix"])
            aa, bb = positions[a], positions[b]
            # Q(a-b) at independent common reference, with a common logarithm uncertainty.
            rad = px["radius"]
            left, right = c.dose.heads.measure.marginals(px["matrix"])
            params = [
                (np.linalg.eigvalsh(px["matrix"]).min(), rad),
                (np.linalg.eigvalsh(left).min(), np.sqrt(2) * rad),
                (np.linalg.eigvalsh(right).min(), np.sqrt(2) * rad),
                (np.trace(px["matrix"], axis1=1, axis2=2).real.min(), 2 * rad),
            ]
            assert all(f > e + 2e-12 for f, e in params)
            log_error = sum(e / (f - e - 2e-12) for f, e in params) + 1e-11
            err = (
                np.sqrt(8)
                * (
                    np.linalg.norm(aa["matrix"] - bb["matrix"])
                    + aa["radius"]
                    + bb["radius"]
                    + 1e-12
                )
                * log_error
            )
            position_q = reader.change(bb, aa, Kref) + reg.band(-err, err)
            sa, sb = pos.entropy(reader, aa), pos.entropy(reader, bb)
            dst[hi + "_minus_" + lo] = {
                **{k: lab.six.bounds(v) for k, v in prices.items()},
                "head_entropy_change": lab.six.bounds(src_a - src_b),
                "position_charge": lab.six.bounds(position_q),
                "position_entropy_change": lab.six.bounds(sa - sb),
                "area": lab.six.bounds(A.g @ dr),
                "volume": lab.six.bounds(V.g @ dr),
                "volume_adjusted_contraction_nats": lab.six.bounds(corrected),
                "adjusted_minus_connected": lab.six.bounds(
                    corrected - prices["connected"]
                ),
                "adjusted_minus_full": lab.six.bounds(corrected - prices["full"]),
            }
        findings[str(tau)] = dst
    return {
        "seed": SEED,
        "interrogation_probability_in_dense_arm": 1,
        "flight_mixture_fraction": DOSE,
        "native_position_starts": n,
        "head_and_radar_starts": lab.SHOTS,
        "intermediate_yes_of_2pow30": yes,
        "position_counts": raw_positions,
        "head_counts": raw_head,
        "radar_counts": raw_radar,
        "windows": findings,
        "inputs": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "OBSERVATION_DENSITY_PROTOCOL.md",
                HERE / "finite_pairing_dose.json",
                source,
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    data = run()
    path = HERE / "observation_density.json"
    if cli.parse_args().write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), data)
    for k, d in data["windows"].items():
        print(k, d["sparse_minus_dense"])
