"""Propagate native collar comparisons to ordinary-calibrated, signed material recoil."""

import argparse
import json
from pathlib import Path
from itertools import product

import numpy as np
from scipy.fft import irfftn
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.experiments.interval_echo_einstein_station.mobile_collar_source import (
    measure as formation,
)
from ... import measure as green
from ..common_receiver import measure as receiver
from ..energy_source.binding.measure import SCALE

HERE = Path(__file__).parent
INTERIORS = (10, 34, 82, 178)


def central_load(n):
    load = np.zeros((n,) * 3)
    mid = n // 2
    load[mid - 1 : mid + 1, mid - 1 : mid + 1, mid - 1 : mid + 1] = 1 / 8
    return load


def diagonal_visits(n, q):
    """Odd reflection about 0 and n+1, with zero mode cancelling in the image sum."""
    m = 2 * (n + 1)
    lap = 2 * q * (1 - np.cos(2 * np.pi * np.arange(m) / m))
    eigen = lap[:, None, None] + lap[None, :, None] + lap[None, None, : m // 2 + 1]
    eigen[0, 0, 0] = np.inf
    np.reciprocal(eigen, out=eigen)
    p = irfftn(eigen, s=(m, m, m), axes=(0, 1, 2))
    i, j, k = np.ix_(*([2 * np.arange(1, n + 1)] * 3))
    return (
        p[0, 0, 0]
        - p[i, 0, 0]
        - p[0, j, 0]
        - p[0, 0, k]
        + p[i, j, 0]
        + p[i, 0, k]
        + p[0, j, k]
        - p[i, j, k]
    )


def local_source(z):
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    ids = [law.SRC.index((i, j)) for i in (0, 1) for j in (0, 1)]
    embedding = np.eye(law.DIM)[:, ids]
    h = exchange + candidate.KAPPA * contrast
    block = embedding.T @ h @ embedding
    np.testing.assert_allclose(h @ embedding, embedding @ block, atol=3e-14)
    levels, basis = np.linalg.eigh(block)
    np.testing.assert_allclose(levels, [0, 0, 8, 9], atol=1e-13)
    park = formation.local_parts(1, 8)[0, 0]
    np.testing.assert_allclose(
        np.kron(park, np.eye(2)) + np.kron(np.eye(2), park), 16 * np.eye(4)
    )
    calibration = block[1, 1] - block[0, 0]
    rows = []
    for sigma in (1, -1):
        for erase in (0.0, 0.5, 1.0):
            for switches in ("shared", "independent"):
                retained = (1 - erase) ** (1 if switches == "shared" else 2)
                signed = sigma * z * z * retained
                rho = np.diag((1 + signed * np.array([1, -1, -1, 1])) / 4)
                np.testing.assert_allclose(
                    rho.diagonal().reshape(2, 2).sum(axis=0), 0.5
                )
                np.testing.assert_allclose(
                    rho.diagonal().reshape(2, 2).sum(axis=1), 0.5
                )
                for t in (0, 1, 8, 64):
                    u = expm(-1j * t * block)
                    np.testing.assert_allclose(u @ rho @ u.conj().T, rho, atol=8e-14)
                spectral = (basis.T @ rho @ basis).diagonal()
                mismatch = (1 - signed) / 2
                np.testing.assert_allclose(
                    [sum(spectral[:2]), *spectral[2:]],
                    [1 - mismatch, mismatch / 2, mismatch / 2],
                )
                # Select contrast tickets with probability 8/512 and odd parity
                # with 1/512; the remaining tickets never fire on binary loads.
                fired = (8 * mismatch + spectral[-1]) / SCALE
                energy = np.trace(block @ rho)
                np.testing.assert_allclose(fired, energy / SCALE, rtol=2e-14)
                rows.append(
                    dict(
                        orientation=sigma,
                        erasure=erase,
                        switches=switches,
                        energy=float(energy),
                        pulse_probability=float(fired),
                    )
                )
    return float(calibration), rows


def conservative_preparation(z):
    """Post-design control: move the energy to the reader from a second ordinary seam."""
    board, rails, flag = green.rails.layout()
    occupied = set().union(*rails, flag)
    # Choose a four-cycle with two storage edges and two perpendicular transfer edges.
    cycle = next(
        (int(a), int(b), int(c), int(d))
        for a in range(board.n)
        if a not in occupied
        for b in board.neighbours[a]
        if b not in occupied
        for c in board.neighbours[a]
        if c != b and c not in occupied
        for d in board.neighbours[b]
        if d in board.neighbours[c] and d != a and d not in occupied
    )
    law = candidate.native_law
    ex, contrast = law.local_term_parts()
    ids = [law.SRC.index(word) for word in product((0, 1), repeat=2)]
    h = (ex + 8 * contrast)[np.ix_(ids, ids)]
    identity = np.eye(4)
    reader_h, supply_h = np.kron(h, identity), np.kron(identity, h)
    release = expm(-1j * np.pi * h)
    # Reorder (reader L,R, supply L,R) into the two transfer seams.
    perm = np.eye(16).reshape((2,) * 4 + (16,)).transpose(0, 2, 1, 3, 4).reshape(16, 16)
    transfer = perm.T @ np.kron(release, release) @ perm
    total = reader_h + supply_h
    np.testing.assert_allclose(transfer.conj().T @ total @ transfer, total, atol=6e-14)
    # Adjacent equal controls see both arms, so there is no net new energy.
    reference = np.diag((1 + z * z * np.array([1, -1, -1, 1])) / 4)
    initial = np.kron(reference, identity / 4)
    final = transfer @ initial @ transfer.conj().T
    np.testing.assert_allclose(final, np.kron(identity / 4, reference), atol=3e-14)
    delta = final - initial
    reader, supply = [float(np.trace(hh @ delta).real) for hh in (reader_h, supply_h)]
    np.testing.assert_allclose(
        [reader, supply], [17 * z * z / 4, -17 * z * z / 4], atol=4e-14
    )
    return dict(
        native_ports=cycle,
        reader_energy_change=reader,
        supply_energy_change=supply,
        summed_ticket_change=(reader + supply) / SCALE,
    )


def shell(k, rate, field, q):
    n = len(rate)
    lo, hi = n // 2 - k, n // 2 + k  # Interior array indices [lo,hi).
    enclosed = float(rate[lo:hi, lo:hi, lo:hi].sum())
    flux = 0.0
    for axis in range(3):
        f = np.moveaxis(field, axis, 0)[:, lo + 1 : hi + 1, lo + 1 : hi + 1]
        flux += q * (f[lo + 1] - f[lo] + f[hi] - f[hi + 1]).sum()
    np.testing.assert_allclose(flux, enclosed, rtol=2e-12, atol=1e-20)
    edges = 24 * k * k
    return dict(
        k=k,
        edges=edges,
        enclosed_tickets=enclosed,
        outward_flux=float(flux),
        inward_recoil_per_edge=float(1.6 * flux / edges),
        radius_times_inward_recoil=float(1.6 * k * flux / edges),
    )


def first_arrivals(n, point, q):
    """Forward first-passage distribution, with the target and six drains absorbing."""
    live = central_load(n)
    hits = live[point]
    live[point] = 0
    history = {}
    for age in range(1025):
        if age in (0, 4, 16, 64, 256, 1024):
            history[age] = float(hits)
        if age == 1024:
            break
        nxt = (1 - 6 * q) * live
        for axis in range(3):
            v, u = np.moveaxis(nxt, axis, 0), np.moveaxis(live, axis, 0)
            v[1:] += q * u[:-1]
            v[:-1] += q * u[1:]
        hits += nxt[point]
        nxt[point] = 0
        live = nxt
    return history


def guide_panel(n, q, z, energy):
    weights = (q,) * 3
    u = green.poisson(central_load(n), weights)[1:-1, 1:-1, 1:-1]
    diagonal = diagonal_visits(n, q)
    hit = u / diagonal
    assert hit.min() > 0 and hit.max() <= 1
    bias = z * hit
    rate = energy / (2 * SCALE) * bias**2
    charge = bias * np.arctanh(bias)
    np.testing.assert_array_less(
        abs(charge - bias**2), bias**4 / (3 * (1 - bias**2)) + 1e-23
    )
    field = green.poisson(rate, weights)
    ks = sorted({1, n // 2, *(k for k in (2, 4, 8, 16, 32, 64) if k <= n // 2)})
    shells = [shell(k, rate, field, q) for k in ks]
    checks = []
    if n <= 34:
        for index in (n // 2 - 1, n // 2 + 2, n - 2):
            point = (index, n // 2 - 1, n // 2 - 1)
            impulse = np.zeros_like(u)
            impulse[point] = 1
            direct = green.poisson(impulse, weights)[1:-1, 1:-1, 1:-1]
            np.testing.assert_allclose(direct[point], diagonal[point], rtol=4e-14)
            np.testing.assert_allclose(
                (central_load(n) * direct).sum() / direct[point], hit[point]
            )
            history = first_arrivals(n, point, q)
            assert history[1024] <= hit[point] + 2e-14
            if n == 10:
                np.testing.assert_allclose(history[1024], hit[point], atol=3e-14)
            checks.append(
                dict(
                    point_zero_based=point, ever_hit=float(hit[point]), by_tick=history
                )
            )
    return dict(
        interior=n,
        max_hit=float(hit.max()),
        max_bias_squared=float((bias**2).max()),
        total_excess_ticket_rate=float(rate.sum()),
        first_passage_checks=checks,
        shells=shells,
    ), rate


def paired_acquisition(rate, steps):
    """Common random tickets cancel identical events. Count every unmatched pulse and visit."""
    rng = np.random.default_rng(271119)
    # At each address this many common uniform tickets serve each treatment arm.
    proposals = 2**45
    live = rng.binomial(proposals, rate)
    released = int(live.sum())
    totals = np.zeros_like(live)
    for tick in range(receiver.parent.H + 1):
        if not live.any():
            break
        totals += live
        nxt = np.zeros_like(live)
        for j, step in enumerate(steps[:-2]):
            take = rng.binomial(live, 1 / (8 - j))
            live -= take
            axis = np.flatnonzero(step)[0]
            a, b = np.moveaxis(nxt, axis, 0), np.moveaxis(take, axis, 0)
            if step[axis] > 0:
                a[1:] += b[:-1]
            else:
                a[:-1] += b[1:]
        live += nxt
    assert not live.any()
    visits = green.poisson(rate / rate.sum(), (1 / 8,) * 3)[1:-1, 1:-1, 1:-1]
    panels = []
    for probe in (1, 2):
        face = (6, 5, 5, 0)
        setting = receiver.Setting("matched", face, probe, (2, 2, 2))
        arms = receiver.count_panel(rng, totals, setting)
        expected_reward = receiver.signed_readout(setting)
        # Check the endpoint probabilities themselves against both local face weights.
        probs = receiver.table(setting)
        from_table = sum(
            (p[words] - 0.5) / 2
            for p, words in zip(probs, receiver.parent.addresses(face), strict=True)
        )
        np.testing.assert_allclose(from_table, expected_reward, atol=2e-14)
        signed = sum(sum(a["right"]) - sum(a["starts"]) / 2 for a in arms) / released
        predicted = float((visits * expected_reward).sum())
        gradient = (visits[5, 4, 4] - visits[6, 4, 4]) / 8
        np.testing.assert_allclose(predicted, -1.6 * gradient)
        # The bound assumes iid starts. Independent Bernoulli proposals at fixed
        # sites are stratified, so use the unconditional per-proposal bound below.
        a = abs(expected_reward) + 2**-33
        remaining = green.poisson(a, (1 / 8,) * 3)[1:-1, 1:-1, 1:-1]
        second = float((visits * (0.25 + 2 * a * remaining)).sum())
        expected_count = proposals * rate.sum()
        x = np.log(4 / 0.01)
        # Bound the independent whole-walk rewards divided by expected releases.
        width = (
            np.sqrt(2 * second * x / expected_count)
            + 2 * receiver.parent.H * x / (3 * expected_count)
            + 120 * 2**-33
            + 1e-11
        )
        observed = signed * released / expected_count
        assert abs(observed - predicted) < width and observed + width < 0
        panels.append(
            dict(
                probe=probe,
                arms=arms,
                signed_per_expected_excess=observed,
                forecast=predicted,
                simultaneous_99_width=float(width),
            )
        )
    return dict(
        proposals_per_address=proposals,
        expected_releases=float(proposals * rate.sum()),
        unmatched_releases=released,
        visits=int(totals.sum()),
        last_drain=tick,
        panels=panels,
    )


def run(sizes=INTERIORS):
    from matched_einstein_limit.collar_identifiability.refinement import admission

    steps, weights = green.measured_steps(None)
    np.testing.assert_array_equal(weights, [1 / 8] * 3)
    source = {}
    for strength in (0, 8):
        panel = admission.material_panel(1, 1, strength)
        record = next(r for r in panel["records"] if r["time"] == 2)
        source[strength] = record
    z = source[8]["d"][-1]
    assert abs(source[0]["d"][-1]) < 2e-14
    np.testing.assert_allclose(
        admission.balanced_charge(
            np.array(source[8]["q"]), np.array(source[8]["d"]), 1
        ),
        z * np.arctanh(z),
    )
    energy, local = local_source(z)
    result = dict(
        input_return=source,
        ordinary_unlike_minus_equal=energy,
        energy_ticket_scale=SCALE,
        parked_archive_energy=16,
        local=local,
        conservative_transfer=conservative_preparation(z),
        guides=[],
    )
    for n in sizes:
        panel, rate = guide_panel(n, weights[0], z, energy)
        result["guides"].append(panel)
        if n == 10:
            result["paired_acquisition"] = paired_acquisition(rate, steps)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    output = HERE / "results.json"
    if args.write:
        output.write_text(json.dumps(result, indent=2) + "\n")
    else:
        from numpy.testing import assert_allclose

        # Numerical libraries can change the last bits of the Green transform.
        def check(old, new):
            if isinstance(old, dict):
                assert old.keys() == new.keys()
                for key in old:
                    check(old[key], new[key])
            elif isinstance(old, list):
                assert len(old) == len(new)
                for left, right in zip(old, new, strict=True):
                    check(left, right)
            elif isinstance(old, float):
                assert_allclose(old, new, rtol=5e-12, atol=2e-16)
            else:
                assert old == new

        check(json.loads(output.read_text()), json.loads(json.dumps(result)))
