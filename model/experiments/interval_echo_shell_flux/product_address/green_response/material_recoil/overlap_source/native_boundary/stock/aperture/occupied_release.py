"""Free vacancy and excess-load trajectories on continuously occupied written seams."""

import argparse
from hashlib import sha256
from itertools import product
import json
from math import comb, factorial
from fractions import Fraction

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path

from . import thermal_pair as lab
from . import thermal_ordinary as standard

HERE = lab.HERE
SECTORS = list(product(lab.STRENGTHS, repeat=2))
CUTS = (0, 1, 2, 4, 6)
AGES = (1, 4, 16, 64)
N = 2**52
RISK = np.log(2e6 / 0.01)
EPS = np.sqrt(RISK / (2 * N))
CONDITIONAL = np.array([lab.arm_weights(h, SECTORS) for h in range(4)])
MIX = np.array([1, -1, -1, 1])


def apparatus(board, edges, load, clamp, placed):
    """Read the relative single-defect generator against the filled reference."""
    cut = set(map(tuple, edges))
    adj = board.adjacency().toarray().astype(int)
    live = adj.copy()
    marks = []
    for u, v in edges:
        live[u, v] = live[v, u] = 0
        marks.append(np.isin(np.arange(board.n), [u, v]).astype(int))
    active = np.arange(board.n) != clamp
    # Twice the native energy relative to uniform filling. The half exchange
    # and the contrast are both retained when a mobile defect meets a clamp.
    k = lab.motion.K
    diag = live.sum(1) * (1 + 2 * k * (load - 1) ** 2) + 2 * k * (load**2 - 1) * (
        6 - board.degree
    )
    if placed != 1:

        def local(c):
            return (
                int(load != c) - int(1 != c) + 2 * k * ((load - c) ** 2 - (1 - c) ** 2)
            )

        diag += adj[clamp] * (local(placed) - local(1))
    twice = np.diag(diag) - live
    # Independent native-column compression. Subtract each contact's reading
    # on the reference configuration before assembling the generator.
    _, _, raw, _ = lab.contact.stock.kernels()
    ref = np.ones(board.n, dtype=int)
    if clamp >= 0:
        ref[clamp] = placed
    direct = np.diag(4 * k * (load**2 - 1) * (6 - board.degree))
    for u, v in board.seams:
        if (u, v) in cut:
            continue
        vacuum = 3 * ref[u] + ref[v]
        words = np.array([3 * load + ref[v], 3 * ref[u] + load])
        direct[np.ix_([u, v], [u, v])] += raw[np.ix_(words, words)] - raw[
            vacuum, vacuum
        ] * np.eye(2, dtype=int)
    idx = np.ix_(active, active)
    np.testing.assert_array_equal(direct[idx], 2 * twice[idx])
    return twice[idx], np.array(marks)[:, active]


def chart():
    board = lab.motion.cold.federation(0)
    old = json.loads((HERE / "pair_contact.json").read_bytes())
    edges = old["physical_seams"]
    adj = board.adjacency().tolil()
    for u, v in edges:
        adj[u, v] = adj[v, u] = 0
    r = (
        shortest_path(
            adj.tocsr(), indices=np.ravel(edges), unweighted=True, directed=False
        )
        .min(0)
        .astype(int)
    )
    return board, old, edges, r, int(np.flatnonzero(r == 1)[0])


def sector_prior(old):
    table = np.array(old["accepted_weight"])
    p = table.ravel() / table.sum()
    q = np.outer(table.sum(1), table.sum(0)).ravel() / table.sum() ** 2
    return p, q


def material_ops(h, marks, inside):
    o = np.diag(inside.astype(float))
    current = h @ o - o @ h
    plain = current @ h - h @ current
    terms = np.array([current * m - current * m[:, None] for m in marks])
    assert np.max(abs(terms[:, ~(inside[:, None] ^ inside)])) == 0
    return plain, terms


def exact_onset(twice, marks, rad, max_order=30):
    """Integer product rule for the derivatives, before averaging any phases."""
    n = len(rad)
    off = twice.copy()
    np.fill_diagonal(off, 0)
    rows, cols = np.nonzero(off)
    assert np.all(off[rows, cols] == -1)
    start = np.ones(n, dtype=object)
    mixed = np.zeros((max_order // 2 + 1, len(CUTS)), dtype=object)
    coeff = {0: -2, 32: 1, 33: 1}
    width = 0
    hopping = np.abs(off).sum(1)
    for a, b in SECTORS:
        diag = (twice.diagonal() + 2 * (a * marks[0] + b * marks[1])).astype(object)
        width = max(width, int(max(diag + hopping) - min(diag - hopping)))
        vectors = [start]
        for _ in range(max_order):
            nxt = diag * vectors[-1]
            np.subtract.at(nxt, rows, vectors[-1][cols])
            vectors.append(nxt)
        for m in range(2, max_order + 1, 2):
            density = sum(
                ((-1) ** (m // 2 + j) * comb(m, j)) * vectors[j] * vectors[m - j]
                for j in range(m + 1)
            )
            assert density.sum() == 0
            for k, r in enumerate(CUTS):
                mixed[m // 2, k] += coeff[a] * coeff[b] * density[rad <= r].sum()
    spectrum_width = Fraction(width, 2)
    out = {}
    for k, r in enumerate(CUTS):
        nonzero = np.flatnonzero(mixed[:, k])
        if len(nonzero):
            m = 2 * int(nonzero[0])

            def coefficient(j):
                return Fraction(mixed[j // 2, k], 14 * 16 * n * 2**j * factorial(j))

            leading = coefficient(m)
            for exponent in range(24):
                tick = Fraction(1, 2**exponent)
                remainder = sum(
                    abs(coefficient(j)) * tick ** (j - m)
                    for j in range(m + 2, max_order + 1, 2)
                )
                # The next odd derivative at zero vanishes. Center the generator
                # spectrally for the global norm of the even Taylor remainder.
                omitted = max_order + 2
                remainder += (
                    spectrum_width**omitted
                    * tick ** (omitted - m)
                    / (14 * factorial(omitted))
                )
                if remainder < abs(leading):
                    break
            assert remainder < abs(leading)
            out[str(r)] = dict(
                order=m,
                derivative=str(leading * factorial(m)),
                certified_through=str(tick),
                signed_at_endpoint=[
                    float((leading - remainder) * tick**m),
                    float((leading + remainder) * tick**m),
                ],
            )
        else:
            out[str(r)] = dict(zero_through=max_order)
    return out


def evolve(twice, marks, rad):
    n = len(rad)
    h0 = twice / 2
    regions = [rad <= r for r in CUTS]
    ops = [material_ops(h0, marks, v) for v in regions]
    for k in range(1, len(ops)):
        assert not ops[k][1].any()
    position = []
    response = []
    interference = []
    stamp = 0.0
    start = np.ones(n) / np.sqrt(n)
    for a, b in SECTORS:
        h = h0 + np.diag(a * marks[0] + b * marks[1])
        e, v = np.linalg.eigh(h)
        # Cross-check the finite eigenprogram and its probability readouts.
        stamp = max(
            stamp,
            float(np.linalg.norm(h @ v - v * e)),
            float(np.linalg.norm(v.T @ v - np.eye(n))),
        )
        amps = v.T @ start
        waves = v @ (amps[:, None] * np.exp(-1j * np.outer(e, AGES)))
        p = abs(waves) ** 2
        avg = (v * amps) @ np.sinc((e[:, None] - e) * AGES[-1] / np.pi)
        avg = np.sum(avg * (v * amps), axis=1)
        p = np.column_stack([p, avg])
        assert p.min() > -1e-12
        np.testing.assert_allclose(p.sum(0), 1, atol=2e-12, rtol=0)
        position.append(
            [np.bincount(rad, np.maximum(q, 0), minlength=max(rad) + 1) for q in p.T]
        )
        hpsi = h @ waves
        acceleration = 2 * (abs(hpsi) ** 2 - (waves.conj() * (h @ hpsi)).real)
        acc = []
        direct_read = []
        for inside, (base, terms) in zip(regions, ops, strict=True):
            readings = np.stack(
                [
                    np.sum(waves.conj() * (op @ waves), axis=0).real
                    for op in (base, *terms)
                ]
            )
            np.testing.assert_allclose(
                np.array([1, a, b]) @ readings, acceleration[inside].sum(0), atol=8e-12
            )
            norm, flags = lab.probe.instrument(
                sparse.csr_matrix(base + a * terms[0] + b * terms[1]), waves
            )
            signed = norm * (flags[:, 0] - flags[:, 1])
            np.testing.assert_allclose(signed, acceleration[inside].sum(0), atol=8e-12)
            direct_read.append(signed)
            acc.append(readings)
        response.append(acc)
        interference.append(direct_read)
    assert stamp < 2e-10
    return np.array(position), np.array(response), np.array(interference), ops, stamp


def counts(rng, probs):
    p = np.maximum(probs, 0)
    return rng.multinomial(N, p / p.sum()).tolist()


def panel(twice, marks, rad, p, q, delta, rng):
    pos, acc, direct_read, ops, stamp = evolve(twice, marks, rad)
    prior = np.array([p, q]) @ CONDITIONAL
    sector = np.array(SECTORS)
    conditional = np.einsum("hs,stx->htx", CONDITIONAL, pos)
    sharp_tables = np.concatenate([conditional, np.einsum("hs,stx->htx", prior, pos)])
    bank = np.array([[counts(rng, v) for v in row] for row in sharp_tables])
    source_shift = {}
    for t, tag in enumerate((*map(str, AGES), "uniform_0_64")):
        by_cut = {}
        for r in CUTS:
            freq = bank[:, t, : r + 1].sum(1) / N
            predicted = lab.probe.mix_bands(
                delta, MIX @ freq[:4] + np.array([-4, 4]) * EPS
            )
            observed = freq[5] - freq[4] + np.array([-2, 2]) * EPS
            assert max(predicted[0], observed[0]) <= min(predicted[1], observed[1])
            by_cut[str(r)] = dict(predicted=predicted.tolist(), held=observed.tolist())
        source_shift[tag] = by_cut
    # Exchange head registers between independent copies after the flight.
    # The averaged coefficients are exactly what the retained marginal gives.
    same = acc[:, :, 0] + np.einsum("se,sket->skt", sector, acc[:, :, 1:])
    baseline = np.einsum("hs,skt->hkt", prior, same)
    reset = []
    for w in prior:
        averaged = np.einsum("s,sket->ket", w, acc)
        reset.append(np.einsum("e,ket->kt", np.r_[1, w @ sector], averaged))
    delta_association = baseline - np.array(reset)
    np.testing.assert_allclose(delta_association[:, 1:], 0, atol=5e-14)
    # Terminal coherence flags read the complete acceleration of the freely
    # evolved preparation; the material operator uses only native hopping,
    # the fixed chart and at most two hops inside or outside each cut.
    acc_tables = {}
    for k, r in enumerate(CUTS):
        base, terms = ops[k]
        largest = max(
            np.max(np.sum(abs(base + a * terms[0] + b * terms[1]), axis=1))
            for a, b in SECTORS
        )
        # Signed two-output read: split unused probability equally after
        # the signed interference measurement. This keeps every trial.
        probs = (1 + direct_read[:, k] / largest) / 2
        assert probs.min() >= 0 and probs.max() <= 1
        sharp = CONDITIONAL @ probs
        banks = rng.binomial(N, np.vstack([sharp, prior @ probs]))
        freq = banks / N
        diff = 2 * largest * (freq[5] - freq[4])
        held = np.array([diff - 4 * largest * EPS, diff + 4 * largest * EPS])
        measured_mixed = 2 * largest * (MIX @ freq[:4])
        forecast = [
            lab.probe.mix_bands(delta, m + np.array([-8, 8]) * largest * EPS).tolist()
            for m in measured_mixed
        ]
        assert all(
            max(f[0], lo) <= min(f[1], hi)
            for f, lo, hi in zip(forecast, *held, strict=True)
        )
        acc_tables[str(r)] = dict(
            norm=float(largest),
            positive_flags=banks.tolist(),
            held=held.T.tolist(),
            predicted=forecast,
        )
    shuffle = {}
    for k, r in enumerate(CUTS):
        norm = acc_tables[str(r)]["norm"]
        # Fresh source heads sampled from the *same arm* after the flight.
        # Discarding the association preserves both material and head marginals.
        released = rng.binomial(N, (1 + np.array(reset)[:, k] / norm) / 2)
        retained = np.array(acc_tables[str(r)]["positive_flags"])[4:]
        change = 2 * norm * (retained - released) / N
        err = 4 * norm * EPS
        np.testing.assert_array_less(abs(change - delta_association[:, k]), err)
        shuffle[str(r)] = dict(
            shuffled_positive_flags=released.tolist(),
            retained_minus_shuffled=np.stack(
                [change - err, change + err], axis=-1
            ).tolist(),
        )
    return dict(
        shuffle=shuffle,
        shell_counts=bank.tolist(),
        ages=[*AGES, "uniform_0_64"],
        cuts=list(CUTS),
        releases=source_shift,
        acceleration=acc_tables,
        reset_acceleration_difference=(np.array(reset)[1] - reset[0]).tolist(),
        association_contribution=delta_association.tolist(),
        spectral_residual=stamp,
        exact_onset=exact_onset(twice, marks, rad),
    )


def run():
    board, old, edges, rad, clamp = chart()
    # Check all output columns and the commuting vacuum/occupation head
    # blocks, including phase spectra originally read by the ordinary clock.
    lab.contact.blocks()
    delta = lab.probe.log_factor(old["source"])
    p, q = sector_prior(old)
    rng = np.random.default_rng(6993126)
    output = {}
    for load, name in ((0, "vacancy"), (2, "excess")):
        panels = {}
        for label, station, c in (
            ("free", -1, 1),
            ("reserved", clamp, 1),
            ("ordinary", clamp, 2),
        ):
            twice, marks = apparatus(board, edges, load, station, c)
            panels[label] = panel(
                twice, marks, rad[np.arange(board.n) != station], p, q, delta, rng
            )
        ratios = {}
        for t in (*map(str, AGES), "uniform_0_64"):
            idx = panels["ordinary"]["ages"].index(int(t) if t.isdigit() else t)
            row = {}
            for r in CUTS:
                cal = (
                    sum(panels["ordinary"]["shell_counts"][4][idx][: r + 1]) / N
                    - sum(panels["reserved"]["shell_counts"][4][idx][: r + 1]) / N
                )
                row[str(r)] = standard.quotient(
                    panels["reserved"]["releases"][t][str(r)]["held"],
                    [cal - 2 * EPS, cal + 2 * EPS],
                )
            ratios[t] = row
        output[name] = dict(panels=panels, ordinary_equivalent=ratios)
    return dict(
        sha256={
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in ("OCCUPIED_RELEASE.md", "pair_contact.json")
        },
        measured_displacement=delta.tolist(),
        ordinary_clamp=clamp,
        shots_per_bank=N,
        source_edges=edges,
        results=output,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    out = run()
    path = HERE / "occupied_release.json"
    if opt.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), out
        )
    for label, v in out["results"].items():
        print(
            label,
            v["panels"]["free"]["exact_onset"],
            v["ordinary_equivalent"]["uniform_0_64"],
        )
