"""Identical written-Hamiltonian comparison of erasure and calibrated ordinary source preparation."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from .. import source
from ..free_probe import measure as fp

HERE = Path(__file__).resolve().parent


def pure(x):
    return np.outer(x, x)


def stocks():
    eye = np.eye(2)
    fair = eye / 2
    paired = np.diag([0, 0.5, 0.5, 0])
    erased = np.eye(4) / 4
    standards = {"target": ((0, paired), (0, erased))}
    feedback = paired.copy()
    feedback[1, 2] = feedback[2, 1] = -0.25
    standards["phase_feedback"] = ((0, paired), (0, feedback))
    for c in (0, 1):
        mem = pure(eye[c])
        lo = np.kron(mem, pure(eye[1 - c]))
        hi = np.kron(mem, fair)
        standards[f"sharp_coin_{c}"] = ((0, lo), (0, hi))
        for s in (1, -1):
            phi = np.array([1, s]) / np.sqrt(2)
            standards[f"fixed_coin_{c}_orbit_{s}"] = (
                (0, lo),
                (0, np.kron(mem, pure(phi))),
            )
    for u in (1, -1):
        mem = pure(np.array([1, u]) / np.sqrt(2))
        for v in (1, -1):
            orb = pure(np.array([1, v]) / np.sqrt(2))
            standards[f"exchange_{u}_{v}"] = (
                (0, np.kron(mem, orb)),
                (0, np.kron(mem, fair)),
            )
    for c in (0, 1):
        mem = pure(eye[c])
        standards[f"mass_coin_{c}"] = (
            (mem, 0),
            (18 / 19 * mem, np.kron(mem, fair) / 19),
        )
    for label, a, b in (
        ("balanced_sharp", "sharp_coin_0", "sharp_coin_1"),
        ("balanced_coherent", "fixed_coin_0_orbit_1", "fixed_coin_1_orbit_-1"),
        ("balanced_mass", "mass_coin_0", "mass_coin_1"),
    ):
        standards[label] = tuple(
            tuple((x + y) / 2 for x, y in zip(t1, t2, strict=True))
            for t1, t2 in zip(standards[a], standards[b], strict=True)
        )
    result = {}
    for name, arms in standards.items():
        out = []
        for r0, r1 in arms:
            mats = [
                np.zeros((d, d)) if np.isscalar(r) else r for d, r in ((2, r0), (4, r1))
            ]
            assert abs(sum(np.trace(w) for w in mats) - 1) < 1e-12
            assert min(np.linalg.eigvalsh(w).min() for w in mats) > -1e-12
            out.append(mats)
        result[name] = out
    return result


def trace(rho, mat):
    return np.einsum("ij,...ji->...", rho, mat).real


def preparations(seam, origins):
    groups = {}
    for n in (0, 1):
        for extra in [None, *origins]:
            ports = [()] if n == 0 else [(p,) for p in seam]
            groups[n, extra] = [
                (tuple(sorted((*w, *(() if extra is None else (extra,))))), m)
                for m in (0, 2)
                for w in ports
            ]
    return groups


def local_energy(board, seam, keys):
    h = source.formation.emission.four_h(8) / 4
    tags = [1 + i for i, (_, d) in enumerate(source.LAW.LABELS) if abs(d) == 1]
    local_ids = [
        9 * tags[m] + source.LAW.SRC.index(tuple(int(p in w) for p in seam))
        for w, m in keys
    ]
    _, projectors = source.energy_projectors(8)
    local = h[np.ix_(local_ids, local_ids)]
    shape = local.shape
    shares = np.zeros((board.n, *shape))
    currents = np.zeros((len(board.seams), *shape), complex)
    # All states in one block have the same fixed spectator outside the
    # written seam, including the extra far port when it is present.
    for ii, (w, _) in enumerate(keys):
        b = np.isin(np.arange(board.n), w)
        shares[:, ii, ii] += 8 * (6 - board.degree) * b
        for edge in board.seams:
            if edge != seam and b[edge[0]] != b[edge[1]]:
                shares[list(edge), ii, ii] += 4.25
    for p in seam:
        shares[p] += local / 2
    # Coherence connects distinct endpoint words only across the written seam.
    a = seam[0]
    occupation = np.array([int(a in w) for w, _ in keys])
    currents[board.seams.index(seam)] = 1j * (
        local * occupation[None, :] - occupation[:, None] * local
    )
    return (
        shares,
        currents,
        np.stack(
            [(proj / den)[np.ix_(local_ids, local_ids)] for _, proj, den in projectors]
        ),
    )


def read_operators(h, bits, ids):
    es, v = eigh(h.toarray())
    coeff = v[ids, :].T
    counts = np.repeat(bits, 4, axis=0)

    def effect(psi):
        return np.einsum("pi,px,pj->xij", psi.conj(), counts, psi, optimize=True)

    timed = np.stack(
        [
            effect(v @ (np.exp(-1j * t * (es - es[0]))[:, None] * coeff))
            for t in fp.TIMES
        ]
    )
    result = {}
    for tol in (1e-8, 1e-10):
        scale = np.split(np.arange(len(es)), np.flatnonzero(np.diff(es) > tol) + 1)
        mean = np.zeros_like(timed[0])
        for ix in scale:
            mean += effect(v[:, ix] @ coeff[ix])
        if result:
            np.testing.assert_allclose(result["mean"], mean, rtol=0, atol=2e-12)
        result = dict(
            mean=mean,
            energies=np.array([es[ix].mean() for ix in scale]),
            weights=np.array([coeff[ix].T @ coeff[ix] for ix in scale]),
        )
    result["number"] = np.concatenate((timed, result.pop("mean")[None]))
    return result


def calibrated():
    board, seam, dist, origins = fp.make_board()
    prep = preparations(seam, origins)
    blocks = {}
    for total in range(3):
        h, words, bits = fp.generator(board, seam, total, 8, True)
        for (n, extra), keys in prep.items():
            if n + int(extra is not None) != total:
                continue
            ids = [4 * words.index(w) + m for w, m in keys]
            ops = read_operators(h, bits, ids)
            shares, current, energy = local_energy(board, seam, keys)
            np.testing.assert_allclose(
                shares.sum(axis=0), h[ids, :][:, ids].toarray(), atol=1e-13
            )
            ops.update(shares=shares, current=current, local_spectrum=energy)
            blocks[n, extra] = ops
    return board, dist, origins, blocks


def aggregate(arms, extra, blocks):
    # Keep energies on the same extensive Hamiltonian, even for a preparation
    # that mixes different conserved load-number sectors.
    obs = {}
    for op in ("number", "shares", "current", "local_spectrum"):
        obs[op] = sum(trace(m, blocks[n, extra][op]) for n, m in enumerate(arms))
    hist = {}
    for n, m in enumerate(arms):
        b = blocks[n, extra]
        for e, w in zip(b["energies"], trace(m, b["weights"]), strict=True):
            key = str(round(e, 8))
            hist[key] = hist.get(key, 0.0) + w
    obs["hist"] = hist
    return obs


def difference(a, b):
    return sum(abs(a.get(e, 0.0) - b.get(e, 0.0)) for e in a.keys() | b.keys()) / 2


def run():
    board, dist, origins, blocks = calibrated()
    data = {}
    target_reference = {
        p: [aggregate(arm, p, blocks)["hist"] for arm in stocks()["target"]]
        for p in origins
    }
    for name, (low, high) in stocks().items():
        lo = aggregate(low, None, blocks)
        hi = aggregate(high, None, blocks)
        probes = []
        for p in origins:
            lx = aggregate(low, p, blocks)
            hx = aggregate(high, p, blocks)
            contrast = (hx["number"] - hi["number"]) - (lx["number"] - lo["number"])
            np.testing.assert_allclose(contrast.sum(axis=1), 0, atol=5e-12)
            np.testing.assert_allclose(
                hx["shares"] - lx["shares"], hi["shares"] - lo["shares"], atol=3e-14
            )
            probes.append(
                dict(
                    port=p,
                    port_response=contrast.tolist(),
                    radial=(contrast @ dist).tolist(),
                    with_extra_energy_TV=[
                        difference(z["hist"], target_reference[p][i])
                        for i, z in enumerate((lx, hx))
                    ],
                    with_extra_occupations=[z["number"].tolist() for z in (lx, hx)],
                )
            )
        np.testing.assert_allclose(lo["current"], 0, atol=0)
        np.testing.assert_allclose(hi["current"], 0, atol=0)
        data[name] = dict(
            source_only_energy_by_port=[z["shares"].tolist() for z in (lo, hi)],
            energy_gain=float(sum(hi["shares"] - lo["shares"])),
            local_sector_weights=[z["local_spectrum"].tolist() for z in (lo, hi)],
            source_only_port_occupations=[z["number"].tolist() for z in (lo, hi)],
            total_energy_distributions=[z["hist"] for z in (lo, hi)],
            probes=probes,
        )
    target = data["target"]
    for name, row in data.items():
        row["total_energy_TV_from_target"] = [
            difference(a, b)
            for a, b in zip(
                row["total_energy_distributions"],
                target["total_energy_distributions"],
                strict=True,
            )
        ]
    return dict(
        protocol_sha256=sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        followup_protocol_sha256=sha256(
            (HERE / "COMMON_BASELINE.md").read_bytes()
        ).hexdigest(),
        times=[*fp.TIMES, "infinite_time_mean"],
        local_energies=[0, 8, 9, 32, 33],
        standards=data,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    d = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(d, allow_nan=False, separators=(",", ":")) + "\n")
    else:
        fp.check(json.loads(path.read_text()), d)
    for name, r in d["standards"].items():
        print(
            name,
            r["energy_gain"],
            [np.round(s["radial"], 6) for s in r["probes"]],
            r["total_energy_TV_from_target"],
        )
