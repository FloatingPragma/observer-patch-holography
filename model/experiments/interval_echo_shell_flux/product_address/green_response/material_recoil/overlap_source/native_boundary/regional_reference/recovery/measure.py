"""Flagged/unflagged information loss and native Bell routing on the two-gateway cut."""

import argparse
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag, eigh, expm
from scipy.optimize import brentq

from .. import measure as region
from ...matched_probe import measure as matched
from ...matched_probe import time_four
from ...matched_probe.signal_time import signal as timing
from ....closed_region.measure import layout
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
PSET = (0, 1 / 8, 1 / 4, 1 / 2, 3 / 4, 1)
ANGLES = np.pi * np.arange(5) / 8


def preparation():
    e = np.eye(4)
    b = (e[1] + e[2]) / np.sqrt(2)
    d = (e[1] - e[2]) / np.sqrt(2)
    u = (e[0] - e[3]) / np.sqrt(2)
    bright = np.outer(b, b)
    out = {}
    out["P"] = matched.stocks()["target"][0][1]
    out["E"] = matched.stocks()["target"][1][1]
    for k, theta in enumerate(ANGLES):
        v = np.cos(theta) * d + 1j * np.sin(theta) * u
        out["FG" + str(k)] = 0.25 * bright + 0.75 * np.outer(v, v.conj())
    np.testing.assert_allclose(
        out["FG0"], matched.stocks()["phase_feedback"][1][1], atol=1e-15
    )
    return out, b, d, u


def bypass():
    board = federation(1)
    physical, _, _, _ = layout(2)
    a, b = [physical[j] for j in matched.fp.make_board()[1]]
    todo = deque([a])
    prev = {a: None}
    while b not in prev:
        node = todo.popleft()
        for z in sorted(board.neighbours[node]):
            if z not in prev and {z, node} != {a, b}:
                prev[z] = node
                todo.append(z)
    path = [b]
    while path[-1] != a:
        path.append(prev[path[-1]])
    path.reverse()
    assert set(path) - set(physical)
    hs = matched.source.LAW.seam_term(8)
    trit, _ = matched.source.LAW._load_ops()
    np.testing.assert_allclose(expm(-1j * np.pi * hs), trit, atol=2e-14)
    # On the full one-unit path a relay swap permutes two basis addresses.
    size = len(path)
    perm = np.eye(size)
    for i in range(size - 2):
        perm[[i, i + 1]] = perm[[i + 1, i]]
    original = np.eye(size)[:, [0, -1]]
    for t in ANGLES:
        near = np.eye(size, dtype=complex)
        near[-2:, -2:] = expm(-2j * t * hs[np.ix_([1, 3], [1, 3])])
        relayed = perm.T @ near @ perm @ original
        expected = np.exp(-17j * t) * (
            np.cos(t) * original + 1j * np.sin(t) * original[:, ::-1]
        )
        np.testing.assert_allclose(relayed, expected, atol=1e-14)
    return dict(
        parent_source_ports=[a, b],
        parent_bypass=path,
        outside_induced_pair=[p for p in path if p not in physical],
        pulses_per_dark_route=2 * (size - 2) + 1,
        conditional_probability=0.75,
    )


def margins(a):
    t = a.reshape(3, 4, 3, 4)
    return np.trace(t, axis1=1, axis2=3), np.trace(t, axis1=0, axis2=2)


def embed(r):
    out = np.zeros((12, 12), dtype=complex)
    ix = [4, 8, 6, 10]
    out[np.ix_(ix, ix)] = r
    return out


def reference(beta):
    _, seam, _, h, _, channels, _ = region.folded()
    en, v = eigh(h)
    pr = np.exp(-beta * (en - en.min()))
    thermal = (v * (pr / pr.sum())) @ v.T
    ix = [4 * p + m for p in seam for m in range(4)]
    vac = sum(
        thermal[4 * p : 4 * p + 4, 4 * p : 4 * p + 4]
        for p in range(24)
        if p not in seam
    )
    return block_diag(vac, thermal[np.ix_(ix, ix)]), channels


def modular(r):
    v, u = eigh(r)
    assert min(v) > 0
    return (u * (-np.log(v))) @ u.conj().T


def expectation(r, a):
    return float(np.trace(r @ a).real)


def entropy(r):
    ev = np.linalg.eigvalsh(r)
    assert ev.min() > -1e-13
    ev = ev[ev > 1e-14]
    return float(-ev @ np.log(ev))


def group(theta, d, u):
    # In side/coin order the dark vector is d and the unmatched one u, with
    # the sign register repeated. Every other Bell direction is fixed.
    generator = np.outer(d, u) + np.outer(u, d)
    perm = [0, 2, 1, 3]
    gen = generator[np.ix_(perm, perm)]
    return block_diag(np.eye(4), np.kron(expm(1j * theta * gen), np.eye(2))), gen


def analysis(beta, states, d, u):
    sigma, channels = reference(beta)
    sl, sm = margins(sigma)
    sigmap = np.kron(sl, sm)
    ks = modular(sigma)
    kl = modular(sl)
    km = modular(sm)
    separate = np.kron(kl, np.eye(4)) + np.kron(np.eye(3), km)
    rho = {n: embed(v) for n, v in states.items()}
    out = {}
    for n, r in rho.items():
        rl, rm = margins(r)
        q = np.kron(rl, sm)
        d0 = expectation(r, ks) - entropy(r)
        dload = expectation(rl, kl) - entropy(rl)
        out[n] = dict(
            entropy=entropy(r),
            input_D=d0,
            load_D=dload,
            modular=expectation(r, ks),
            connected_modular=expectation(r, ks - separate),
            thinning=[],
        )
        for p in PSET:
            rp = (1 - p) * r + p * q
            sp = (1 - p) * sigma + p * sigmap
            kp = modular(sp)
            ur = d0 - (expectation(rp, kp) - entropy(rp))
            fr = p * (d0 - dload)
            # Keep the actual flag cost H_B(p) out of this conditional
            # modular expectation. It cancels in all source differences.
            flagged_linear = p * (expectation(r, ks) - expectation(rl, kl))
            unflagged_linear = expectation(r, ks) - expectation(rp, kp)
            ev, v = eigh(kp)
            groups = np.split(
                np.arange(len(ev)), np.flatnonzero(np.diff(ev) > 1e-9) + 1
            )
            outcomes = [
                float(np.trace(v[:, gg].conj().T @ rp @ v[:, gg]).real) for gg in groups
            ]
            if 0 < p < 1:
                flagged_ref = block_diag((1 - p) * sigma, p * sl)
                flagged_rho = block_diag((1 - p) * r, p * rl)
                direct = d0 - (
                    expectation(flagged_rho, modular(flagged_ref))
                    - entropy(flagged_rho)
                )
                np.testing.assert_allclose(direct, fr, atol=3e-14)
            out[n]["thinning"].append(
                dict(
                    p=p,
                    flagged_loss=fr,
                    unflagged_loss=ur,
                    flagged_conditional_modular_difference=flagged_linear,
                    unflagged_modular_difference=unflagged_linear,
                    unflagged_K_levels=[float(ev[gg].mean()) for gg in groups],
                    unflagged_K_probabilities=outcomes,
                    thinned_state_spectrum=np.linalg.eigvalsh(rp).tolist(),
                )
            )
            assert ur + 1e-12 >= fr >= -1e-12
    # All-angle equality is algebraic. The Gibbs reference on the occupied
    # regional block is a linear combination of the four exact native channel
    # projectors. Partial traces of EACH channel are scalar on side and coin.
    _, gen = group(0, d, u)
    twice = np.rint(2 * gen).astype(int)
    for four in channels.values():
        g = np.kron(twice, np.eye(2, dtype=int))
        assert np.array_equal(g @ four, four @ g)
        arr = four.reshape(2, 2, 2, 2, 2, 2)
        mem = np.trace(arr, axis1=0, axis2=3).reshape(4, 4)
        assert np.array_equal(mem, np.kron(np.eye(2, dtype=int), mem[:2, :2]))
        sides = np.trace(four.reshape(2, 4, 2, 4), axis1=1, axis2=3)
        assert np.array_equal(sides, np.eye(2, dtype=int) * sides[0, 0])
    f = rho["FG0"]
    for k, theta in enumerate(ANGLES):
        w, _ = group(theta, d, u)
        r = rho["FG" + str(k)]
        rl, rm = margins(r)
        fl, fm = margins(f)
        np.testing.assert_allclose(w @ f @ w.conj().T, r, atol=2e-14)
        np.testing.assert_allclose(w @ sigma @ w.conj().T, sigma, atol=2e-14)
        np.testing.assert_allclose(w @ sigmap @ w.conj().T, sigmap, atol=2e-14)
        np.testing.assert_allclose(
            w @ np.kron(fl, sm) @ w.conj().T, np.kron(rl, sm), atol=2e-14
        )
        np.testing.assert_allclose(rm, fm, atol=3e-15)
        for a, b in zip(
            out["FG0"]["thinning"], out["FG" + str(k)]["thinning"], strict=True
        ):
            for key in a:
                np.testing.assert_allclose(a[key], b[key], atol=1e-13)
    return dict(
        beta=beta,
        reference_load_marginal=np.real(sl).tolist(),
        reference_memory_marginal=np.real(sm).tolist(),
        states=out,
    )


def propagation(states):
    _, dist, origins, blocks = matched.calibrated()
    base = states["P"]
    flight = {}
    for n, r in states.items():
        rows = []
        for p in origins:
            op = blocks[1, p]["number"] - blocks[1, None]["number"]
            rows.append((matched.trace(r - base, op) @ dist).tolist())
        flight[n] = dict(launches=rows, mean=np.mean(rows, axis=0).tolist())
    # The two neutral hands commute with Q=N_first−source_coin. Bell phase
    # coherences between unmatched words or between matched and unmatched
    # sectors cannot affect them. G's unmatched incoherent corner is the
    # difference between twice E and P. Hence its independent full clock
    # waveform is obtainable from existing physical-count calibration arms.
    phys, _ = timing.get_arms()
    banks, _ = timing.bank(phys)
    ix = {n: timing.NAMES.index(n) for n in timing.NAMES}
    clock = {}
    counts = {}
    shots = 2**30
    rng = np.random.default_rng(501142)
    family = len(states) * 2 * 3 * 2
    epsilon = np.sqrt(np.log(2 * family / 0.01) / (2 * shots))
    q = np.array([1, 0, 0, -1])
    keep = q[:, None] == q[None, :]
    actual = {
        key: [np.zeros((2, 2)), (r * keep).real]
        for key, r in zip(timing.NAMES, states.values(), strict=True)
    }
    for r in states.values():
        np.testing.assert_allclose((r * keep).imag, 0, atol=1e-15)
    direct, direct_check = timing.bank(actual)
    delta_max = 0.0
    for col, (n, r) in enumerate(states.items()):
        cs = []
        for k, b in enumerate(banks):
            f, weights, _ = b
            if n in ("P", "E"):
                w = weights[:, ix["paired" if n == "P" else "erased"]]
            else:
                theta = ANGLES[int(n[2:])]
                # Construct the response for G independently from registered
                # positive bright and dark decomposition: bright clock
                # read is 3P-2F, since P=(B+D)/2, F=(B+3D)/4.
                wp = weights[:, ix["paired"]]
                we = weights[:, ix["erased"]]
                wf = weights[:, ix["feedback"]]
                wg = 0.25 * (3 * wp - 2 * wf) + 0.75 * (2 * we - wp)
                w = np.cos(theta) ** 2 * wf + np.sin(theta) ** 2 * wg
            tt = np.arange(321) / 32
            zz = np.cos(tt[:, None] * f) @ w
            fd, wd, _ = direct[k]
            raw = wd[:, col]
            independent = np.cos(tt[:, None] * fd) @ raw
            delta_max = max(delta_max, float(np.max(abs(independent - zz))))
            j = np.flatnonzero(zz[:-1] * zz[1:] < 0)
            cs.append(
                [
                    float(brentq(lambda t: np.cos(t * f) @ w, tt[i], tt[i + 1]))
                    for i in j
                ]
            )
            for i, t in enumerate(cs[-1]):
                ends = np.array([t - 0.001, t + 0.001])
                ps = (1 + np.cos(ends[:, None] * fd) @ raw) / 2
                yes = rng.binomial(shots, ps)
                oriented = (yes / shots - 0.5) * (-1) ** i
                assert oriented[0] > epsilon and oriented[1] < -epsilon
                counts[f"{n}/{k}/{i + 1}"] = dict(
                    times=ends.tolist(),
                    same=yes.tolist(),
                    opposite=(shots - yes).tolist(),
                )
        clock[n] = cs
    assert delta_max < 1e-9
    intervals = {}
    for n in states:
        events = np.array(clock[n]) - np.array(clock["P"])
        delta = events[0] - events[1]
        intervals[n] = dict(
            near_minus_far_change=delta.tolist(),
            intervals=[(delta + v).tolist() for v in (-0.004, 0.004)],
        )
    cert = time_four.run()
    er = cert["contrasts"]["target"]["radius_billionths"]
    fb = cert["contrasts"]["phase_feedback"]["radius_billionths"]
    # Rational signed-density identity in each Q block, before any flight.
    twiceG = np.array([[3, 0, 0, -3], [0, 1, 1, 0], [0, 1, 1, 0], [-3, 0, 0, 3]])
    np.testing.assert_allclose(twiceG / 8, states["FG4"], atol=1e-15)
    np.testing.assert_array_equal(
        twiceG * keep, 12 * states["E"] - 4 * matched.stocks()["phase_feedback"][1][1]
    )
    lower = F(3 * er[0] - fb[1], 2 * 10**9)
    upper = F(3 * er[1] - fb[0], 2 * 10**9)
    enclosure = dict(
        lower=str(lower),
        upper=str(upper),
        outward_billionths=[(lower * 10**9).__floor__(), (upper * 10**9).__ceil__()],
    )
    return dict(
        times=[0, 1, 4, 16, 64, "time_mean"],
        flight=flight,
        clock_near_far_events=clock,
        clock_contrasts=intervals,
        independent_clock_max_error=max(delta_max, direct_check),
        certified_G_time_four=enclosure,
        count_protocol_sha256=sha256((HERE / "CHECKS.md").read_bytes()).hexdigest(),
        shot_count=shots,
        count_banks=counts,
        simultaneous_probability_half_width=float(epsilon),
    )


def run():
    states, b, d, u = preparation()
    return dict(
        protocol_sha256=sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        preparation=bypass(),
        reference=[analysis(f * np.log(8 / 7), states, d, u) for f in (0.5, 1, 2)],
        probes=propagation(states),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        region.check(json.loads(path.read_text()), result)
    print(result["preparation"])
    for n, r in result["reference"][1]["states"].items():
        print(
            n,
            r["input_D"],
            r["connected_modular"],
            [
                (z["flagged_loss"], z["unflagged_loss"])
                for z in r["thinning"]
                if z["p"] == 0.5
            ],
            result["probes"]["flight"][n]["mean"][2],
        )
