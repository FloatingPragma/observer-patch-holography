"""Independent finite-budget, transfer-clock, and post-run fork checks (pytest -s)."""

from decimal import Decimal as D
from fractions import Fraction as F
from functools import cache
from hashlib import sha256
from itertools import combinations, product
from math import log, pi, sqrt

import mpmath as mp
import numpy as np
import pytest
from scipy.optimize import brentq

from . import raw_verify as raw


N = 2**39
T = np.exp2(np.arange(12))/256
KS = tuple(k*s for k in (.5, 1, 2, 4, 8, 16, 64, 128) for s in (1, -1))+(0,)
RADIUS = sqrt(log(2*204/.01)/(2*N))+5.1e-14


def branch_error(record, admitted, fail):
    assert record["failures"].keys() == fail.keys()
    actual = [D(record["accepted"]), *[D(record["failures"][key]) for key in fail]]
    relative = D(0)
    for value, expected in zip(actual, (admitted, *fail.values()), strict=True):
        assert value.is_finite() and value >= 0 and bool(value) == bool(expected)
        if expected:
            relative = max(relative, abs(value/expected-1))
    assert abs(sum(actual)-1) < D("3e-16")
    return relative


def test_420_branch_masses_without_binary64_underflow():
    rows = raw.read_json(raw.HERE/"budgets.json")["budgets"]
    cases = list(product((20, 80, 320), (1, 2), raw.CAPS, ("matched", "bare_128")))
    assert [(r["cells"], r["prefix"], r["cap"], r["arm"]) for r in rows] == cases
    with mp.workdps(90):
        # The bare vacuum two-level pulse at 128 has gap 128 sqrt(2).
        bare = D(str(mp.sin(mp.pi/mp.sqrt(8))**2/2))
    full = raw.fuel(128, False, 0)
    worst = [D(0)]*3
    for r in rows:
        E = len(raw.board((20, 80, 320).index(r["cells"]))[0])
        args = E, r["prefix"], r["cap"], -1
        if r["arm"] == "matched":
            worst[0] = max(worst[0], branch_error(r, *raw.ledger(*args, .5)))
        else:
            for j, p in enumerate((bare, full), 1):
                worst[j] = max(worst[j], branch_error(r, *raw.ledger(*args, p)))
        if r["cells"] == 320 and r["prefix"] == 2 and r["cap"] in (1, 8) and r["arm"] == "matched":
            assert D(r["accepted"]).adjusted() == {1: -3037, 8: -18}[r["cap"]]
            with pytest.raises(AssertionError):
                branch_error(dict(r, accepted="0E+0"), *raw.ledger(*args, .5))
    print("420 budget masses, max relative differences (matched, bare analytic, bare 198-state):",
          ", ".join(f"{x:.3E}" for x in worst))
    assert worst[0] < D("6e-16") and worst[1] < D("2e-12") and worst[2] < D("4e-11")


def sampler(i, k):
    # H/k on the transmitting edge is [[1,-1],[-1,1]]/2, from the native write.
    with mp.workdps(60):
        p = [(1-mp.cos(mp.mpf(k)*float(t)))/2 for t in T]
        frozen = [float(mp.nint(x*10**13)/10**13) for x in p]
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([44130, i])))
    return rng.binomial(N, frozen).tolist(), np.array([float(x) for x in p])


def rising_bank(counts):
    """Invert the unaliased prefix by root finding, retaining every preceding constraint."""
    assert len(counts) == len(T) and all(type(c) is int and 0 <= c <= N for c in counts)
    low, high = 0., 128.  # Includes the zero arm. This is the archive's conservative envelope.
    for t, hits in zip(T, counts):
        if high*t > pi:
            return None
        p = hits/N

        def root(y):
            y = min(max(0, y), 1)
            return brentq(lambda x: np.sin(t*x/2)**2-y, 0, pi/t, xtol=1e-14)

        low, high = max(low, root(p-RADIUS)), min(high, root(p+RADIUS))
        if low > high:
            return None
        if 8*hits > N:
            return np.clip(root(p), low, high), low, high, t
    return None


def test_every_clock_count_seed_and_simultaneous_rising_interval():
    archive = raw.read_json(raw.HERE/"clock.json")
    assert archive["protocol_sha256"] == sha256((raw.HERE/"CLOCK_TEST.md").read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    assert archive["confidence_radius"] == pytest.approx(RADIUS, rel=2e-15, abs=0)
    assert archive["shots_per_appointment"] == N and archive["appointments"] == list(T)
    assert [r["strength_in_simulator"] for r in archive["runs"]] == list(KS)
    np.testing.assert_array_equal(raw.local()[4][:, [4, 2]], np.eye(6, dtype=int)[:, [4, 2]] @ [[1, -1], [-1, 1]])
    max_error, max_empirical = 0., 0.
    for i, (row, k) in enumerate(zip(archive["runs"], KS, strict=True)):
        counts, p = sampler(i, k)
        assert row["counts"] == counts
        max_empirical = max(max_empirical, np.max(abs(np.asarray(counts)/N-p)))
        reading = rising_bank(counts)
        if k == 0:
            assert reading is None and counts == [0]*12
            assert row == {"strength_in_simulator": 0, "counts": counts, "calibration": None}
        else:
            value, low, high, t = reading
            assert low <= abs(k) <= high and t == 1/abs(k)
            max_error = max(max_error, np.max(abs(np.subtract(reading, row["calibration"]))))
            if i == 0:
                # Altering a bank after the stop cannot change the decoder; it does break custody.
                bad = [*counts[:-1], counts[-1]+1]
                assert rising_bank(bad) == reading and bad != sampler(i, k)[0]
                assert counts != sampler(i+1, k)[0]
    assert max_empirical < RADIUS and max_error < 2e-12
    assert rising_bank([N]*12) is None  # Outside the first unaliased interval.
    print(f"204 counts exact; max clock interval discrepancy {max_error:.3g}, "
          f"max empirical probability deviation {max_empirical:.3g} (radius {RADIUS:.3g})")


@cache
def spectrometer():
    tree = raw.component(0)
    donor = raw.read_json(raw.DONOR/"analyzers.json")
    assert sha256((raw.DONOR/"analyzers.json").read_bytes().replace(b"\r\n", b"\n")).hexdigest() == (
        "bca73380e3f4960b8b4aa1268023573845e9f40463c8545af35ed550e7f17a69")
    theta = np.array([p["angles"] for p in donor["pulses"]])
    assert theta.shape == (32, 96, len(donor["ports"])) and np.isfinite(theta).all()
    # A singleton eigenspace fixes its analyzer up to sign. At the only repeated
    # eigenvalue use the exact signed differences in the historical chart order.
    k2, geom = np.linalg.eigh(tree.lap)
    repeated = next(b for b in raw.bands(k2) if len(b) > 1)
    selected = raw.spectroscopy()
    assert repeated.tolist() == [x for x, _ in selected] == [14, 15]
    for x, (a, b) in selected:
        geom[:, x] = 0
        geom[[tree.words.index(a), tree.words.index(b)], x] = [1/np.sqrt(2), -1/np.sqrt(2)]
    np.testing.assert_allclose(geom.T @ geom, np.eye(32), rtol=0, atol=2e-12)
    physical = np.zeros((32, 96, 32))
    for j, p in enumerate(donor["ports"]):
        physical[:, :, tree.sites == p] = theta[:, :, j, None]
    values, modes = np.linalg.eigh(tree.h2)
    return tree, donor["step"], physical, values, modes, geom


def drive(ratio, *, wrong_phase=False, silent_column=None):
    """One trajectory per historical pulse, integrating the free term in its eigenbasis."""
    _, step, waits, w, V, _ = spectrometer()
    U = (V*np.exp((-4j*step*ratio)*w)) @ V.T
    out = np.zeros((32, 32), complex)
    for j, word in enumerate(waits):
        v = np.eye(32, dtype=complex)[0]
        if j != silent_column:
            phases = np.exp(-1j*(abs(ratio) if wrong_phase else ratio)*word)
            for gates in phases:
                v = gates*(U @ v)
        out[:, j] = v
    return out


def test_512_words_simultaneous_fidelity_and_false_green_controls():
    tree, step, waits, w, V, rays = spectrometer()
    base = drive(1)
    # Free segment and diagonal port phase: Lipschitz constants for 96 products.
    costs = 96*step*np.linalg.norm(4*tree.h2, ord=2) + abs(waits).max(axis=2).sum(axis=1)
    leakage = np.linalg.norm(base-rays*np.sum(rays*base, axis=0), axis=0)
    threshold = np.min((sqrt(.02)-leakage-1e-8)/costs)
    err_fid = err_home = err_bound = actual_ratio = certified_ratio = norm_drift = err_state_bound = 0.
    min_bound, min_fid, min_home = 1., 1., 1.
    for row, k in zip(raw.read_json(raw.HERE/"clock.json")["runs"][:-1], KS[:-1], strict=True):
        val, low, high, _ = rising_bank(row["counts"])
        radius = max(val-low, high-val)/val
        certified_ratio, actual_ratio = max(certified_ratio, radius), max(actual_ratio, abs(abs(k)/val-1))
        actual = drive(k/val)
        distance = np.linalg.norm(actual-(base if k > 0 else base.conj()), axis=0)
        errors = radius*costs
        assert np.all(distance <= errors+2e-11) and radius < threshold
        err_state_bound = max(err_state_bound, abs(max(errors)-row["maximum_state_error"]))
        bounds = 1-(leakage+errors+1e-8)**2
        fid = abs(np.sum(rays*actual, axis=0))**2
        assert min(bounds) > .98 and np.all(fid+2e-11 >= bounds)
        assert len(row["fidelity"]) == 32
        err_fid = max(err_fid, max(abs(fid-row["fidelity"])))
        err_bound = max(err_bound, abs(min(bounds)-row["minimum_bound"]))
        later = V @ (np.exp(-.5j*37*k*w)[:, None]*(V.T @ actual))
        drift = max(np.max(abs(np.sum(abs(s)**2, axis=0)-1)) for s in (actual, later))
        norm_drift = max(norm_drift, drift)
        # This diagnostic depends on the propagator. Each route must satisfy the allowance.
        assert 0 <= row["maximum_column_norm_squared_drift"] < 5e-12
        home = np.sum(abs(later[np.isin(tree.sites, [9, 11])][:, [14, 15]])**2, axis=0)
        err_home = max(err_home, max(abs(home-row["returns_37"])))
        min_bound, min_fid, min_home = min(min_bound, min(bounds)), min(min_fid, min(fid)), min(min_home, min(home))
        if k == -8:
            # Real targets license complex conjugating both the free and addressed terms.
            assert min(abs(np.sum(rays*drive(k/val, wrong_phase=True), axis=0))**2) < .1
            damaged = abs(np.sum(rays*drive(k/val, silent_column=0), axis=0))**2
            assert min(damaged[[14, 15]]) > .98 and damaged[0] < .04
    print(f"512 conditional fidelities: min {min_fid:.12f}, certified {min_bound:.12f}, "
          f"min home {min_home:.12f}; max discrepancies fidelity {err_fid:.3g}, home {err_home:.3g}, "
          f"floor {err_bound:.3g}, state bound {err_state_bound:.3g}; norm drift {norm_drift:.3g}")
    print(f"Relative magnitude error {actual_ratio:.6g}, certified radius {certified_ratio:.6g}, "
          f"maximum radius sufficient for .98 {threshold:.6g}")
    assert err_fid < 5.2e-11 and err_home < 4e-12 and err_bound < 2e-12
    assert err_state_bound < 2e-11 and norm_drift < 4e-12


def test_post_run_fork_census_and_minima():
    donor = {tuple(pair) for _, pair in raw.spectroscopy()}
    maximum = 0.
    forks = raw.read_json(raw.HERE/"budgets.json")["forks"]
    assert [r["cells"] for r in forks] == [20, 80, 320]
    for level, rec in enumerate(forks):
        tree = raw.component(level)
        _, deg, _ = raw.board(level)
        propagator = raw.propagation(level, 8, 37)[0]
        # Edges are exact native column transitions and are keyed by the one coin that changes.
        children, transfers = {}, {}
        for i, (_, coins) in enumerate(tree.states):
            children[i] = []
            for j in np.flatnonzero(tree.h2[:, i] == -1):
                changes = tree.states[j][1]
                if len(changes) != len(coins)+1:
                    continue
                children[i].append(j)
                e, = changes ^ coins
                transfers.setdefault(e, []).append((i, j))
        count = eigen = below = 0
        other = []
        for i, sons in children.items():
            if deg[tree.sites[i]] != 5:
                continue
            for a, b in combinations(sons, 2):
                pair = tree.words[a], tree.words[b]
                contacts = [next(iter(tree.states[j][1] ^ tree.states[i][1])) for j in (a, b)]
                ret = []
                for sign in (-1, 1):
                    v = np.eye(len(tree.states), dtype=int)[i]
                    for q, yes in zip((0, 1, 2, 0, 2), (sign == -1, False, False, False, False)):
                        hit = np.zeros_like(v)
                        if q == 2:
                            hit[tree.sites == tree.sites[i]] = 2*v[tree.sites == tree.sites[i]]
                        else:
                            for u, s in transfers[contacts[q]]:
                                hit[u] = hit[s] = v[u]+v[s]
                        v = hit if yes else 2*v-hit
                    sharp = np.zeros_like(v)
                    sharp[[a, b]] = [1, sign]
                    np.testing.assert_array_equal(v, sharp*v[a])
                    assert abs(v[a]) == 8  # denominator 2**5 from five queries
                    applied = tree.h2 @ sharp
                    if sign == -1:
                        exact = np.array_equal(applied, sharp*applied[a])
                        eigen += int(exact)
                    ports = np.isin(tree.sites, tree.sites[[a, b]])
                    ret.append(sum(abs((propagator @ sharp)[ports])**2)/2)
                count += 1
                if pair not in donor:
                    below += int(ret[0] < .94)
                    other.append((ret[0], ret[1], pair, exact))
        assert (rec["cells"], rec["compiled_both_signs"], rec["minus_exact_eigenrays"], rec["other_forks_below_94"]) == (
            20*4**level, count, eigen, below)
        lowest = min(other)
        witness = rec["lowest_other"]
        assert witness["histories"] == [list(p) for p in lowest[2]] and witness["donor"] is False
        assert witness["minus_eigenray"] == lowest[3]
        maximum = max(maximum, abs(lowest[0]-witness["minus_conditional_home"]), abs(lowest[1]-witness["plus_conditional_home"]))
        edges = len(raw.board(level)[0])
        for target in (F(47, 50), F(49, 50)):
            pair = [next(n for n in range(1, 97) if 1-F(7, 8)**n-edges*failure**n-F(2+n, 2**n) >= target)
                    for failure in (F(1, 2), F(3, 4))]
            assert pair == rec["minimum_caps_from_uniform_bounds"][str(target)]
    print(f"Post-run fork minima: max absolute discrepancy {maximum:.3g}")
    assert maximum < 3e-12
