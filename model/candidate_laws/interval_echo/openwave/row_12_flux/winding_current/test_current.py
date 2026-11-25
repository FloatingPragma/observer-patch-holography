import json

import numpy as np
import pytest
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import expm_multiply

from . import band, calibrate, measure as cur


def sec(n):
    return cur.BinarySector(cur.belt.ring(10), n)


def doubled(sector, cut):
    """Assemble the independent two-pointer Hamiltonian directly contact by contact."""
    h = np.zeros((2 * len(sector.states),) * 2)
    for i, word in enumerate(sector.states):
        occupied = set(word)
        for a, b in sector.car.seams:
            charge = int(a in occupied) + int(b in occupied)
            crosses = abs((a - cut) % 10 - (b - cut) % 10) > 5
            for pointer in [0, 1]:
                ind = 2 * i + pointer
                if charge == 1:
                    j = sector.index[tuple(sorted(occupied ^ {a, b}))]
                    h[ind, ind] += 0.5
                    h[2 * j + int(pointer != crosses), ind] -= 0.5
                if charge == 2:
                    h[ind, ind] -= 2 * cur.KAPPA
    return h


@pytest.mark.parametrize("n", [3, 4])
def test_real_cut_instrument_and_translated_cuts(n):
    sector = sec(n)
    source = int(cur.belt.selected(sector)[0])
    vacuum = np.zeros(len(sector.states) * 2)
    vacuum[2 * source] = 1
    base = cur.probabilities(sector, "native", 0)[1]
    for cut in [0, 4]:
        state = expm_multiply(csc_matrix(-1j * doubled(sector, cut)), vacuum)
        direct = abs(state.reshape(-1, 2)) ** 2
        theoretical = cur.probabilities(sector, "native", cut)[1]
        np.testing.assert_allclose(direct, theoretical, atol=1e-12, rtol=0)
        shifted = base.copy()
        for idx, word in enumerate(sector.states):
            if sum(x < cut for x in (*word, *sector.states[source])) % 2:
                shifted[idx] = shifted[idx, ::-1]
        np.testing.assert_allclose(direct, shifted, atol=1e-12, rtol=0)


@pytest.mark.parametrize("n,current,paths", [(3, -1, 18), (4, 2, 32)])
def test_sharp_tour_and_first_actual_odd_amplitude(n, current, paths):
    sector = sec(n)
    sharp = cur.belt.selected(sector)
    source = int(sharp[0])
    edges = sector.K[np.ix_(sharp, sharp)]
    assert np.all(
        np.count_nonzero(np.triu(edges, 1) + np.triu(edges, 1).T, axis=1) == 2
    )
    for cut in range(10):
        tr = cur.tour(sector, cut)
        assert len({r[0] for r in tr}) == 10
        assert sum(r[2] for r in tr) == current
        assert sum(r[3] for r in tr) == 3
    # Derivatives of the direct propagator. Up to order four no odd cycle
    # exists; at order five all shortest-hop amplitudes have the same sign.
    adj = np.rint(-2 * doubled(sector, 0)).astype("int64")
    v = np.zeros(len(adj), dtype="int64")
    v[2 * source] = 1
    for i in range(5):
        v = adj @ v
        assert v[2 * source + 1] == (paths if i == 4 else 0)


def test_all_calibrations_and_uncharged_controls():
    receipt = json.loads((cur.HERE / "calibration.json").read_text())
    assert len(receipt["banks"]) == 1800
    assert len({tuple(r["edge"]) for r in receipt["banks"]}) == 20
    words = receipt["local_words"]
    for r in receipt["banks"]:
        a, b = r["edge"]
        cut = r["cut"]
        tagged = abs((a - cut) % 10 - (b - cut) % 10) > 5
        total = 0
        transferred = 0
        for j, bit, count in r["occupied_bins"]:
            moved = r["source"] != words[j]
            total += count
            transferred += count * moved
            assert bit == int(moved and tagged)
        assert total == calibrate.SHOTS
        if r["source"] != [0, 0]:
            allowed = {0.0: 0.0, 0.5: 0.5, 1.0: 1.0}[r["wait_pi_units"]]
            assert abs(transferred / total - allowed) < 0.02
    with np.load(cur.HERE / "counts.npz") as f:
        for n in (3, 4):
            bank = f[str(n)]
            assert np.all(bank.sum(axis=(-1, -2)) == cur.SHOTS)
            assert bank[cur.MODES.index("unwired"), :, :, :, 1].sum() == 0
            assert bank[cur.MODES.index("parked"), :, :, :, 1].sum() == 0
        source = cur.belt.selected(sec(4))[0]
        assert np.all(f["4"][1, :, :, source, 1] == 0)
        assert np.all(f["4"][0, :, 1:, source, 1] > 0)


def test_exact_spectral_enclosures_and_clock_histograms():
    r = json.loads((cur.HERE / "band.json").read_text())
    with np.load(cur.HERE / "clock.npz") as banks:
        for row in r["rows"]:
            n, mode = row["n"], row["mode"]
            sector = sec(n)
            h = sector.hamiltonian(0 if mode == "exchange_only" else cur.KAPPA)
            ids = cur.belt.selected(sector) if mode == "sharp" else np.arange(len(h))
            for b in [0, 1]:
                hb = h * cur.counter_signs(sector, 0) if b else h
                _, _, e, cert = band.resolution(sector, hb[np.ix_(ids, ids)])
                recorded = row["branches"][b]["certificate"]
                # Eigenspaces may have a different basis on another LAPACK;
                # both exact-integer matrix certificates enclose the same levels.
                np.testing.assert_allclose(
                    np.array(recorded["eigenvalue_numerators"])
                    / recorded["dyadic_scale"],
                    e,
                    rtol=0,
                    atol=recorded["energy_error"] + cert["energy_error"],
                )
                bins, counts = banks[f"{n}_{mode}_{b}"]
                lo, hi = row["clock_window"]
                assert sum(counts) == band.SHOTS
                assert counts[(lo <= bins) & (bins <= hi)].sum() == row["in_window"][b]
            if n == 4 and mode == "native":
                assert row["certified_lower"] > 0.0085
                assert min(x["gap_below"] for x in row["branches"]) > 13.96
            if n == 4 and mode == "sharp":
                np.testing.assert_allclose(
                    row["branches"][0]["energy"],
                    row["branches"][1]["energy"],
                    atol=1e-12,
                )
                assert row["certified_lower"] == 0


def test_clock_sampler_uses_finite_fourier_kernel(monkeypatch):
    rng = np.random.default_rng(31)
    monkeypatch.setattr(band, "SHOTS", 2**16)
    x = 12.37
    bins, count = band.qpe_samples(np.array([x * band.STEP]), np.ones(1), rng)
    # For a fractional frequency the two closest exact Fourier probabilities
    # appear on their respective detector bins.
    for j in [12, 13, 14, 11]:
        theory = (np.sinc(x - j) / np.sinc((x - j) / band.M)) ** 2
        assert abs(count[bins == j][0] / band.SHOTS - theory) < 0.009
