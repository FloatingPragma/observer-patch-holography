"""Development observations on the unguarded inherited electrical material."""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from math import acos, asin, cos, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.signal import find_peaks
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh, expm_multiply
from scipy.special import ellipk

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model.carrier import federation

from ..electron_field.fan_recoil.fresh_kick.measure import divide, interval, phase
from ..row_20_electric.experiment import guide

HERE = Path(__file__).resolve().parent
SHOTS = 2**42
SCAN_SHOTS = 2**34
SEED = 30420608
LINKS = 8192
STEP = 1 / 16
LENGTHS = (4, 6, 8)
TICKS = 2**32
DRAW_PAD = 5e-9
SIGNED_PHASE_ENERGY_ERROR = 2e-8


def count(rng, values, shots=SHOTS):
    p = np.array(values, float)
    assert p.min() > -1e-12 and abs(p.sum() - 1) < 1e-12
    q = np.rint(np.clip(p, 0, 1) * TICKS).astype(np.int64)
    q[q.argmax()] += TICKS - q.sum()
    assert np.max(abs(q / TICKS - p)) < DRAW_PAD
    return rng.multinomial(shots, q / TICKS)


def xy(rng, z, gate=1):
    return [
        count(rng, [gate * (1 + value) / 2, gate * (1 - value) / 2, 1 - gate]).tolist()
        for value in (z.real, z.imag)
    ]


def physical_basis(n):
    width = 2 * n
    words = [
        sum(1 << i for i in locations) for locations in combinations(range(width), n)
    ]
    addresses = {w: i for i, w in enumerate(words)}
    bits = ((np.array(words)[:, None] >> np.arange(width)) & 1).astype(np.int8)
    return words, addresses, bits


def matrix(words, lookup, bits, degrees, reference, endpoint, cut=False):
    """All native internal blank contacts including violating dimer words."""
    width = len(degrees)
    ports = 2 * bits
    source = law.seam_term(8)
    same = source[law.SRC.index((0, 0)), law.SRC.index((0, 0))]
    unlike = source[law.SRC.index((0, 2)), law.SRC.index((0, 2))]
    move = source[law.SRC.index((0, 2)), law.SRC.index((2, 0))]
    assert (same, unlike, move) == (0, 32.5, -0.5)
    joined = np.column_stack(
        (
            2 * np.ones(len(bits), dtype=int),
            ports,
            endpoint * np.ones(len(bits), dtype=int),
        )
    )
    equal = joined[:, 1:] == joined[:, :-1]
    diagonal = (8 * (ports - reference) ** 2) @ (6 - degrees)
    permitted = np.ones(width + 1, dtype=bool)
    if cut:
        permitted[2:-1:2] = False
    diagonal = diagonal + unlike * (~equal[:, permitted]).sum(axis=1)
    rows = list(range(len(bits)))
    cols = rows.copy()
    entries = diagonal.astype(float).tolist()
    for j in range(width - 1):
        if permitted[j + 1]:
            for col in np.flatnonzero(bits[:, j] != bits[:, j + 1]):
                rows.append(lookup[words[col] ^ (3 << j)])
                cols.append(int(col))
                entries.append(move)
    h = coo_matrix((entries, (rows, cols)), shape=(len(bits), len(bits))).tocsr()
    assert (h - h.T).nnz == 0
    return h, equal, diagonal


def eigen(h, k, tag):
    if h.shape[0] < 1000:
        values, v = eigh(h.toarray(), subset_by_index=(h.shape[0] - k, h.shape[0] - 1))
    else:
        v0 = np.random.default_rng(tag).standard_normal(h.shape[0])
        values, v = eigsh(h, k=k, which="LA", v0=v0, tol=2e-14)
    e = values[::-1]
    v = v[:, ::-1]
    assert np.linalg.norm(h @ v - v * e, ord="fro") < 1e-8
    assert np.linalg.norm(v.T @ v - np.eye(k), ord="fro") < 1e-9
    return e, v


def register(gap, rng, gate=1, first_power=-4):
    """Fit one signed frequency without using a furnished line wavelength."""
    previous = None
    banks = []
    for t in (2**b for b in range(first_power, 5)):
        bank = xy(rng, np.exp(1j * t * gap), gate)
        band = phase(*bank, 4 * DRAW_PAD / gate)
        n = 0 if previous is None else round((sum(previous) * t - sum(band)) / (4 * pi))
        this = [(q + n * 2 * pi) / t for q in band]
        if previous is not None:
            assert (previous[1] - previous[0]) * t < 2 * pi
            this = [max(this[0], previous[0]), min(this[1], previous[1])]
        assert this[0] < this[1]
        previous = this
        banks.append({"t": t, "counts": bank})
    return {
        "count_banks": banks,
        "phase_frequency_interval": [
            previous[0] - SIGNED_PHASE_ENERGY_ERROR,
            previous[1] + SIGNED_PHASE_ENERGY_ERROR,
        ],
    }


def calibrate(rng):
    clock = register(1, rng)
    scale = clock["phase_frequency_interval"]
    count_transfer = count(rng, [np.sin(1 / 2) ** 2, np.cos(1 / 2) ** 2])
    hop = [asin(sqrt(v)) for v in interval(int(count_transfer[0]), SHOTS, DRAW_PAD)]
    hop = divide(hop, scale)
    unlike = xy(rng, np.exp(1j * 32.5 / 32))
    cost = divide([32 * q for q in phase(*unlike, 4 * DRAW_PAD)], scale)
    return {
        "clock": clock,
        "isolated_transfer_counts": count_transfer.tolist(),
        "exchange_interval": hop,
        "isolated_unlike_xy": unlike,
        "unlike_interval": cost,
    }


def uniform_model(h, w):
    """Antiferromagnetic uniform XXZ. Wall translation uses two port distances."""
    gamma = np.arccosh(w / h)
    modulus = brentq(
        lambda m: np.pi * ellipk(1 - m) / ellipk(m) - gamma, 1e-10, 1 - 1e-10
    )
    scale = 2 * h * ellipk(modulus) * np.sinh(gamma) / np.pi
    return float(scale * sqrt(1 - modulus)), float(scale**2 * modulus / 4)


def scout(rng, gap, spinon_weight, ref_admission):
    # Randomized sharp one-wall source. Preserve spectral-window refusals.
    weights = spinon_weight / len(spinon_weight) * ref_admission
    p = float(weights.sum())
    z = np.exp(1j * np.outer(STEP * np.arange(LINKS), gap)) @ weights
    bins = np.array(
        [count(rng, [(p + q) / 2, (p - q) / 2, 1 - p], SCAN_SHOTS) for q in z.real]
    )
    read = (bins[:, 0] - bins[:, 1]) / SCAN_SHOTS
    fft = np.abs(np.fft.rfft(np.hanning(LINKS) * read, 2**18))
    freq = 2 * pi * np.fft.rfftfreq(2**18, STEP)
    peaks, _ = find_peaks(fft)
    peaks = peaks[(freq[peaks] > 31) & (freq[peaks] < 34) & (fft[peaks] > 20)]
    assert len(peaks) == len(gap)
    logs = np.log(fft[np.array([peaks - 1, peaks, peaks + 1])])
    correction = (logs[0] - logs[2]) / (2 * (logs[0] - 2 * logs[1] + logs[2]))
    measured = (peaks + correction) * freq[1]
    assert np.max(abs(measured - gap)) < 0.001
    return measured, bins


def wavelength(rng, amplitudes, gate):
    # Coarse unique-wall flag is taken without resolving which bond. Native
    # dimer quarter exchanges compare adjacent wall-position words within it.
    p = float(np.vdot(amplitudes, amplitudes).real)
    banks = []
    diffs = []
    for parity in (0, 1):
        a = amplitudes[parity:-1:2]
        b = 1j * amplitudes[parity + 1 :: 2]
        plus = float(np.vdot(b + 1j * a, b + 1j * a).real / 2)
        minus = float(np.vdot(a + 1j * b, a + 1j * b).real / 2)
        bank = count(
            rng,
            [
                gate * plus,
                gate * minus,
                gate * (p - plus - minus),
                gate * (1 - p),
                1 - gate,
            ],
        )
        retained = sum(bank[:3])
        # The exact band is isolated and each line has neighbours separated
        # by more than 0.10, giving a common residual-based frame allowance.
        guard = 5e-7 + 6 * DRAW_PAD / gate
        x = interval(int(bank[0]), int(retained), guard)
        y = interval(int(bank[1]), int(retained), guard)
        diffs.append([x[0] - y[1], x[1] - y[0]])
        banks.append(bank.tolist())
    c = [sum(p) / 2 for p in zip(*diffs, strict=True)]
    assert -1 < c[0] < c[1] < 1
    q = [pi - acos(x) for x in c]
    density = count(rng, [*(gate * abs(amplitudes) ** 2), gate * (1 - p), 1 - gate])
    held = density[:-2] / sum(density[:-2])
    template = np.sin((np.arange(len(amplitudes)) + 1) * sum(q) / 2) ** 2
    template /= sum(template)
    return {
        "counts_right_left_unpaired_non_single_wall_upstream_refusals": banks,
        "measured_wave_number_from_band_minimum": q,
        "unique_wall_given_line": p,
        "held_wall_density_counts_then_other_and_upstream_refused": density.tolist(),
        "maximum_held_wall_density_error": float(np.max(abs(held - template))),
    }


def flight(rng, h, cut, equal, wall_ids):
    n = len(wall_ids) - 1
    a = n // 2 - 1
    assert n == 8
    sites = np.flatnonzero(equal[wall_ids[a]])
    pole = 2 * (a + 3)  # fixed local comparator outside both initial words
    assert len(sites) == 1 and pole not in (2 * a, 2 * (a + 1))
    initial = np.zeros((h.shape[0], 2), complex)
    initial[wall_ids[a], :] = 1 / sqrt(2)
    initial[wall_ids[a + 1], :] = np.array([1j, -1j]) / sqrt(2)
    phase_cohorts = [
        xy(rng, 2 * initial[wall_ids[a], s].conj() * initial[wall_ids[a + 1], s])
        for s in range(2)
    ]
    output = {}
    for name, op in (("native", h), ("cut", cut)):
        # A common offset avoids long arbitrary phase from the static stock.
        op = op.copy()
        op.setdiag(op.diagonal() - op[wall_ids[a], wall_ids[a]])
        off = float(np.trace(initial.conj().T @ (op @ initial)).real)
        assert abs(off) < 1e-12
        evol = expm_multiply(
            -1j * op,
            initial,
            start=0,
            stop=4,
            num=5,
            traceA=-1j * float(op.diagonal().sum()),
        )
        rows = []
        for time in (0, 1, 2, 4):
            data = []
            for s in range(2):
                v = abs(evol[time, :, s]) ** 2
                retained = v[wall_ids]
                other = 1 - retained.sum()
                meter = float(v @ equal[:, pole])
                counts = count(rng, [*retained, other])
                pcount = count(rng, [meter, 1 - meter])
                center = float(counts[:-1] @ np.arange(n + 1) / counts[:-1].sum())
                data.append(
                    {
                        "wall_address_counts_then_other": counts.tolist(),
                        "remote_equal_bond_counts": pcount.tolist(),
                        "center_given_one_wall": center,
                        "ideal_remote_equal_bond": meter,
                    }
                )
            delta = data[0]["center_given_one_wall"] - data[1]["center_given_one_wall"]
            picked = [i["ideal_remote_equal_bond"] for i in data]
            if name == "cut":
                assert abs(picked[0] - picked[1]) < 2e-12
            remote_a, remote_b = (
                interval(z["remote_equal_bond_counts"][0], SHOTS, 1e-9 + DRAW_PAD)
                for z in data
            )
            rows.append(
                {
                    "wait": time,
                    "sign_arms": data,
                    "centroid_difference": delta,
                    "remote_sign_difference_interval": [
                        remote_a[0] - remote_b[1],
                        remote_a[1] - remote_b[0],
                    ],
                }
            )
        output[name] = rows
    assert abs(output["native"][-1]["centroid_difference"]) > 1
    return {
        "source_wall_addresses": [2 * a, 2 * (a + 1)],
        "remote_physical_bond": pole,
        "phase_labels": [pi / 2, -pi / 2],
        "independent_initial_wall_phase_counts": phase_cohorts,
        "flights": output,
    }


def endpoint_control(rows):
    """Fix each rival on its own endpoint records before forecasting interior lines."""

    def examine(coordinate, exponent, result):
        coordinates = [
            [
                coordinate(k)
                for k in row["spatial_read"]["measured_wave_number_from_band_minimum"]
            ]
            for row in rows
        ]
        energies = [
            [e**exponent for e in row["line_clock"]["frequency_clock_units"]]
            for row in rows
        ]
        x_left, x_right = coordinates[0], coordinates[-1]
        e_left, e_right = energies[0], energies[-1]
        for row, x, e in zip(
            rows[1:-1], coordinates[1:-1], energies[1:-1], strict=True
        ):
            weight = divide(
                [x[0] - x_left[1], x[1] - x_left[0]],
                [x_right[0] - x_left[1], x_right[1] - x_left[0]],
            )
            predicted = [
                e_left[0] + (e_right[0] - e_left[1]) * weight[0],
                e_left[1] + (e_right[1] - e_left[0]) * weight[1],
            ]
            row[result] = [e[0] - predicted[1], e[1] - predicted[0]]

    examine(lambda q: 1 - cos(q), 1, "endpoint_calibrated_cosine_foil_error")
    examine(
        lambda q: 2 - 2 * cos(q), 2, "endpoint_calibrated_lattice_squared_shell_error"
    )
    examine(lambda q: q * q, 2, "endpoint_calibrated_continuum_squared_shell_error")


def arm(rng, n, ports, degrees, reference, bench):
    words, ids, bits = physical_basis(n)
    h, equal, diag = matrix(words, ids, bits, degrees, reference, 2)
    host, eqhost, _ = matrix(words, ids, bits, degrees, reference, 0)
    walls = equal.sum(axis=1)
    assert np.all(walls % 2 == 1)
    wids = np.flatnonzero(walls == 1)
    wids = wids[np.argsort(equal[wids].argmax(axis=1))]
    np.testing.assert_array_equal(equal[wids].argmax(axis=1), np.arange(n + 1) * 2)
    assert len(wids) == n + 1 and np.ptp(diag[wids]) == 0
    control = h[wids][:, wids].toarray()
    np.testing.assert_array_equal(
        control,
        np.diag(np.full(n + 1, diag[wids][0]))
        + np.diag(np.full(n, -0.5), 1)
        + np.diag(np.full(n, -0.5), -1),
    )
    other = np.setdiff1d(np.arange(len(bits)), wids)
    q = h[other][:, other]
    upper = float(
        np.max(
            q.diagonal() + np.asarray(abs(q).sum(axis=1)).ravel() - abs(q.diagonal())
        )
    )
    escape = h[other][:, wids].getnnz(axis=0)
    assert min(escape) > 0
    ep, vec = eigen(h, n + 1, n + reference)
    eh, vh = eigen(host, 2, 100 + n + reference)
    neu = np.flatnonzero(eqhost.sum(axis=1) == 0)
    assert len(neu) == 1 and eh[0] - eh[1] > 10
    ref_accepted = round(float(abs(vh[neu[0], 0]) ** 2), 11)
    host_read = register(round(float(eh[0]), 9), rng, ref_accepted, first_power=-10)
    assert (
        host_read["phase_frequency_interval"][0]
        < eh[0]
        < host_read["phase_frequency_interval"][1]
    )
    gaps = np.round(eh[0] - ep, 9)
    assert np.min(gaps) > 31 and np.max(gaps) < 34 and upper < eh[0] - 34
    assert min(np.diff(gaps)) > 0.1
    wall_frame = np.round(vec[wids], 11)
    weights = (abs(wall_frame) ** 2).sum(axis=0)
    line_window, bins = scout(rng, gaps, weights, ref_accepted)
    clock = bench["clock"]["phase_frequency_interval"]
    rows = []
    for j, picked in enumerate(line_window):
        hit = np.flatnonzero(abs(gaps - picked) < 0.004)
        assert len(hit) == 1 and hit[0] == j
        gate = float(weights[j] * ref_accepted / (n + 1))
        timed = register(gaps[j], rng, gate)
        timed["frequency_clock_units"] = divide(
            timed["phase_frequency_interval"], clock
        )
        spatial = wavelength(rng, wall_frame[:, j], gate)
        qlo, qhi = spatial["measured_wave_number_from_band_minimum"]
        clo, chi = timed["frequency_clock_units"]
        rlo, rhi = bench["xxz_rest_mass_interval"]
        vlo, vhi = bench["xxz_velocity_second_order_interval"]
        lattice = [
            sqrt(rlo**2 + 4 * vlo * np.sin(qlo / 2) ** 2),
            sqrt(rhi**2 + 4 * vhi * np.sin(qhi / 2) ** 2),
        ]
        continuum = [sqrt(rlo**2 + vlo * qlo**2), sqrt(rhi**2 + vhi * qhi**2)]
        rows.append(
            {
                "scouted_frequency": float(picked),
                "spatial_read": spatial,
                "line_clock": timed,
                "uniform_xxz_lattice_error": [clo - lattice[1], chi - lattice[0]],
                "uniform_xxz_continuum_error": [clo - continuum[1], chi - continuum[0]],
            }
        )
    out = {
        "dimers": n,
        "reference_load": reference,
        "physical_ports": ports,
        "physical_degrees": degrees.tolist(),
        "configurations": len(bits),
        "one_wall_addresses": len(wids),
        "nonprojected_word_couplings_per_one_wall": escape.tolist(),
        "off_one_wall_block_upper": upper,
        "host_reference_gate": ref_accepted,
        "independent_host_phase_read": host_read,
        "lines": rows,
    }
    endpoint_control(rows)
    if n == 8 and reference == 0:
        cut, _, _ = matrix(words, ids, bits, degrees, reference, 2, cut=True)
        out["direct_phase_reversal"] = flight(rng, h, cut, equal, wids)
    return out, bins


def run():
    rng = np.random.default_rng(SEED)
    bench = calibrate(rng)
    h, w = bench["exchange_interval"], bench["unlike_interval"]
    nominal = uniform_model(sum(h) / 2, sum(w) / 2)
    # The corner-range supplement is a development guard, rather than a
    # derivative certificate for the elliptic integral under all settings.
    corners = [uniform_model(a, b) for a in h for b in w]
    bench["xxz_rest_mass_interval"] = [
        min(c[0] for c in corners) - 0.0001,
        max(c[0] for c in corners) + 0.0001,
    ]
    bench["xxz_velocity_second_order_interval"] = [
        min(c[1] for c in corners) - 0.0001,
        max(c[1] for c in corners) + 0.0001,
    ]
    bench["xxz_prediction_center"] = list(nominal)
    board = federation(2)
    _, source, _ = guide(board)
    arrays = {}
    results = []
    for n, reference in (*[(k, 0) for k in LENGTHS], (8, 1)):
        ports = source[1 : 2 * n + 1]
        deg = np.asarray(board.degree)[ports]
        np.testing.assert_array_equal(deg, np.tile([6, 6, 5, 5], n // 2))
        assert all(
            b in board.neighbours[a] for a, b in zip(ports[:-1], ports[1:], strict=True)
        )
        row, raw = arm(rng, n, ports, deg, reference, bench)
        row["scout_count_array"] = "wall_" + str(n) + "_" + str(reference)
        arrays[row["scout_count_array"]] = raw
        results.append(row)
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "bench": bench,
        "bank_originals": SHOTS,
        "detector_probability_ticks": TICKS,
        "detector_quantization_error_allowance": DRAW_PAD,
        "scout_originals": SCAN_SHOTS,
        "scout_wait_spacing": STEP,
        "scout_bins": ["phase_plus", "phase_minus", "host_or_wave_window_refused"],
        "runs": results,
    }, arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result, arrays = run()
    result = json.loads(json.dumps(result, allow_nan=False))
    path = HERE / "records.json"
    raw = HERE / "scout_counts.npz"
    if args.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        np.savez_compressed(raw, **arrays)
    else:
        from ..electron_field.replay import compare

        compare(result, json.loads(path.read_text()), "native_wall")
        with np.load(raw) as receipt:
            assert set(receipt.files) == set(arrays)
            for name, array in arrays.items():
                np.testing.assert_array_equal(array, receipt[name])
        print("Verified unguarded electric material wall counts and spectra.")


if __name__ == "__main__":
    main()
