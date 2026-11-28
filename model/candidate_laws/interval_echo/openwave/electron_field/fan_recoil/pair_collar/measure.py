"""Paired sign, formation and recoil records against the same filled reference."""

import argparse
from hashlib import sha256
import json
from math import acos, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal, expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from ...replay import compare
from ..fresh_kick.measure import LOG, difference, divide, interval, phase
from ..fresh_kick.propagation.measure import (
    SHOTS,
    counts,
    difference_of_pointer,
    prepare,
    symbol,
    translation,
)
from .fixture import (
    CELLS,
    SOURCE,
    assert_joint_columns,
    geometry,
    signed_generator,
    source_vectors,
)

HERE = Path(__file__).resolve().parent
SEED = 101010455001
ANGULAR_PAD = 1e-9
WAIT = 128


def xy(rng, z):
    a = np.array([z.real, z.imag])
    assert max(abs(a)) < 1 + 2e-13
    a = np.clip(a, -1, 1)
    return [counts(rng, [(1 + x) / 2, (1 - x) / 2]) for x in a]


def frequency(rng, signal):
    previous = [0, 128]
    banks = []
    for tick in range(-6, 12):
        t = 2**tick
        bank = xy(rng, signal(t))
        band = phase(*bank, ANGULAR_PAD)
        assert (previous[1] - previous[0]) * t < pi
        advance = round((sum(previous) * t - sum(band)) / (4 * pi))
        candidate = [(q + 2 * pi * advance) / t for q in band]
        previous = [max(previous[0], candidate[0]), min(previous[1], candidate[1])]
        assert previous[0] < previous[1]
        banks.append({"wait": t, "pointer_counts": bank})
    return {"phase_interval": previous, "banks": banks}


def source_test(rng, paths):
    h, u, blank, written = source_vectors()
    signs = (
        "left_vacancy",
        "left_excess",
        "neutral_sham",
        "coherent_even",
        "coherent_odd",
    )
    ins = [written[:, i] for i in range(3)]
    ins += [(written[:, 0] + s * written[:, 1]) / sqrt(2) for s in (1, -1)]
    states = [u @ v for v in ins]
    load_bins = [law.SRC.index(word) for word in ((0, 2), (2, 0), (1, 1))]
    reports = {}
    for name, v, z in zip(signs, ins, states, strict=True):
        observer = [float(abs(z[i::9]) @ abs(z[i::9])) for i in load_bins]
        assert sum(observer) > 1 - 1e-14
        # The exact source label is total-neutral even when q_L and q_R differ.
        after_free = expm(-1j * h) @ z
        inv = u @ after_free
        exact_no_erasure = float(
            abs(after_free[load_bins[-1] :: 9]) @ abs(after_free[load_bins[-1] :: 9])
        )
        if name != "neutral_sham":
            assert exact_no_erasure < 1e-24
        written_return = float(
            abs(inv[load_bins[-1] :: 9]) @ abs(inv[load_bins[-1] :: 9])
        )
        assert written_return > 1 - 2e-14
        reports[name] = {
            "negative_left_positive_left_or_neutral_counts": counts(rng, observer),
            "vacuum_to_pair_native_free_probability": float(
                np.linalg.norm(blank[:, :2].T @ expm(-1j * h) @ v) ** 2
            ),
            "free_neutral_probability_after_creation": exact_no_erasure,
            "addressed_echo_return_to_filled_probability": written_return,
            "isolated_echo_mean_work": float(
                np.vdot(z, h @ z).real - np.vdot(v, h @ v).real
            ),
        }
    for name in signs:
        assert reports[name]["vacuum_to_pair_native_free_probability"] < 1e-24
    even, odd = states[-2:]
    gaps = [float(np.vdot(z, h @ z).real) for z in (even, odd)]
    assert all(
        np.linalg.norm(h @ z - e * z) < 1e-13
        for z, e in zip((even, odd), gaps, strict=True)
    )
    reference = frequency(rng, lambda t: np.exp(1j * t * (gaps[1] - gaps[0])))
    independent_clock = reference["phase_interval"]
    for label, z, g in zip(signs[-2:], (even, odd), gaps, strict=True):
        read = frequency(rng, lambda t: np.exp(1j * t * g))
        read["calibrated_contact_work"] = divide(
            read["phase_interval"], independent_clock
        )
        reports[label]["contact_work_read"] = read
    gate = expm(-1j * pi * law.seam_term(8))
    table = {
        pair: law.SRC[int(np.argmax(abs(gate[:, i])))] for i, pair in enumerate(law.SRC)
    }
    for i, (a, b) in enumerate(law.SRC):
        assert table[a, b] == (b, a)
        np.testing.assert_allclose(
            gate[law.SRC.index((b, a)), i], 1, rtol=0, atol=5e-14
        )
    for sig in (-1, 1):
        material = {p: 1 for arm in paths for p in arm}
        for q, p in zip((sig, -sig), SOURCE, strict=True):
            material[p] = 1 + q
        for arm, q in zip(paths, (sig, -sig), strict=True):
            for a, b in zip(arm[:3], arm[1:4], strict=True):
                material[a], material[b] = table[material[a], material[b]]
            assert all(material[p] == 1 for p in arm[:3]) and material[arm[3]] == 1 + q
        assert sorted(v for v in material.values() if v != 1) == [0, 2]
    return {
        "collar_modes": reports,
        "clock": reference,
        "source_ports": list(SOURCE),
        "physical_staging_full_swap_time": pi,
    }


def standing(index, sign, reference):
    k = pi * index / CELLS
    energy, u, _ = symbol(k, 0.5, 16 * sign * (1 - reference))
    j = np.arange(CELLS)
    f = np.exp(1j * k * j[:, None]) * u + np.exp(-1j * k * (j[:, None] + 1)) * u[::-1]
    a = f.ravel()
    return prepare(a / np.linalg.norm(a)), energy + 48


def independent_k(rng, vector):
    data = translation(rng, vector, 1, 1e-9)
    endpoints = vector.reshape(CELLS, 3)[[0, -1]]
    pp = float(
        np.sum(
            abs(endpoints[:, 1]) ** 2 + abs(endpoints[:, 0] + endpoints[:, 2]) ** 2 / 2
        )
    )
    pm = float(np.sum(abs(endpoints[:, 0] - endpoints[:, 2]) ** 2 / 2))
    b = counts(rng, [pp, pm, 1 - pp - pm])
    ends = difference_of_pointer(b, 1e-9)
    cal = [
        p + q / 2 for p, q in zip(data["x"]["quadrature_interval"], ends, strict=True)
    ]
    assert cal[0] > -1 and cal[0] < 1 and cal[1] > -1
    return {
        "three_contact_counts": data,
        "end_parity_counts": b,
        "cosine_interval": cal,
        "wave_number_interval": [acos(min(1, cal[1])), acos(cal[0])],
    }


def enveloped(k, u):
    env = sqrt(2 / 32) * np.sin(np.arange(1, 32) * pi / 32)
    x = np.zeros((CELLS, 3), complex)
    x[8:39] = env[:, None] * np.exp(1j * k * np.arange(8, 39))[:, None] * u
    return prepare(x.ravel())


def mobile_packet(k, sign):
    _, u, _ = symbol(k, 0.5, 16 * sign)
    return enveloped(k, u)


def release(vector, h):
    e, v = eigh_tridiagonal(h.diagonal(), h.diagonal(1))
    return v @ (np.exp(-1j * WAIT * (e - e[0])) * (v.T @ vector))


def spatial_motion(rng, h, cut, sign, clock):
    result = []
    rsync = WAIT * np.array(clock)
    for direction in (1, -1):
        a = mobile_packet(direction * pi / 2, sign)
        wave = translation(rng, a, 1, 1e-9)
        if direction > 0:
            assert wave["wave_number_interval"][0] > 0
        else:
            assert wave["wave_number_interval"][1] < 0
        shots = {}
        for label, packet in (
            ("initial", a),
            ("released", release(a, h)),
            ("cut", release(a, cut)),
        ):
            p = (abs(packet.reshape(CELLS, 3)) ** 2).sum(axis=1)
            bank = counts(rng, p)
            accepted = sum(bank[:-1])
            center = float(np.arange(CELLS) @ np.array(bank[:-1]) / accepted)
            error = (CELLS - 1) * (sqrt(LOG / (2 * accepted)) + 1e-9)
            shots[label] = {
                "owner_counts_then_refusal": bank,
                "centroid": [center - error, center + error],
            }
        v_est = divide(
            difference(shots["released"]["centroid"], shots["initial"]["centroid"]),
            rsync,
        )
        cut_v = divide(
            difference(shots["cut"]["centroid"], shots["initial"]["centroid"]), rsync
        )
        assert v_est[0] * direction > 0 and v_est[1] * direction > 0
        assert cut_v[0] < 0 < cut_v[1]
        result.append(
            {
                "direction": direction,
                "phase_readout": wave,
                "count_banks": shots,
                "signed_velocity": v_est,
                "cut_velocity": cut_v,
            }
        )
    return result


def joint_signed_counts(rng, board, paths):
    # Orthogonal ordered sign sectors of the coherent source are compatible
    # with both terminal sign reports. Unread source parity need not be sampled.
    # The coherent source receives the same pulse sequence on both signs.
    # Its equal-slot input is a comparison across bands, separate from the
    # branch-selective Ramsey and group-motion protocols.
    neutral_packet = enveloped(pi / 2, np.ones(3) / sqrt(3))
    data = {}
    for arm, ports in enumerate(paths):
        for sign in (-1, 1):
            generator = signed_generator(board, ports[3:], sign, 0)
            detector = release(neutral_packet, generator)
            far = float(np.vdot(detector[75:], detector[75:]).real)
            data[arm, sign] = np.array([1 - far, far])
    p = np.r_[np.kron(data[0, -1], data[1, 1]), np.kron(data[0, 1], data[1, -1])] / 2
    bank = counts(rng, p)
    exact = np.asarray(p).reshape(2, 2, 2)
    draws = np.asarray(bank[:-1]).reshape(2, 2, 2)
    for arm, block in enumerate(draws):
        sums = block.sum()
        minus = block[1].sum() if arm == 0 else block[:, 1].sum()
        plus = block[:, 1].sum() if arm == 0 else block[1].sum()
        assert (
            interval(int(minus), int(sums), 1e-9)[1]
            < interval(int(plus), int(sums), 1e-9)[0]
        )
    return {
        "order": [
            "left_sign_minus_or_plus",
            "left_far_false_or_true",
            "right_far_false_or_true",
        ],
        "wait": WAIT,
        "far_owner_index_at_least": 25,
        "common_packet_initial_far": float(
            np.vdot(neutral_packet[75:], neutral_packet[75:]).real
        ),
        "joint_counts_then_refusal": bank,
        "exact_joint_probabilities": exact.tolist(),
    }


def one_arm(rng, board, ports, q, reference, clock):
    h = signed_generator(board, ports, q, reference)
    lines = {}
    for index in (0, 2, 8, 24, 40) if reference == 0 else (0, 2, 24):
        a, e = standing(index, q, reference)
        # This exact bulk-reflection standing vector is physically prepared.
        # The residual controls using the phase of a one-line wave in readout.
        measured = float(np.vdot(a, h @ a).real)
        assert (
            abs(e - measured) < 2e-12 and np.linalg.norm(h @ a - measured * a) < 3e-12
        )
        read = frequency(rng, lambda t: np.exp(1j * t * measured))
        band = divide(read["phase_interval"], clock)
        lines[str(index)] = {
            "wave": independent_k(rng, a),
            "formation": read,
            "formation_interval_neutral_clock": band,
        }
    e0 = lines["0"]["formation_interval_neutral_clock"]
    for key, line in lines.items():
        if key == "0":
            continue
        e = line["formation_interval_neutral_clock"]
        k = line["wave"]["wave_number_interval"]
        de = difference(e, e0)
        line["recoil_interval"] = de
        line["squared_shell_coefficient_interval"] = divide(
            [e[0] ** 2 - e0[1] ** 2, e[1] ** 2 - e0[0] ** 2], [x**2 for x in k]
        )
    return h, lines


def overlap(ab, cd):
    return max(ab[0], cd[0]) <= min(ab[1], cd[1])


def trial(rng, board, paths, left_charge, reference, clock):
    arm = []
    for q, ports in zip((left_charge, -left_charge), paths, strict=True):
        h, records = one_arm(rng, board, ports[3:], q, reference, clock)
        read = {"sign": q, "standing_records": records}
        if reference == 0:
            cut = signed_generator(board, ports[3:], q, reference, cut=True)
            read["motion"] = spatial_motion(rng, h, cut, q, clock)
        arm.append(read)
    # Read an independent two-arm Ramsey product, from actual staged signs.
    # The projector/detector comparing each leg gets its own counted cohort.
    pair = {}
    for index in (0, 24):
        energy = []
        for q, ports in zip((left_charge, -left_charge), paths, strict=True):
            selected, _ = standing(index, q, reference)
            h = signed_generator(board, ports[3:], q, reference)
            energy.append(float(np.vdot(selected, h @ selected).real))
        reading = frequency(rng, lambda t: np.prod(np.exp(1j * t * np.array(energy))))
        value = divide(reading["phase_interval"], clock)
        left, right = (
            a["standing_records"][str(index)]["formation_interval_neutral_clock"]
            for a in arm
        )
        expected = [a + b for a, b in zip(left, right, strict=True)]
        assert overlap(value, expected)
        pair[str(index)] = {
            "measured": value,
            "forecast_from_separate_arms": expected,
            "read": reading,
        }
    a, b = sorted(arm, key=lambda p: p["sign"])
    ratios = {}
    for key in a["standing_records"]:
        x, y = a["standing_records"][key], b["standing_records"][key]
        assert overlap(
            x["wave"]["wave_number_interval"], y["wave"]["wave_number_interval"]
        )
        if key == "0":
            assert overlap(
                x["formation_interval_neutral_clock"],
                y["formation_interval_neutral_clock"],
            ) == (reference == 1)
        else:
            assert overlap(x["recoil_interval"], y["recoil_interval"]) == (
                reference == 1
            )
            assert overlap(
                x["squared_shell_coefficient_interval"],
                y["squared_shell_coefficient_interval"],
            ) == (reference == 1)
            ratios[key] = divide(y["recoil_interval"], x["recoil_interval"])
    if reference == 0:
        ar, bp = (p["motion"][0]["signed_velocity"] for p in (a, b))
        assert bp[0] > ar[1] > 0
    result = {
        "reference_load": reference,
        "left_charge": left_charge,
        "arm_records": arm,
        "positive_to_negative_recoil_ratios": ratios,
        "pair_clock": pair,
    }
    return result


def run():
    rng = np.random.default_rng(SEED)
    board, paths, protected = geometry()
    source = source_test(rng, paths)
    clock = source["clock"]["phase_interval"]
    results = []
    for reference in (0, 1):
        assert_joint_columns(board, paths, reference)
        for left in (-1, 1):
            results.append(trial(rng, board, paths, left, reference, clock))
    for i in (0, 2):
        for a, b in zip(
            sorted(results[i]["arm_records"], key=lambda p: p["sign"]),
            sorted(results[i + 1]["arm_records"], key=lambda p: p["sign"]),
            strict=True,
        ):
            for key, x in a["standing_records"].items():
                assert overlap(
                    x["formation_interval_neutral_clock"],
                    b["standing_records"][key]["formation_interval_neutral_clock"],
                )
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "source": source,
        "paths_with_staging": paths,
        "first_successful_protected_prefix": protected,
        "runs": results,
        "coherent_source_joint_arrival": joint_signed_counts(rng, board, paths),
        "trials_per_bank": SHOTS,
        "rng_seed": SEED,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    d = json.loads(json.dumps(run(), allow_nan=False))
    p = HERE / "records.json"
    if args.write:
        p.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    else:
        compare(d, json.loads(p.read_text()), "pair_collar")
        print(
            "Verified common-collar creation, both formation signs, recoil, and pair-clock counts."
        )


if __name__ == "__main__":
    main()
