"""Energy-resolved free emission into two initially unit-filled owners."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from model.experiments.interval_echo_gauge.flag_measure.measure import meter

from .bound import certificate, receiver_envelope
from .apparatus import frame_columns, geometry_record, local_frames, read_local, sector

HERE = Path(__file__).resolve().parent
SHOTS = 2**22
TIMES = (0, 8, 32, 64, 128, 192, 256, 512)
STRENGTH = 8


def samples(p, rng):
    p = np.asarray(p)
    if abs(float(p.sum()) - 1) > 1e-9:
        raise ValueError("incomplete observer outcome law")
    # Quantize probabilities below one part in 10^12 before seeded sampling.
    bins = np.rint(p.ravel() * 10**12).astype(np.int64)
    bins[bins.argmax()] += 10**12 - int(bins.sum())
    return rng.multinomial(SHOTS, bins / 10**12).reshape(p.shape)


def spectra():
    blocks, frames = local_frames(STRENGTH)
    levels = [np.diag(f.T @ h @ f) for h, f in zip(blocks, frames, strict=True)]
    return frames, levels


def calibrate(frames, rng):
    """Fresh regional calibrations with actual load and exchange-parity outcomes."""
    measured = []
    for kind, f in zip(("vacancy", "neutral_pair"), frames, strict=True):
        arms = []
        for src in f.T:
            # Terminal swap parity, obtained by the private-phase, half-swap decoder.
            pminus = float(abs(src[1] - src[0]) ** 2 / 2)
            parity = samples([1 - pminus, pminus], rng)
            position = samples(abs(src) ** 2, rng)
            # Count squared contrasts on every physical active bond. A boundary
            # bond reaching the filled buffer contributes half its unequal-load flag.
            tags = np.asarray(meter())
            rails = 8 if kind == "vacancy" else 6
            inner = 0 if kind == "vacancy" else 3
            orders = ((0, 1), (1, 0)) if kind == "vacancy" else ((0, 2), (2, 0))
            diagonal = []
            for order in orders:
                word = np.ones(rails, dtype=int)
                word[inner : inner + 2] = order
                differences = np.diff(tags[word])
                boundary = np.delete(differences, inner)
                diagonal.append(
                    STRENGTH * int(differences @ differences)
                    + np.count_nonzero(boundary) / 2
                )
            e = (position @ diagonal + parity[1]) / SHOTS
            arms.append(
                dict(
                    position=position.tolist(),
                    swap_plus_minus=parity.tolist(),
                    measured_energy=round(e, 9),
                )
            )
        measured.append(arms)
    return measured


def factor(n, strength, cut, which, frames, bits):
    charged = n == 8
    words, h = sector(n, 2 if charged else 1, 1, strength, cut)
    if charged:
        initial, _ = frame_columns(words, [(0, 1), (3, 4)], frames)
        final, _ = frame_columns(words, [(0, 1), (5, 6)], frames)
        column = 2 * bits[0] + bits[1]
        readers = [(0, 1, 2, (0, 1), frames[0]), (5, 6, 7, (0, 2), frames[1])]
    else:
        initial, _ = frame_columns(words, [(1, 2)], [frames[1]])
        final, _ = frame_columns(words, [(3, 4)], [frames[1]])
        column = bits[-1]
        readers = [(3, 4, 5, (0, 2), frames[1])]
    incoming = initial[:, column]
    target = final[:, column]
    if which == "reverse":
        incoming = target
        if charged:
            readers[-1] = (3, 4, 2, (0, 2), frames[1])
        else:
            readers[0] = (1, 2, 0, (0, 2), frames[1])
    ev, v = eigh(h)
    coef = v.T @ incoming
    by_time = []
    for t in TIMES:
        y = v @ (np.exp(-1j * ev * t) * coef)
        by_time.append(read_local(y, words, readers))
    hp = initial.T @ h @ initial
    hd = final.T @ h @ final

    def moments(x):
        z = h @ x
        mean = float(x @ z)
        return [
            round(mean, 11),
            round(float(np.vdot(z - mean * x, z - mean * x).real), 11),
        ]

    result = dict(
        dimension=len(words),
        source_energy_and_variance=moments(initial[:, column]),
        target_energy_and_variance=moments(target),
        independent_chamber_energy_difference=round(float(abs(hp - hd).max()), 11),
    )
    return by_time, result


def run_case(name, strength, frames, bits, rng, cut=False, reverse=False):
    which = "reverse" if reverse else "forward"
    a, ma = factor(8, strength, 4 if cut else None, which, frames, bits)
    b, mb = factor(6, strength, None, which, frames, bits)
    rows = []
    for t, pa, pb in zip(TIMES, a, b, strict=True):
        joint = pa[:, :, None] * pb[None, None, :]
        counts = samples(joint, rng)
        rows.append(
            dict(
                exposure=t,
                three_read_counts=counts.tolist(),
                matching_count=int(counts[bits]),
                matching_probability=round(float(joint[bits]), 11),
                all_three_complete=int(counts[:2, :2, :2].sum()),
            )
        )
    return dict(
        setting=name,
        strength=strength,
        source_levels=list(bits),
        charged_factor=ma,
        second_factor=mb,
        records=rows,
    )


def without_surplus(frames, rng):
    """Same neutral pair and route after filling the single surplus vacancy."""
    w, h = sector(8, 1, 1, STRENGTH)
    s, _ = frame_columns(w, [(3, 4)], [frames[1]])
    ev, v = eigh(h)
    y = v @ (np.exp(-128j * ev) * (v.T @ s[:, 1]))
    p = read_local(y, w, [(5, 6, 7, (0, 2), frames[1])])
    return dict(
        neutral_ground_excited_miss=samples(p, rng).tolist(),
        probabilities=np.round(p, 11).tolist(),
    )


def no_pairs(frames, rng):
    w, h = sector(8, 1, 0, STRENGTH)
    s = np.zeros(len(w))
    for j in (0, 1):
        word = [1] * 8
        word[j] = 0
        s[w.index(tuple(word))] = frames[0][j, 0]
    ev, v = eigh(h)
    y = v @ (np.exp(-128j * ev) * (v.T @ s))
    pa = read_local(y, w, [(0, 1, 2, (0, 1), frames[0]), (5, 6, 7, (0, 2), frames[1])])
    assert np.all(pa[:, :2] == 0)
    pb = np.array([0, 0, 1])
    return samples(pa[:, :, None] * pb, rng).tolist()


def all_receiver_words(frames, rng):
    """Lower charge read joined to every nontrivial zero-charge word in owner 6."""
    words, h = sector(8, 2, 1, STRENGTH)
    s = frame_columns(words, [(0, 1), (3, 4)], frames)[0][:, 3]
    ev, v = eigh(h)
    cols = []
    for vacancy in range(5, 8):
        for e in range(5, 8):
            if e == vacancy:
                continue
            column = np.zeros(len(words))
            for c in (0, 1):
                w = [1] * 8
                w[c] = w[vacancy] = 0
                w[e] = 2
                column[words.index(tuple(w))] = frames[0][c, 0]
            cols.append(column)
    records = []
    for t in TIMES:
        y = v @ (np.exp(-1j * ev * t) * (v.T @ s))
        p = abs(np.asarray(cols) @ y) ** 2
        records.append(
            dict(
                exposure=t,
                six_words_and_miss=samples([*p, 1 - p.sum()], rng).tolist(),
                emission_probability=round(float(p.sum()), 12),
            )
        )
    return records


def run():
    frames, levels = spectra()
    rng = np.random.default_rng(8297314)
    settings = [
        run_case("two_odd", 8, frames, (0, 1, 1), rng),
        run_case("two_even", 8, frames, (0, 0, 0), rng),
        run_case("raised_charge", 8, frames, (1, 1, 1), rng),
        run_case("reversed", 8, frames, (0, 1, 1), rng, reverse=True),
        run_case("first_bridge_cut", 8, frames, (0, 1, 1), rng, cut=True),
        run_case("half_contrast_during_flight", 4, frames, (0, 1, 1), rng),
    ]
    return dict(
        shots_per_setting=SHOTS,
        geometry=geometry_record(),
        native_sources_and_readers=[x.round(10).tolist() for x in frames],
        daughter_levels=[x.round(11).tolist() for x in levels],
        charged_and_first_neutral_calibration=calibrate(frames, rng),
        second_neutral_calibration=calibrate(frames, rng)[1],
        settings=settings,
        no_seeded_pair=no_pairs(frames, rng),
        first_wing_filled_surplus_control=without_surplus(frames, rng),
        raised_charge_drop_bound=certificate(),
        complete_neutral_receiver_bound=receiver_envelope(),
        raised_charge_all_receiver_words=all_receiver_words(frames, rng),
        outcomes="charged, first neutral, second neutral; each index is lower, upper, or miss",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    else:
        saved = json.loads((HERE / "results.json").read_text())
        assert result == saved
        print("Free-emission records and all three local calibrations reproduced.")


if __name__ == "__main__":
    main()
