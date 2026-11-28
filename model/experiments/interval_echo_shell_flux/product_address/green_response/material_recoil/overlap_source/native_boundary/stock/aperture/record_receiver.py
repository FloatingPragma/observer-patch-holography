"""Use source counts and ordinary calibrations to predict a sparse repair intervention."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import expm

from . import pair_probe as probe
from . import record_density as source

HERE = source.HERE
contact = source.contact
TRIALS = 2**48


def fringe_calibration(clocks):
    phase = np.mean(
        [complex(*(2 * np.asarray(row) / contact.SHOTS - 1)) for row in clocks.values()]
    )
    error = 2 * np.sqrt(contact.RISK / contact.SHOTS)
    distance = abs(1 - phase)
    return np.array([max(0, distance - error) ** 2 / 2, (distance + error) ** 2 / 2])


def sharp_bright(dwell=contact.TAU):
    b, d = contact.blocks()["full"][2]
    z = expm(1j * dwell * b / 4) @ expm(-1j * dwell * (b + d) / 4)
    eyes = np.eye(len(b))
    # Every accepted sharp word has the same conditional visibility within its magnitude.
    np.testing.assert_allclose(z.diagonal(), np.full(len(b), z[0, 0]))
    paths = (eyes, z)
    return np.array(
        [(1 + np.kron(a, h.conj().T)[0, 0].real) / 2 for a in paths for h in paths]
    )


def release_calibration(records, gain):
    out = {}
    for cut, path in records["64"].items():
        if not cut.endswith("_undisturbed"):
            continue
        mixed = {}
        for arm, (norm, banks) in path.items():
            cells = {
                tag: probe.stats(banks[tag], norm) for tag in ("00", "01", "10", "11")
            }
            mixed[arm] = (
                cells["00"] + cells["11"] - cells["01"][::-1] - cells["10"][::-1]
            )
        out[cut] = mixed["kinetic"] + probe.mix_bands(gain, mixed["contrast"])
    return out


def run(bank, clocks, records, tape):
    calibrated = fringe_calibration(clocks)
    cuts = release_calibration(records, contact.timing.lab.ordinary_gain(tape) / 4)
    _, classes = source.classes()
    panel = [(len(w), source.counts(w)) for w in classes]
    rng = np.random.default_rng(982451653)
    bright = sharp_bright()
    error = np.sqrt(2 * contact.RISK / TRIALS)
    cases = {}
    for menu, row in bank["cases"]["1/16"].items():
        counts = row["counts"]
        odds = np.array(counts["reference_log_odds"])
        assert odds[0] > 0
        # An observer supplies only a signed modular charge and log-odds interval.
        delta = probe.mix_bands(-np.array(counts["source_nats_per_attempt"]), 1 / odds)
        forecast = probe.mix_bands(calibrated, delta)
        old, new, _ = source.recover(
            source.combine(panel, source.input_weights(source.LEVELS[0]), menu)
        )

        def receiver(accepted):
            # Unaccepted attempts carry a rejection flag, hence zero bright count.
            return int(
                rng.binomial(TRIALS, float(np.asarray(accepted, float) @ bright))
            )

        counts_bright = {
            k: receiver(v)
            for k, v in (("native", old), ("omitted", old), ("repaired", new))
        }

        def difference(k):
            return (counts_bright[k] - counts_bright["native"]) / TRIALS + np.array(
                [-error, error]
            )

        passive = difference("omitted")
        held = difference("repaired")
        true = float(
            np.array([float(b - a) for a, b in zip(old, new, strict=True)]) @ bright
        )
        assert forecast[0] <= true <= forecast[1] and held[0] <= true <= held[1]
        assert passive[0] <= 0 <= passive[1]
        cases[menu] = dict(
            delta_from_charge=delta.tolist(),
            bright=counts_bright,
            passive=passive.tolist(),
            from_charge_fringe=forecast.tolist(),
            held_fringe=held.tolist(),
            from_charge_release={
                k: probe.mix_bands(delta, v).tolist() for k, v in cuts.items()
            },
        )
    sparse = bank["cases"]["1/16"]["none"]
    dense = bank["cases"]["1/4"]["none"]
    # A source-controlled shutter mixes native and repaired copies at unchanged reference.
    thin = sparse["receiver_nats_per_attempt"] / dense["receiver_nats_per_attempt"]
    matched = thin * float(source.Q(dense["delta_both"]))
    assert 0 < thin < 1 and matched > 0 > float(source.Q(sparse["delta_both"]))
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "RECORD_DENSITY.md",
                HERE / "record_density.json",
                HERE / "pair_contact.json",
                HERE / "pair_probe.json",
                HERE / "mediated_release.json",
            )
        },
        trials=TRIALS,
        fringe_calibration=calibrated.tolist(),
        release_calibration={k: v.tolist() for k, v in cuts.items()},
        cases=cases,
        matched_charge=dict(
            nats=sparse["receiver_nats_per_attempt"],
            quarter_fill_repair_probability=thin,
            quarter_fill_delta=matched,
            sixteenth_fill_delta=float(source.Q(sparse["delta_both"])),
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()

    def read(name):
        return json.loads((HERE / (name + ".json")).read_bytes())

    out = run(
        read("record_density"),
        read("pair_contact")["clock"],
        read("pair_probe")["ordinary_banks"],
        read("mediated_release"),
    )
    path = HERE / "record_receiver.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        contact.timing.lab.lab.cold.apparatus.check(read("record_receiver"), out)
    print(out["matched_charge"])
    for menu, row in out["cases"].items():
        print(
            menu,
            row["from_charge_fringe"],
            row["held_fringe"],
            row["from_charge_release"],
        )
