"""Count resolved material excitations and predict a distant quarter-kick response."""

import argparse
import json
from pathlib import Path

import numpy as np

from .. import measure as ring

HERE = Path(__file__).resolve().parent
SHOTS = 2**26
WAITS = (0, 4, 8, 16, 32, 64, 128)
SEED = 970427


def aperture(phase, labels):
    """Mean exp(i*j*phase) over j=0,...,labels-1, including the zero limit."""
    x = np.remainder(np.asarray(phase) / (2 * np.pi) + 0.5, 1.0) - 0.5
    return np.exp(1j * np.pi * (labels - 1) * x) * np.sinc(labels * x) / np.sinc(x)


def sharp_cold(block):
    trans, e, modes, _, _ = block
    sharp = np.asarray(trans.getrow(0).toarray()).ravel().conj()
    amp = (modes.conj().T @ sharp) * aperture(-e / 512, 2048) ** 3
    success = float(np.vdot(amp, amp).real)
    amp /= np.sqrt(success)
    return amp, success


def signed_flight(source, blocks):
    """Return the two complete word amplitudes after independent local kick signs."""
    nwords = len(source)
    flipped = source[np.arange(nwords) ^ 1]
    launches = np.stack(
        (source - 1j * flipped, source + 1j * flipped), axis=1
    ) / np.sqrt(2)
    np.testing.assert_allclose(np.linalg.norm(launches, axis=0), 1, atol=1e-12)
    waves = [np.zeros_like(launches) for _ in WAITS]
    for trans, e, modes, _, _ in blocks:
        coefficients = modes.conj().T @ (trans.conj().T @ launches)
        for i, wait in enumerate(WAITS):
            waves[i] += trans @ (
                modes @ (np.exp(-1j * wait * e)[:, None] * coefficients)
            )
    return launches, waves


def local_bank(waves, shift, success, rng):
    flipped = np.arange(len(waves)) ^ (1 << shift)
    means = np.sum(waves.conj() * waves[flipped], axis=0).real
    counts = []
    for prob in np.clip((1 + means) / 2, 0, 1):
        admitted = int(rng.binomial(SHOTS, success))
        positive = int(rng.binomial(admitted, prob))
        counts.append([SHOTS - admitted, positive, admitted - positive])
    return counts


def run_panel(panel, owners):
    n = panel["dimer_count"]
    kappa = panel["kappa"]
    rng = np.random.default_rng([SEED, n, int(kappa)])
    _, blocks = ring.spectrum(n, kappa)
    trans0, e0, modes0, _, _ = blocks[0]
    cold, success = sharp_cold(blocks[0])
    distance = float(np.linalg.norm(cold[1:]))
    cold_words = trans0 @ modes0
    source = cold_words @ cold
    flipped_frame = cold_words[np.arange(1 << n) ^ 1]
    probabilities, calibration = [], []
    for k, (trans, e, modes, _, line) in enumerate(blocks):
        flip = modes.conj().T @ (trans.conj().T @ flipped_frame)
        frequency = ring.read_clock(panel["clocks"][k]["plus_counts"])
        delta = e[:, None] - e0[None, :] - frequency
        clock = aperture(delta / 128, 8192)
        p = float(np.linalg.norm((clock * flip) @ cold) ** 2)
        total = float(np.linalg.norm(flip @ cold) ** 2)
        probabilities.extend([p, max(0.0, total - p)])
        leakage = max(abs(np.delete(clock[:, 0], line)) ** 2)
        calibration.append(
            dict(
                k=k,
                frequency=frequency,
                line_acceptance=float(abs(clock[line, 0]) ** 2),
                other_acceptance_ceiling=float(leakage),
            )
        )
    np.testing.assert_allclose(sum(probabilities), 1, atol=2e-12)
    original = np.array([1 - success, *(success * np.asarray(probabilities))])
    counts = rng.multinomial(SHOTS, original / original.sum())
    launches, waves = signed_flight(source, blocks)
    words = np.arange(1 << n)
    bits = (words[:, None] >> np.arange(n)) & 1
    diag = kappa * (
        bits.sum(axis=1) + 2 * (bits != np.roll(bits, 1, axis=1)).sum(axis=1)
    )
    distal = []
    for wait, wave in zip(WAITS, waves, strict=True):
        disabled = np.exp(-1j * wait * diag)[:, None] * launches
        for shift in (2, n // 2):
            assert owners[shift] != owners[0]
            distal.append(
                dict(
                    wait=wait,
                    address=shift,
                    owner=owners[shift],
                    live=local_bank(wave, shift, success, rng),
                    exchange_off=local_bank(disabled, shift, success, rng),
                )
            )
    return dict(
        n=n,
        kappa=kappa,
        origin_owner=owners[0],
        cold_success=success,
        cold_trace_distance=distance,
        clock=calibration,
        refused=int(counts[0]),
        pole_and_complement=counts[1:].reshape(n, 2).tolist(),
        local=distal,
    )


def read(path):
    return json.loads(path.read_text())


def replay():
    prior = read(HERE.parent / "records.json")
    owners = {
        len(b["stock"]["owners"]): b["stock"]["owners"]
        for b in prior["boards"]
        if b["cells"] == 80 and b["stock"]
    }
    return dict(
        seed=SEED,
        starts=SHOTS,
        panels=[
            run_panel(p, owners[p["dimer_count"]])
            for p in prior["panels"]
            if p["cells"] == 80
        ],
    )


def decode(panel, shots=SHOTS):
    admitted = shots - panel["refused"]
    epsilon = np.sqrt(np.log(2048 / 0.01) / (2 * admitted))
    # The earlier full-spectrum enclosures place every excited level above 17
    # on the sharp-word energy scale. Bound the cold leakage for ANY input.
    pass_lower = admitted / shots - np.sqrt(np.log(2048 / 0.01) / (2 * shots))
    cold_bound = (1 / (2048 * np.sin(17 / 1024))) ** 3 / np.sqrt(pass_lower)
    error = epsilon + cold_bound
    contrast = np.asarray(panel["pole_and_complement"])[:, 0] / admitted
    leak = np.array([r["other_acceptance_ceiling"] for r in panel["clock"]])
    accept = np.array([r["line_acceptance"] for r in panel["clock"]])
    # Each ordered energy error is below 8e-9. For the squared filter the
    # global Lipschitz bound is (labels - 1) * step, less than 64.
    filter_error = 64 * 2 * 8e-9 + 1e-10
    leak += filter_error
    accept -= filter_error
    omega = np.array([r["frequency"] for r in panel["clock"]])
    lo = np.maximum(0.0, contrast - error - leak)
    hi = np.minimum(1.0, (contrast + error) / accept)
    discrepancies, flat_errors, signals, disabled = [], [], [], []
    for row in panel["local"]:

        def difference(name):
            a, b = np.array(row[name])
            return a[1] / sum(a[1:]) - b[1] / sum(b[1:])

        delay, site = row["wait"], row["address"]
        spatial = 2 * np.pi * np.arange(panel["n"]) * site / panel["n"]
        predicted = np.sin(spatial - omega * delay) @ contrast
        flat = np.sin(spatial - omega[0] * delay) @ contrast
        signal = difference("live")
        discrepancies.append(abs(signal - predicted))
        flat_errors.append(abs(signal - flat))
        signals.append(abs(signal))
        disabled.append(abs(difference("exchange_off")))
    return dict(
        min_scaled_weight_lower=float(panel["n"] * lo.min()),
        max_scaled_weight_upper=float(panel["n"] * hi.max()),
        summed_weight_lower=float(lo.sum()),
        max_pole_error=max(discrepancies),
        max_flat_error=max(flat_errors),
        max_local_response=max(signals),
        max_exchange_off_response=max(disabled),
    )


def rounded(value):
    if isinstance(value, dict):
        return {k: rounded(v) for k, v in value.items()}
    if isinstance(value, list):
        return [rounded(v) for v in value]
    return round(value, 11) if isinstance(value, float) else value


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    dest = parser.parse_args().output
    data = rounded(replay())
    if dest:
        dest.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        assert data == read(HERE / "records.json")
    for p in data["panels"]:
        print(p["n"], p["kappa"], decode(p))
