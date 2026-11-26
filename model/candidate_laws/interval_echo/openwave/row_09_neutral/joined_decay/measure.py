"""Read signed ordinary-energy transfers from one distant dimer on an intact mount."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from ..passive_decay.measure import SHOTS, WAITS, hamiltonians, hardware, physical

N_CLOCK = 16384
PERIOD = 128
SEED = 391107
HERE = Path(__file__).resolve().parent


def core_pulse(frame):
    """Exchange the outer isolated-core levels with an explicit real phase convention."""
    columns = frame[:, [0, -1]].copy()
    markers = np.argmax(abs(columns), axis=0)
    columns *= np.sign(columns[markers, np.arange(2)])
    difference = columns[:, 0] - columns[:, 1]
    return np.eye(len(frame)) - np.outer(difference, difference)


def delay_kernel(delta, bins, labels=N_CLOCK):
    """Coherent sum over the time at which the distant Z pulse is applied."""
    phase = bins[:, None, None] / labels - delta / PERIOD
    phase = (phase + 0.5) % 1.0 - 0.5
    return (
        np.exp(1j * np.pi * (labels - 1) * phase)
        * np.sinc(labels * phase)
        / np.sinc(phase)
    )


def instruments(energies, probe, cutoff, bins):
    delta = energies[None, :] - energies[:, None]
    assert np.ptp(energies) < PERIOD / 2
    effects = []
    for band in (bins, -bins):
        paths = delay_kernel(delta, band) * probe
        maps = paths.reshape(-1, len(energies))
        effects.append(maps.conj().T @ maps)
    excluded = probe * ~((delta > 0) & (delta < cutoff))
    null_maps = (delay_kernel(delta, bins) * excluded).reshape(-1, len(energies))
    null_bound = np.linalg.eigvalsh(null_maps.conj().T @ null_maps)[-1]
    effects.append(np.eye(len(energies)) - sum(effects))
    assert min(np.linalg.eigvalsh(effect)[0] for effect in effects) > -1e-12
    return effects, float(null_bound)


def draw(rng, probabilities):
    p = np.asarray(probabilities)
    assert abs(p.sum() - 1) < 1e-10 and p.min() > -1e-10
    p = np.maximum(0, p)
    return rng.multinomial(SHOTS, p / p.sum()).tolist()


def run():
    board, dimers, heads = hardware()
    assert len(dimers) == 4
    source, _, split, attachment = hamiltonians(4)
    actual = physical(board, dimers, heads)
    np.testing.assert_allclose(
        actual - actual[0], np.diag(split + attachment) - split[0, 0]
    )
    source_levels, core = eigh(source)
    cutoff = (source_levels[-1] - source_levels[0]) / 10
    assert cutoff > 2
    local_pulse = np.kron(core_pulse(core), np.eye(4))
    far_z = np.diag(1 - 2 * (np.arange(16) % 2))
    far_x = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    read_ports = np.kron(core, np.kron(np.eye(2), far_x))
    bins = np.arange(N_CLOCK // (4 * PERIOD), 2 * N_CLOCK // PERIOD + 1)
    rng = np.random.default_rng(SEED)
    panels = []
    for name, h, pulsed in [
        ("intact pulse", split + attachment, True),
        ("intact idle", split + attachment, False),
        ("cut pulse", split, True),
    ]:
        energies, levels = eigh(h)
        ground = levels[:, 0]
        prepared = local_pulse @ ground if pulsed else ground
        amplitudes = levels.T @ prepared
        effects, bound = instruments(energies, levels.T @ far_z @ levels, cutoff, bins)
        flight, absorption = [], []
        for wait in WAITS:
            traveling = np.exp(-1j * wait * energies) * amplitudes
            p = abs(read_ports.T @ (levels @ traveling)).reshape(4, 2, 2) ** 2
            flight.append(
                {
                    "wait": wait,
                    "core_by_far_x": np.reshape(
                        draw(rng, p.sum(axis=1).ravel()), (4, 2)
                    ).tolist(),
                }
            )
            if wait in (0, 16):
                probs = [
                    np.vdot(traveling, effect @ traveling).real for effect in effects
                ]
                absorption.append({"wait": wait, "loss_gain_other": draw(rng, probs)})
        panels.append(
            {
                "arm": name,
                "far_probe_light_free_bound": round(bound, 11),
                "pulse_energy_above_joined_ground": round(
                    float(abs(amplitudes) ** 2 @ (energies - energies[0])), 9
                ),
                "flight": flight,
                "absorption": absorption,
            }
        )
    return {
        "shots": SHOTS,
        "seed": SEED,
        "clock_labels": N_CLOCK,
        "energy_period": PERIOD,
        "dimers": [list(pair) for pair in dimers],
        "absorption_band": [0.25, 2.0],
        "light_cutoff": round(float(cutoff), 9),
        "panels": panels,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    dest = parser.parse_args().output
    record = run()
    if dest:
        dest.write_text(json.dumps(record, indent=2) + "\n")
    else:
        assert record == json.loads((HERE / "records.json").read_text())
        print(
            "Verified the intact, idle and cut count banks and false-absorption bounds."
        )
