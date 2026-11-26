"""Test transmission and signed light absorption four active links from a heavy source."""

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from ..joined_decay.measure import core_pulse, delay_kernel, draw, N_CLOCK
from ..passive_decay.measure import (
    SHOTS,
    WAITS,
    physical,
    compare_records,
    hamiltonians,
)
from ..passive_decay.material import columns, edge
from .material import route

HERE = Path(__file__).resolve().parent


def rebuild(board, dimers, heads, cut, crossing=False):
    n = len(dimers)
    primary, secondary = dimers[-2][0], dimers[-1][0]
    deposits = dict(heads)
    if cut:
        deposits[edge(primary, secondary)] = 32
    basis = np.arange(2**n)
    bits = (basis[:, None] >> np.arange(n - 1, -1, -1)) & 1
    a, b = bits[:, :2].T
    chain = bits[:, [0, *range(2, n)]]
    diag = 32 * a + 33 * b - 32 * a * b - np.sum(chain[:, :-1] * chain[:, 1:], axis=1)
    if crossing:
        c, d = bits[:, 2:4].T
        diag += (2 * c - 1) * d
    if cut:
        diag += chain[:, -2] * chain[:, -1]
    all_terms = physical(board, dimers, deposits)
    np.testing.assert_array_equal(diag, all_terms - all_terms[0])
    h = np.diag(diag.astype(float))
    for i in range(n):
        h[basis, basis ^ (1 << i)] -= 0.5
    return h


def far_energy(e, v, z, stocks, cutoff, period, labels):
    bins = np.arange(labels // (4 * period), 2 * labels // period + 1)
    delta = e[None, :] - e[:, None]
    probe = v.T @ z @ v
    allowed = (delta > 0) & (delta < cutoff)
    false = np.zeros((len(e), len(e)), complex)
    probabilities = np.zeros((2, len(stocks)))
    # Small independent batches keep the entire 64-state instrument in memory.
    for start in range(0, len(bins), 16):
        band = bins[start : start + 16]
        for sign in (1, -1):
            maps = delay_kernel(delta, sign * band, labels=labels) * probe
            if sign == 1:
                omitted = (maps * (~allowed)).reshape(-1, len(e))
                false += omitted.conj().T @ omitted
            for i, input_amplitudes in enumerate(stocks):
                signal = maps @ input_amplitudes
                probabilities[int(sign < 0), i] += np.linalg.norm(signal) ** 2
    bound = np.linalg.eigvalsh(false)[-1]
    assert np.all(probabilities >= 0) and np.all(probabilities.sum(axis=0) < 1 + 1e-10)
    return probabilities, float(bound)


def run(lengths=(4, 3), labels=N_CLOCK, seed=26092612, crossing=False):
    board, dimers, heads, search, geometry = route(lengths, crossing)
    s, *_ = hamiltonians(4)
    se, sv = eigh(s)
    cutoff = (se[-1] - se[0]) / 10
    n = len(dimers)
    nr = 2 ** (n - 2)
    pulse = np.kron(core_pulse(sv), np.eye(nr))
    z = np.diag(1 - 2 * (np.arange(2**n) % 2))
    x_basis = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    read = np.kron(sv, np.kron(np.eye(nr // 2), x_basis))
    rng = np.random.default_rng(seed)
    period = 128
    hfull = rebuild(board, dimers, heads, False, crossing)
    hcut = rebuild(board, dimers, heads, True, crossing)
    while max(np.ptp(np.linalg.eigvalsh(h)) for h in (hfull, hcut)) >= period / 2:
        period *= 2
    panels = []
    for arm, h, excited in [
        ("joined_pulse", hfull, True),
        ("joined_idle", hfull, False),
        ("cut_pulse", hcut, True),
    ]:
        e, v = eigh(h)
        ground = v[:, 0]
        initial = pulse @ ground if excited else ground
        amp = v.T @ initial
        flight = []
        energy = []
        spectra = []
        for t in WAITS:
            evolved = np.exp(-1j * t * e) * amp
            p = abs(read.T @ (v @ evolved)).reshape(4, nr // 2, 2) ** 2
            joint = draw(rng, p.sum(axis=1).ravel())
            flight.append(
                dict(wait=t, core_by_far_x=np.asarray(joint).reshape(4, 2).tolist())
            )
            if t in (0, 16):
                spectra.append(evolved)
        probabilities, bound = far_energy(e, v, z, spectra, cutoff, period, labels)
        for t, p in zip((0, 16), probabilities.T, strict=True):
            energy.append(dict(wait=t, loss_gain_other=draw(rng, [*p, 1 - p.sum()])))
        panels.append(
            dict(
                arm=arm,
                flight=flight,
                energy=energy,
                injected_energy=round(float(abs(amp) ** 2 @ (e - e[0])), 9),
                false_loss_bound=round(bound, 11),
                joined_spectrum=np.round(e - e[0], 9).tolist(),
            )
        )
    # The receiver-alone formation spectrum uses the same four dimers. Pinning
    # a=b=0 contributes only the source scalar and leaves all of their terms.
    cold = eigh(hfull[:nr, :nr], eigvals_only=True)
    return dict(
        shots=SHOTS,
        seed=seed,
        search=search,
        dimers=dimers,
        geometry=geometry,
        native_column_residual=round(columns(), 13),
        energy_period=period,
        clock_labels=labels,
        receiver_formation=float(cold[1] - cold[0]),
        parent_formation=float(se[3] - se[0]),
        far_cut_contact=edge(dimers[-2][0], dimers[-1][0]),
        acceptance=[0.25, 2.0],
        light_cutoff=float(cutoff),
        panels=panels,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    dest = parser.parse_args().output
    rec = run()
    if dest:
        dest.write_text(json.dumps(rec, indent=1) + "\n")
    else:
        compare_records(json.loads((HERE / "records.json").read_text()), rec)
        print(
            "Verified the extended material, all near and distal controls, and signed loss records."
        )
