"""Read ordinary energy through the far dimer before any collective branch enrollment."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh, svdvals
from scipy.sparse.linalg import eigsh
from . import weak_strand as ws
from .measure import HERE
from ..joined_decay.measure import core_pulse
from ..passive_decay.measure import compare_records, hamiltonians
from ..passive_decay.readers import circuit


PERIOD = 64
LABELS = 32768
DESIGN = "STRAND_CLOCK.md"


def instrument(h, port, n, limit):
    ev, v = eigh(h)
    gaps = ev - ev[0]
    assert max(gaps) < PERIOD / 2
    bit = ws.fields(n)[0][:, port]
    phase = 1 - 2 * bit
    element = v.T @ (phase * v[:, 0])
    bins = np.fft.fftfreq(LABELS) * PERIOD
    delays = 2 * np.pi * np.arange(LABELS) / PERIOD
    kernel = np.fft.ifft(np.exp(-1j * np.outer(delays, gaps)), axis=0) * element
    low = (bins > 0) & (bins < limit)
    affordable = (gaps > 1e-10) & (gaps <= limit)
    false = svdvals(kernel[low][:, ~affordable])[0] ** 2
    masks = [bins == 0, low, bins >= limit, bins < 0]
    effects = []
    for mask in masks:
        selected = kernel[mask]
        effects.append(selected.conj().T @ selected)
    whole = np.sum(effects, axis=0)
    effects = [np.eye(len(ev)) - whole, *effects]
    return effects, v, gaps, affordable, float(false)


def observe(effects, frame, wave, core_frame):
    reduced = np.einsum("ab,acx->bcx", core_frame, wave.reshape(4, 2, -1)) @ frame
    result = np.stack(
        [np.einsum("bci,ij,bcj->b", reduced.conj(), f, reduced).real for f in effects],
        axis=1,
    )
    assert result.min() > -2e-11
    np.testing.assert_allclose(result.sum(), 1, atol=2e-11)
    return result


def inventory(b, dimers, heads):
    rings = dimers[3:]
    sites = set(sum(dimers, ()))
    report = []
    order = np.arange(4)
    for a, z in rings:
        std = next((p for p in (a, z) if b.degree[p] == 5), None)
        if std is not None:
            kind = "private"
            phase_by_load = np.array([0, 1]) ** 2
        else:
            pin = next(
                (
                    pair
                    for pair in heads
                    if bool(set(pair) & {a, z}) and len(set(pair) - sites) == 1
                ),
                None,
            )
            assert pin is not None
            std = next(iter(set(pin) & {a, z}))
            kind = f"written unit pin {pin}"
            # Normalize the written unit pin's calibrated energies by its
            # actual forbidden-word cost, then subtract from the unloaded word.
            energies = np.array([heads[pin], 0])
            phase_by_load = (energies[0] - energies) / heads[pin]
        pointers = []
        for x, y in ((0, 0), (0, 1), (1, 0), (1, 1)):
            # The target rail is delivered to the common standard by the native
            # full blank dimer swap if necessary; the returned rail is restored.
            phase = phase_by_load[x] + phase_by_load[y] - 1
            amps = np.fft.ifft(np.exp(-2j * np.pi * order * phase / 4))
            pointers.append(int(np.argmax(abs(amps))))
            np.testing.assert_allclose(amps, np.eye(4)[(x + y - 1) % 4], atol=2e-13)
        report.append(dict(dimer=[a, z], standard=kind, port=std, binary=pointers))
    return report


def acquire():
    _, core_frame = eigh(hamiltonians(4)[0])
    core_e = eigh(hamiltonians(4)[0], eigvals_only=True)
    limit = core_e[-1] - core_e[2]
    pulse = circuit(core_pulse(core_frame), 32)
    settings = {}
    for n in (6, 8):
        b, dimers, heads, _ = ws.ring_layout(n, True, False)
        native, _ = ws.constructed(dimers, b, heads, "native")
        _, v = eigsh(native, k=1, which="SA", tol=5e-14, v0=np.ones(native.shape[0]))
        idle = v[:, 0]
        active = (pulse @ idle.reshape(4, -1)).ravel()
        flights = {}
        for name in ("native", "source_cut", "no_reversal"):
            h, coupling = ws.constructed(dimers, b, heads, name)
            after = ws.propagate(h, np.column_stack((active, idle)), 32)
            flights[name] = [
                after,
                [-np.vdot(s, coupling * s).real for s in after.T],
                float(np.vdot(active, h @ active).real - np.vdot(idle, h @ idle).real),
            ]
        observations = {}
        for port in (0, n - 1):
            effects, frame, gaps, bands, false = instrument(
                ws.strand(n).toarray(), port, n, limit
            )
            far = effects[2]
            near_phase = 1 - 2 * ws.fields(n)[0][:, 0]
            known = frame.T @ (near_phase * frame[:, 0])
            marks = known * bands
            rows = {}
            for name, (after, _, _) in flights.items():
                banks = {}
                for label, vec in zip(("pulse", "idle"), after.T, strict=True):
                    result = observe(effects, frame, vec, core_frame)
                    coeff = (
                        np.einsum("ab,acx->bcx", core_frame, vec.reshape(4, 2, -1))
                        @ frame
                    )
                    outside = coeff * ~bands
                    # The measured open-system intervention captures all affordable
                    # components in orthogonal parked records; the replaced state
                    # is the strand vacuum. All other amplitudes retain coherence.
                    no_packet = np.einsum(
                        "bci,ij,bcj->b", outside.conj(), far, outside
                    ).real
                    parked = (abs(coeff[:, :, bands]) ** 2).sum(axis=(1, 2))
                    with_reset = no_packet + far[0, 0].real * parked
                    banks[label] = dict(
                        core_by_bin=result.tolist(), reset_core_soft=with_reset.tolist()
                    )
                rows[name] = banks
            observations[str(port)] = dict(
                false_upper=false,
                known_soft=float(np.vdot(known, far @ known).real),
                known_in_band_soft=float(np.vdot(marks, far @ marks).real),
                ground_soft=float(far[0, 0].real),
                rows=rows,
            )
        settings[str(n)] = dict(
            inventory=inventory(b, dimers, heads),
            results=observations,
            pulse_and_idle_cut_work={
                k: [float(x) for x in w] for k, (_, w, _) in flights.items()
            },
            pulse_work={k: w for k, (_, _, w) in flights.items()},
        )
    return dict(
        protocol_sha256=sha256((HERE / DESIGN).read_bytes()).hexdigest(),
        clock=dict(
            labels=LABELS,
            period=PERIOD,
            duration=(LABELS - 1) * 2 * np.pi / PERIOD,
            affordable_upper=limit,
            columns=[
                "vacuum_refused",
                "zero",
                "affordable",
                "higher_positive",
                "negative",
            ],
        ),
        settings=settings,
    )


def decode(data):
    for n, row in data["settings"].items():
        end = row["results"][str(int(n) - 1)]
        print(
            n,
            "false bound",
            end["false_upper"],
            "known",
            end["known_soft"],
            "ground",
            end["ground_soft"],
        )
        for key, banks in end["rows"].items():
            for arm, rec in banks.items():
                x = np.array(rec["core_by_bin"])
                reset = np.array(rec["reset_core_soft"])
                print(
                    key,
                    arm,
                    "soft, joint lower, with reset, conditioned false limit",
                    np.round(
                        [
                            x[:, 2].sum(),
                            x[1:3, 2].sum(),
                            reset[1:3].sum(),
                            x[1:3].sum() * end["false_upper"],
                        ],
                        9,
                    ),
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    result = acquire()
    decode(result)
    if args.acquire:
        with args.acquire.open("x") as f:
            json.dump(result, f, allow_nan=False, separators=(",", ":"))
            f.write("\n")
    else:
        compare_records(json.loads((HERE / "strand_clock.json").read_text()), result)
