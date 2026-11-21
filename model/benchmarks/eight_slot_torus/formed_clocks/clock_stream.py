"""Unresolved ordinary births measured by one installed age-controlled clock factory."""
from __future__ import annotations

import argparse
import io
import hashlib
import json
from pathlib import Path

import numpy as np

from . import clock_factory as factory
from . import clock_loader as loader
from . import clock_records
from .. import formation, law
from .loaded_clock import PORTS


HERE = Path(__file__).parent
PROTOCOL = HERE / "clock_factory_protocol.json"
ARCHIVE = HERE / "clock_stream.npz"
MENUS = ("installed", "passive", "capture-off")
OUTCOMES = ("+x", "+y", "-x", "-y", "local below loading age", "marked elsewhere", "unmarked")
SHOTS = 1 << 20
SEED = 391173
PUBLISHED = "10fd939b3b7662506ff698cc535d7c66930fa5f8e78345709f714e1ac0040d98"


def ordinary_history(shape, root, focus, ticks, angle):
    """Retain every native first write; later histories stay orthogonal on the wires."""
    a, b = loader.ordinary_feeds(shape, root, *focus)
    free = np.multiply.outer(a, b)
    writes = np.empty((ticks, *shape))
    unmarked = np.empty(ticks+1)
    rail_error = 0.
    for t in range(ticks+1):
        unmarked[t] = np.vdot(free, free).real
        if t == ticks:
            break
        free, born = formation.first_writes(free, angle)
        writes[t] = abs(born)**2
        step0, step1 = (law.step_contact(free, phi) for phi in loader.SOURCE.rail_contact_phases)
        rail_error = max(rail_error, float(np.max(abs(step0-step1))))
        free = step0
    return writes, unmarked, rail_error


def _translated_overlaps(waves):
    """Overlap at every origin with each of the two normalized signed loops."""
    even = (waves[..., 0]-waves[..., 2] + np.roll(waves[..., 1]-waves[..., 3], (-1,-1), (1,2))) / 2
    odd = (np.roll(waves[..., 1]-waves[..., 2], -1, axis=1)
           + np.roll(waves[..., 0]-waves[..., 3], -1, axis=2)) / 2
    return ((even+odd)/np.sqrt(2), (even-odd)/np.sqrt(2))


def impulse(shape, ticks, installed):
    """One actual birth with a balanced rail preparation, including all premature reads."""
    states = np.zeros((2, *shape, law.SLOTS), complex)
    states[(slice(None), 0, 0, 0, 0)] = law.EQUAL/np.sqrt(2)
    ports = np.empty((ticks+1, *shape, 4))
    gates = factory.uniform_schedule() if installed else {}
    dark_overlap = 0.
    for age in range(ticks+1):
        if not installed:
            for overlaps in _translated_overlaps(states):
                dark_overlap = max(dark_overlap, float(abs(overlaps).max()))
        for i, phase in enumerate(PORTS):
            ports[age, ..., i] = np.sum(abs(states[0]+phase*states[1])**2, axis=-1)/4
        if age < ticks:
            states = factory.advance(states, age, gates)
    return ports, dark_overlap


def _convolve(weight, kernel):
    """Convolution only in hardware position. First-wire ages add in probability."""
    axes = (0, 1, 2, 3)
    product = np.fft.fftn(weight, axes=axes)[..., None] * np.fft.fftn(kernel, axes=axes)
    return np.fft.ifftn(product, axes=axes).real


def record(writes, unmarked, impulse_ports, tick):
    full = np.zeros((*writes.shape[1:], len(OUTCOMES)))
    for born_tick, weights in enumerate(writes[:tick]):
        age = tick-born_tick
        local = _convolve(weights, impulse_ports[age])
        if age >= factory.END:
            full[..., :4] += local
        else:
            full[..., 4] += local.sum(axis=-1)
    full[..., 5] = writes[:tick].sum() - full[..., :5].sum(axis=-1)
    full[..., 6] = unmarked[tick]
    return full


def parity_forecast(writes, tick):
    """Observer-port forecast from first-wire powers and two compact routing parities."""
    quadratures = np.zeros((*writes.shape[1:], 4))
    mu = np.array(loader.SOURCE.rail_contact_phases)/2
    delta = mu[1]-mu[0]
    phase_sum = np.r_[0., np.cumsum([law.hop_angle(a, loader.SOURCE.maturation)
                                     for a in range(factory.END, tick)])]
    for born_tick, w in enumerate(writes[:max(0, tick-factory.END+1)]):
        elapsed = tick-born_tick-factory.END
        phase = phase_sum[elapsed]
        even = w + np.roll(w, (1, 1), (0, 1))
        odd = np.roll(w, 1, axis=0) + np.roll(w, 1, axis=1)
        incidence = (np.cos(phase)**2*even + np.sin(phase)**2*odd)/2
        hand = np.exp(1j*delta*(1-elapsed))
        quadratures += incidence[..., None] * (1 + (PORTS.conj()*hand).real)/4
    return quadratures


def produce(protocol):
    shape, root = loader.chart(protocol["stream_carrier"], protocol["stream_root"])
    focus, ticks = protocol["stream_focus"], protocol["stream_ticks"]
    if max(ticks) >= min(protocol["wire_period"], 8192):
        raise ValueError("the first-write reduction ends at the wire or age return")
    mass, free, rail_error = ordinary_history(shape, root, focus, max(ticks), formation.ANGLE)
    null, bypass, null_error = ordinary_history(shape, root, focus, max(ticks), 0.)
    active, _ = impulse(shape, max(ticks), True)
    passive, dark = impulse(shape, max(ticks), False)
    rows = np.stack([[record(m, u, port, t) for t in ticks] for m, u, port
                     in ((mass, free, active), (mass, free, passive), (null, bypass, active))])
    forecast = np.stack([parity_forecast(mass, t) for t in ticks])
    source_norm = np.r_[0., np.cumsum(mass.reshape(len(mass),-1).sum(axis=1))]+free
    failures = []
    def gate(ok, label):
        if not ok:
            failures.append(label)
    gate(max(rail_error, null_error) <= 2e-11, "unmarked rail agreement")
    gate(np.max(abs(source_norm-1)) <= 2e-11 and np.max(abs(bypass-1)) <= 2e-11, "source census")
    gate(rows.min() >= -2e-12 and np.max(abs(rows.sum(axis=-1)-1)) <= 2e-11, "receiver census")
    gate(np.max(abs(rows[0,...,:4]-forecast)) <= 2e-11, "native ports versus parity transfer")
    gate(np.count_nonzero(null) == 0 and np.max(abs(rows[2,...,:6])) <= 2e-11, "disabled capture")
    gate(dark <= 2e-11, "passive projection onto translated loops")
    result = {"writes": mass, "unmarked": free, "words": rows}
    epoch_power = mass.reshape(len(mass),-1).sum(axis=1)
    stats = {"wire_sites_above_1e_12": int(np.count_nonzero(mass.sum(axis=0)>1e-12)),
             "write_epochs_above_1e_12": int(np.count_nonzero(epoch_power>1e-12)),
             "timing_null_tick_2": float(epoch_power[2]),
             "max_parity_forecast_error": float(np.max(abs(rows[0,...,:4]-forecast))),
             "max_source_norm_error": float(abs(source_norm-1).max()),
             "max_unmarked_rail_error": max(rail_error, null_error),
             "max_passive_loop_amplitude": dark, "failures": failures}
    return result, stats


def sample(words):
    counts = np.empty_like(words, dtype=np.int64)
    for key in np.ndindex(words.shape[:-1]):
        probabilities = np.maximum(words[key], 0.)
        probabilities /= probabilities.sum()
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED, *key])))
        counts[key] = rng.multinomial(SHOTS, probabilities)
    return counts


def experiment(path=ARCHIVE, *, acquire=False, seed_replay=False):
    path = Path(path)
    blob = None if acquire else clock_records.snapshot(path, PUBLISHED)
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    result, stats = produce(protocol)
    if acquire:
        blob = clock_records.snapshot(path, PUBLISHED, arrays={**result, "counts": sample(result["words"])})
    with np.load(io.BytesIO(blob), allow_pickle=False) as saved:
        if sorted(saved.files) != sorted([*result, "counts"]):
            raise ValueError("unexpected stream archive layout")
        for key, expected in result.items():
            actual = saved[key]
            if not clock_records.finite_float_record(actual, expected, 2e-11):
                raise ValueError(f"stream replay mismatch in {key}")
        counts = saved["counts"]
        if (counts.dtype != np.int64 or counts.shape != saved["words"].shape or (counts<0).any()
                or (counts>SHOTS).any() or (counts.sum(axis=-1)!=SHOTS).any()
                or (counts[saved["words"] <= 0]!=0).any()):
            raise ValueError("inconsistent stream counts or a forbidden outcome")
        if seed_replay and not np.array_equal(counts, sample(saved["words"])):
            raise ValueError("stream counts disagree with the local seeded acquisition")
    words = result["words"]
    max_cell = float(abs(counts/SHOTS-words).max())
    radius = np.sqrt(np.log(2*words.size/0.001)/(2*SHOTS)) + 2e-11
    if max_cell > radius:
        stats["failures"].append("simultaneous seven-port count envelope")
    summary = []
    for i, tick in enumerate(protocol["stream_ticks"]):
        local = words[0,i,...,:4]
        power = local.sum(axis=-1)
        root = tuple(protocol["stream_root"])
        z = complex(local[root][0]-local[root][2], local[root][1]-local[root][3])
        observed = counts[(0,i,*root)][:4]
        obs_hand = complex(int(observed[0])-int(observed[2]), int(observed[1])-int(observed[3]))
        summary.append({"tick": tick, "all_marked": float(1-result["unmarked"][tick]),
                        "all_loaded": float(power.sum()), "root_loaded": float(power[root]),
                        "root_source_hand": float(np.angle(z)) if power[root] > 0 else None,
                        "root_count_hand": float(np.angle(obs_hand)) if abs(obs_hand)>0 else None})
    stats.update({"appointments": summary, "trials": int(counts.sum()), "max_count_cell_error": max_cell,
                  "family_cell_radius": float(radius), "archive_sha256": hashlib.sha256(blob).hexdigest()})
    return stats


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=ARCHIVE)
    parser.add_argument("--acquire", action="store_true")
    parser.add_argument("--seed-replay", action="store_true")
    args = parser.parse_args()
    output = experiment(args.archive, acquire=args.acquire, seed_replay=args.seed_replay)
    print(json.dumps(output, indent=2))
    raise SystemExit(int(bool(output["failures"])))
