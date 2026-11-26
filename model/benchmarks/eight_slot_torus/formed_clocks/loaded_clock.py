"""All first-wire and terminal records of the locally loaded pair clock.

Each setting is a fresh source. Counts are simulated independent launches.
The model supplies the probabilities; the scorer sees aggregate pointer
words and the independently specified clock forecast.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np

from . import clock_loader as loader
from .. import law
from .clock_modes import contact_fibre


HERE = Path(__file__).parent
BOOK = HERE / "clock_loader_protocol.json"
ARCHIVE = HERE / "loaded_clock_counts.npz"
PUBLISHED_SHA256 = "a96ff88a082ad09a1285c01976677099ecd2142547a1118509174e299029a761"
GRID_SIZE = 1 << 38
SOURCE_TOLERANCE = 2e-11
PORTS = np.array([1, 1j, -1, -1j])
LABELS = ("+x", "+y", "-x", "-y", "born elsewhere in detector chart", "other first wire", "unmarked")
READS = (*loader.MODES, "reference", "dephased")
DELTA_MU = np.diff(loader.SOURCE.rail_contact_phases)[0] / 2


def channel(images, target):
    """Worst entangled-input and Bell fidelities with the ideal encoded rail."""
    overlap = np.einsum("i,ri->r", target.conj().ravel(), images.reshape(2, -1))
    corrected = overlap * np.exp(0.5j * np.array(loader.SOURCE.rail_contact_phases))
    change = corrected[1] - corrected[0]
    p = np.clip(-(corrected[0].conjugate() * change).real / abs(change)**2, 0, 1) if change else 0.
    return {"target_power": (abs(overlap)**2).tolist(),
            "worst_fidelity": float(abs(corrected[0] + p * change)**2),
            "bell_fidelity": float(abs(corrected.mean())**2)}


def mature_powers(carrier, max_wait):
    """Binary powers in every Fourier cell, with all eight slot directions."""
    frequencies = [2*np.pi*np.fft.fftfreq(n) for n in carrier]
    routing = np.stack([contact_fibre([f[i] for i, f in zip(index, frequencies)])[0]
                        for index in np.ndindex(carrier)])
    coins = np.stack([law.contact_coin(loader.IDENTITY, phi).T for phi in loader.SOURCE.rail_contact_phases])
    operator = routing[None] @ coins[:, None]
    powers = []
    for _ in range(max_wait.bit_length()):
        powers.append(operator)
        operator = operator @ operator
    return powers


def at_appointments(states, ages, powers, *, start_age=loader.LOAD_AGE):
    """Literal open hops through maturity, then exact finite-torus unitary powers."""
    current = start_age
    if list(ages) != sorted(set(ages)) or any(not current <= a < 8192 for a in ages):
        raise ValueError("increasing distinct post-load, pre-wrap ages required")
    axes = (1, 2, 3, 4)
    wave = None
    for requested in ages:
        while current < min(requested, loader.SOURCE.maturation[2]):
            states = np.stack([law.step_maturing(s, loader.SOURCE, current, r) for r, s in enumerate(states)])
            current += 1
        gap = requested - current
        if gap:
            if wave is None:
                wave = np.fft.fftn(states, axes=axes).reshape(2, -1, law.SLOTS)
            for bit in range(gap.bit_length()):
                if gap & (1 << bit):
                    wave = (powers[bit] @ wave[..., None])[..., 0]
            states = np.fft.ifftn(wave.reshape(states.shape), axes=axes)
            current = requested
        yield states


def local_ports(components, dephased=False):
    """Four exchange or rail outcomes on co-located amplitudes, including all slots."""
    if dephased:
        return np.full(4, np.vdot(components, components).real / 4)
    probes = components[0] + PORTS[:, None] * components[1]
    return np.sum(abs(probes)**2, axis=1) / 4


def original_words(birth, ports):
    """The local aperture, its complement and both first-wire refusals."""
    word = np.r_[birth["natal"] * ports, birth["natal"] * (1 - ports.sum()), birth["other"], birth["unmarked"]]
    if not np.isfinite(word).all() or min(word) < -2e-11 or abs(word.sum() - 1) > 2e-11:
        raise ArithmeticError("invalid original-launch distribution")
    return word


def _grid(word):
    # A finite categorical source: stable largest-remainder apportionment,
    # with a per-cell error below one grid unit for normalized input.
    scale = GRID_SIZE
    raw = np.maximum(word, 0) * scale
    weights = np.floor(raw).astype(np.int64)
    left = scale - int(weights.sum())
    if not 0 <= left <= len(word):
        raise ArithmeticError("categorical rounding changed the normalization")
    weights[np.argsort(-(raw - weights), kind="stable")[:left]] += 1
    return weights


def design(protocol):
    """Emit fully enumerated original-start probability words in public label order."""
    ages = protocol["terminal_ages"]
    focus = protocol["focus_fractions"]
    rows, keys, channels = [], [], {}
    for k, fixture in enumerate(protocol["fixtures"]):
        carrier, root = loader.chart(fixture["carrier"], fixture["root"])
        sites = loader.corners(carrier, root)
        target = loader.circulation(carrier, root)
        powers = mature_powers(carrier, max(ages) - 128)
        births = [loader.ordinary_launch(carrier, root, *fractions) for fractions in focus]
        for menu in (*loader.MODES, "reference"):
            reference = menu == "reference"
            mode = "exact" if reference else menu
            states, _ = loader.prepare_clock(carrier, root, mode, local_reference=reference)
            if not reference:
                channels[f"{k}/{menu}"] = channel(states, target)
                states /= np.sqrt(2)
            for age, waves in zip(ages, at_appointments(states, ages, powers)):
                for corner, site in enumerate(sites[:1] if reference else sites):
                    components = waves[(slice(None), *site)]
                    for read in (("reference", "dephased") if reference else (menu,)):
                        ports = local_ports(components, dephased=(read == "dephased"))
                        for width, birth in enumerate(births):
                            keys.append((k, READS.index(read), width, corner, age))
                            rows.append(original_words(birth, ports))
    return np.asarray(keys, dtype=np.int32), np.asarray(rows), channels


def counts_for(keys, weights, protocol):
    counts = np.empty_like(weights)
    for i, key in enumerate(keys):
        sequence = np.random.SeedSequence([protocol["count_seed"], *map(int, key)])
        generator = np.random.Generator(np.random.PCG64(sequence))
        counts[i] = generator.multinomial(protocol["trials_per_setting"], weights[i] / GRID_SIZE)
    return counts


def _hand_disk(port_counts, log_factor, trials=None, cell_radius=0.):
    """Disk for the counting law's quadratures, optionally covering the native direction."""
    n = sum(map(int, port_counts))
    if n == 0:
        return 0j, np.inf
    x, y, mx, my = map(int, port_counts)
    z = 2 * complex(x - mx, y - my) / n
    # Each empirical component is a mean of numbers in [-2,2]. The two
    # coordinate bounds imply a disk around z. A disk enclosing zero has
    # all directions. Conditional coverage holds separately for each n.
    radius = np.sqrt(16 * log_factor / n)
    if trials is not None:
        # On the simultaneous count-cell event the grid's local probability
        # g is >= n/trials - 4*cell_radius. The native signed numerator,
        # scaled by 2/g, has the native direction. Its displacement from
        # the grid center is <= 4 sqrt(2) times the per-cell tolerance / g.
        lower = n / trials - 4 * cell_radius
        radius += 4 * np.sqrt(2) * SOURCE_TOLERANCE / lower if lower > 0 else np.inf
    return z, float(radius)


def hand_interval(port_counts, log_factor, *, trials=None, cell_radius=0.):
    """Return the signed count angle and simultaneous angular half-width."""
    z, radius = _hand_disk(port_counts, log_factor, trials, cell_radius)
    return float(np.angle(z)), float(np.arcsin(radius / abs(z)) if radius < abs(z) else np.pi)


def score(keys, words, counts, channels, protocol):
    trials = protocol["trials_per_setting"]
    if (counts.shape != words.shape or counts.dtype.kind not in "iu" or (counts < 0).any()
            or (counts > trials).any()
            or (counts.sum(axis=1) != trials).any() or (counts[words == 0] != 0).any()):
        raise ValueError("inconsistent seven-outcome counts")
    alpha = protocol["family_error"]
    # Split the family error equally between original-start cells and hands.
    grid_cell_radius = np.sqrt(np.log(4 * words.size / alpha) / (2 * trials))
    cell_radius = grid_cell_radius + SOURCE_TOLERANCE
    log_factor = np.log(8 * len(keys) / alpha)
    delta = 27 * 2 * np.sin(law.DELTA / 8)
    mismatches = []
    def gate(ok, name):
        if not ok:
            mismatches.append(name)
    for key, item in channels.items():
        mode = key.split("/")[1]
        if mode == "exact":
            gate(item["worst_fidelity"] >= 1-2e-11, key)
        elif mode == "common":
            gate(item["worst_fidelity"] >= (1-delta**2/2)**2-2e-11, key)
        elif mode == "phase-cut":
            gate(max(item["target_power"]) <= 0.999, key)
        elif mode == "passive":
            gate(max(item["target_power"]) <= 2e-22, key)
    gate(np.all(np.abs(counts / trials - words) <= cell_radius), "categorical coverage")
    max_error = max_radius = aperture_error = phase_error = 0.
    local_intervals = informative = dephased_disks = dephased_include_zero = 0
    max_dephased_visibility = 0.
    for key, word, count in zip(keys, words, counts):
        k, m, w, c, age = map(int, key)
        menu = READS[m]
        if menu == "dephased":
            gate(np.ptp(word[:4]) < 2e-11, f"which-rail {key}")
            center, disk = _hand_disk(count[:4], log_factor, trials, grid_cell_radius)
            dephased_disks += 1
            dephased_include_zero += int(abs(center) <= disk)
            max_dephased_visibility = max(max_dephased_visibility, abs(center))
        if menu not in ("exact", "reference"):
            continue
        a, b = protocol["focus_fractions"][w]
        aperture_error = max(aperture_error, abs(sum(word[:4])/(a*b*0.5)-0.25))
        phase = DELTA_MU * (loader.LOAD_AGE - age + int(menu == "exact"))
        hand = complex(word[0] - word[2], word[1] - word[3])
        phase_error = max(phase_error, abs(np.angle(hand * np.exp(-1j*phase))))
        estimated, radius = hand_interval(count[:4], log_factor, trials=trials, cell_radius=grid_cell_radius)
        error = abs(np.angle(np.exp(1j*(estimated-phase))))
        gate(error <= radius, f"angle {key}")
        max_radius, max_error = max(max_radius, radius), max(max_error, error)
        local_intervals += 1
        informative += int(radius < np.pi)
    gate(aperture_error < 2e-11 and phase_error < 2e-10, "analytic local hand and aperture")
    for fixture in protocol["fixtures"]:
        for fractions in protocol["focus_fractions"]:
            null = loader.ordinary_launch(fixture["carrier"], fixture["root"], *fractions, angle=0.)
            gate(null == {"natal": 0., "other": 0., "unmarked": 1.}, "disabled capture")
    return {"settings": len(keys), "launches": int(counts.sum()), "channels": channels,
            "count_hand_intervals": local_intervals, "informative_count_hands": informative,
            "full_circle_count_hands": local_intervals - informative,
            "dephased_disks": dephased_disks, "dephased_disks_including_zero": dephased_include_zero,
            "max_dephased_count_visibility": max_dephased_visibility, "max_count_phase_error": max_error,
            "max_count_radius": max_radius, "max_aperture_error": aperture_error,
            "max_born_phase_error": phase_error, "max_count_cell_error": float(abs(counts/trials-words).max()),
            "simultaneous_cell_radius": float(cell_radius), "failures": mismatches}


def _records(payload, keys, words, trials):
    with np.load(io.BytesIO(payload), allow_pickle=False) as saved:
        if (sorted(saved.files) != ["counts", "keys", "weights"] or saved["keys"].dtype != np.dtype("<i4")
                or not np.array_equal(keys, saved["keys"])):
            raise ValueError("incorrect setting IDs or missing record columns")
        weights, counts = saved["weights"], saved["counts"]
    if (weights.dtype != np.dtype("<i8") or weights.shape != words.shape or (weights < 0).any()
            or (weights > GRID_SIZE).any() or (weights.sum(axis=1) != GRID_SIZE).any()
            or abs(words - weights/GRID_SIZE).max() > SOURCE_TOLERANCE):
        raise ValueError("archived integer source disagrees with the native replay")
    if (counts.dtype != np.dtype("<i8") or counts.shape != words.shape or (counts < 0).any()
            or (counts > trials).any() or (counts.sum(axis=1) != trials).any()
            or (counts[weights == 0] != 0).any()):
        raise ValueError("inconsistent counts for the archived integer source")
    return weights, counts


def experiment(archive=ARCHIVE, acquire=False, seed_replay=False, *, compare_seed=False):
    protocol = json.loads(BOOK.read_text(encoding="utf-8"))
    archive = Path(archive)
    # Authenticate the very bytes used to decode the records. Replaying a
    # different acquisition requires an explicit seed check and a local regrid.
    if not acquire:
        payload = archive.read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        if digest != PUBLISHED_SHA256 and not seed_replay:
            raise ValueError("archive SHA-256 differs from the published simulated counts")
    keys, words, channels = design(protocol)
    if acquire:
        weights = np.stack([_grid(w) for w in words])
        counts = counts_for(keys, weights, protocol)
        buffer = io.BytesIO()
        np.savez_compressed(buffer, keys=keys, weights=weights, counts=counts)
        payload = buffer.getvalue()
        digest = hashlib.sha256(payload).hexdigest()
        with archive.open("xb") as output:
            output.write(payload)
    weights, counts = _records(payload, keys, words, protocol["trials_per_setting"])
    if not acquire and digest != PUBLISHED_SHA256 and not np.array_equal(weights, [_grid(w) for w in words]):
        raise ValueError("unpublished integer source disagrees with the local acquisition grid")
    result = score(keys, words, counts, channels, protocol)
    if seed_replay or compare_seed:
        difference = counts_for(keys, weights, protocol) - counts
        if seed_replay and np.any(difference):
            raise ValueError("the seeded acquisition disagrees with the archive")
        result["same_seed_replay"] = {
            "differing_rows": int(np.any(difference != 0, axis=1).sum()),
            "max_bin_difference": int(np.max(np.abs(difference))),
        }
    result["archive_sha256"] = digest
    result["count_source"] = ("published simulated acquisition" if digest == PUBLISHED_SHA256
                              else "local simulated acquisition")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=ARCHIVE)
    parser.add_argument("--acquire", action="store_true", help="create counts at an unused path")
    parser.add_argument("--seed-replay", action="store_true", help="redraw the stored counts; also verifies a local acquisition")
    parser.add_argument("--compare-seed", action="store_true",
                        help="report backend-specific multinomial differences on the authenticated archive")
    args = parser.parse_args()
    result = experiment(args.archive, args.acquire, args.seed_replay, compare_seed=args.compare_seed)
    print(json.dumps(result, indent=2))
    return int(bool(result["failures"]))


if __name__ == "__main__":
    raise SystemExit(main())
