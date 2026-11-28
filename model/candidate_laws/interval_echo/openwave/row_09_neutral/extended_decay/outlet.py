"""Release the actual downstream guide loads into a different owner's bare cell."""

import argparse
from collections import defaultdict, deque
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm

from model.candidate_laws.interval_echo import candidate
from ..joined_decay.measure import core_pulse
from ..passive_decay.measure import compare_records, hamiltonians
from ..passive_decay.readers import circuit
from . import wing
from .material import edge, route
from .measure import HERE

M = 4096
PERIOD = 256
SOURCE_WAIT = 32
WAITS = (0, 4, 8, 16, 32, 64)
OUTCOMES = (
    "charge_minus_one",
    "neutral_other",
    "neutral_light",
    "neutral_heavy",
    "charge_plus_one",
)


def outlet_geometry():
    board, dimers, heads, _, _ = route((3,), crossing=True)
    owners = {p // 12 for pair in dimers for p in pair}
    forbidden = {p for pair in heads for p in pair} | {
        p for pair in dimers for p in pair
    }
    queue = deque((p,) for p in sorted(dimers[-1]))
    while queue:
        path = queue.popleft()
        p = path[-1]
        if p // 12 not in owners:
            break
        queue.extend(
            path + (q,)
            for q in sorted(board.neighbours[p])
            if q not in forbidden and q not in path
        )
    else:
        raise ValueError("No clean outlet")
    cell = tuple(range(12 * (p // 12), 12 * (p // 12 + 1)))
    internal = [
        edge(a, b) for a in cell for b in board.neighbours[a] if b in cell and a < b
    ]
    assert not set(cell) & forbidden and len(internal) == 30
    return board, dimers, path, cell, internal


def binary_generator(board, sites, contacts, number):
    """Sum the blank 99-state law in its invariant binary, fixed-inventory sector."""
    words = list(combinations(sites, number))
    lookup = {frozenset(w): j for j, w in enumerate(words)}
    exchange, contrast = candidate.native_law.local_term_parts()
    bond = exchange + candidate.KAPPA * contrast
    entries = candidate.native_law.SRC
    for a, b in ((0, 0), (0, 1), (1, 0), (1, 1)):
        col = entries.index((a, b))
        prediction = np.zeros(len(bond))
        prediction[col] += candidate.KAPPA * (a - b) ** 2 + 0.5
        prediction[entries.index((b, a))] -= 0.5
        np.testing.assert_allclose(bond[:, col], prediction, atol=2e-13)
    h = np.zeros((len(words), len(words)))
    for i, w in enumerate(words):
        occupied = set(w)
        h[i, i] = sum(candidate.KAPPA * (6 - board.degree[p]) for p in w)
        for a, b in contacts:
            bits = (int(a in occupied), int(b in occupied))
            col = entries.index(bits)
            h[i, i] += bond[col, col]
            if bits[0] != bits[1]:
                j = lookup[frozenset(occupied ^ {a, b})]
                h[j, i] += bond[entries.index(bits[::-1]), col]
    return words, h


def material(board, dimers, path, cell, internal, closed):
    e0, e1 = dimers[-1]
    sites = tuple(sorted({e0, e1, *path, *cell}))
    opened = list(zip(path, path[1:]))
    contacts = [edge(e0, e1), *internal, *(opened[1:] if closed else opened)]
    return connected_wing(board, dimers, sites, contacts)


def connected_wing(board, dimers, sites, contacts, number=2):
    e0, e1 = dimers[-1]
    words, mobile = binary_generator(board, sites, contacts, number)
    n = len(words)
    # D retains one load. Its blank internal exchange has off-diagonal -1/2.
    # Its D,E triangle and complementary D-to-unit pin remain open for any
    # occupancy of E's rails, including empty and doubly occupied across E.
    h = np.kron(np.array([[0.0, -0.5], [-0.5, 0.0]]), np.eye(n)) + np.kron(
        np.eye(2), mobile
    )
    for d in (0, 1):
        for j, word in enumerate(words):
            h[d * n + j, d * n + j] += d * (
                32 - 33 * int(e0 in word) - 32 * int(e1 in word)
            )
    return words, h


def receiving_cell(board, cell, internal):
    sectors = {}
    for number in (0, 1, 2):
        words, h = binary_generator(board, cell, internal, number)
        values, vectors = eigh(h)
        sectors[number] = (words, values, vectors)
    inventory_reader(board, cell, sectors)
    return sectors


def relative_inventory_reader(board, cell, sectors, labels, reference):
    """Sum private-phase comparisons of binary loads against a reference inventory."""
    ex, contrast = candidate.native_law.local_term_parts()
    local = ex + candidate.KAPPA * contrast
    order = candidate.native_law.SRC
    swap = np.eye(9)[[order.index((b, a)) for a, b in order]]
    np.testing.assert_allclose(
        expm(-1j * np.pi * local)[:, :9],
        np.vstack((swap, np.zeros((len(local) - 9, 9)))),
        atol=3e-13,
    )
    routes = {}
    for port in cell:
        if board.degree[port] == 5:
            routes[port] = [port]
        else:
            meter = min(
                q for q in board.neighbours[port] if q in cell and board.degree[q] == 5
            )
            routes[port] = [port, meter]
    peaks = []
    for count, (words, _, _) in sectors.items():
        for word in words:
            angles = -2j * np.pi * np.arange(labels) / labels
            phase = np.exp(-angles * reference)
            for port, path in routes.items():
                moved = set(word)
                if len(path) == 2 and ((port in moved) != (path[-1] in moved)):
                    moved.symmetric_difference_update(path)
                phase *= np.exp(angles * int(path[-1] in moved))
            amplitudes = np.fft.ifft(phase)
            target = np.eye(labels)[(count - reference) % labels]
            np.testing.assert_allclose(amplitudes, target, atol=2e-13)
        peaks.append((count - reference) % labels)
    return {str(p): path for p, path in routes.items()}, peaks


def inventory_reader(board, cell, sectors):
    routes, peaks = relative_inventory_reader(board, cell, sectors, 4, 0)
    return dict(port_to_private=routes, inventory_zero_one_two_to_pointer=peaks)


def clock(values, reference, threshold, labels=M, low_edge=0.25):
    energies = values - reference
    assert np.max(abs(energies)) < PERIOD / 2
    bins = np.fft.fftfreq(labels) * PERIOD
    delays = 2 * np.pi * np.arange(labels) / PERIOD
    waves = np.fft.ifft(np.exp(-1j * energies[:, None] * delays), axis=1)
    probabilities = abs(waves) ** 2
    light = (bins >= low_edge) & (bins < threshold)
    heavy = bins >= threshold
    categories = np.column_stack(
        [
            probabilities[:, ~(light | heavy)].sum(1),
            probabilities[:, light].sum(1),
            probabilities[:, heavy].sum(1),
        ]
    )
    np.testing.assert_allclose(categories.sum(1), 1, atol=1e-12)
    return categories


def read_groups(words, sectors, cell):
    cell = set(cell)
    groups = defaultdict(list)
    for d in (0, 1):
        for i, word in enumerate(words):
            inside = tuple(p for p in word if p in cell)
            outside = tuple(p for p in word if p not in cell)
            groups[(len(inside), d, outside)].append((inside, d * len(words) + i))
    result = []
    for (number, _, _), items in groups.items():
        labels, _, vectors = sectors[number]
        rows = dict(items)
        result.append((number, [rows[word] for word in labels], vectors))
    return result


def flight_read(waves, groups, sectors, response):
    # Columns retain the physically read core level and the unobserved C rail.
    weights = {n: np.zeros((4, len(values))) for n, (_, values, _) in sectors.items()}
    for count, rows, vectors in groups:
        projected = vectors.T @ waves[rows]
        weights[count] += (abs(projected.reshape(len(rows), 4, -1)) ** 2).sum(2).T
    bands = np.column_stack(
        [weights[0].sum(1), weights[1] @ response, weights[2].sum(1)]
    )
    np.testing.assert_allclose(bands.sum(), 1.0, atol=2e-12)
    # Only the total weight in a degenerate line is independent of the solver's frame.
    spectra = {}
    for n, (_, values, _) in sectors.items():
        starts = np.r_[0, np.flatnonzero(np.diff(values) > 1e-8) + 1]
        spectra[str(n)] = np.add.reduceat(weights[n], starts, axis=1).round(13).tolist()
    return dict(core_by_outlet=bands.round(13).tolist(), core_by_sector_energy=spectra)


def source_stocks():
    _, _, settings = wing.installation()
    core = hamiltonians(4)[0]
    core_values, cv = eigh(core)
    _, basis = eigh(settings["native"])
    cold = basis[:, 0]
    pulse = np.kron(circuit(core_pulse(cv), 32), np.eye(8)) @ cold
    reader = circuit(cv, 32)
    stocks = {}
    for name, h in settings.items():
        values, v = eigh(h)
        for arm, initial in (("pulse", pulse), ("idle", cold)):
            psi = v @ (np.exp(-1j * SOURCE_WAIT * values) * (v.T @ initial))
            stocks[name, arm] = (np.kron(reader.T, np.eye(8)) @ psi).reshape(4, 2, 2, 2)
    return stocks, float(np.ptp(core_values))


def inject(guide, words, dimers, ports, ready):
    """Embed the unchanged wing and an independent ready ensemble into physical port words."""
    lookup = {frozenset(w): i for i, w in enumerate(words)}
    n = len(words)
    output = np.zeros((2 * n, 4, 2, ready.shape[1]), complex)
    for d in (0, 1):
        for e in (0, 1):
            for q, amp in zip(ports, ready, strict=True):
                i = d * n + lookup[frozenset((q, dimers[-1][1 - e]))]
                output[i] = guide[:, :, d, e, None] * amp
    return output.reshape(2 * n, -1)


def acquire():
    board, dimers, path, cell, internal = outlet_geometry()
    assert len(path) == 3
    sectors = receiving_cell(board, cell, internal)
    stocks, parent_gap = source_stocks()
    _, one, one_frame = sectors[1]
    response = clock(one, one[0], parent_gap / 10)
    calibration = (one_frame**2) @ response
    allowed = (one - one[0] >= 0.25) & (one - one[0] < parent_gap / 10)
    false_light = float(response[~allowed, 1].max())
    spectra = {}
    for n, (_, values, _) in sectors.items():
        starts = np.r_[0, np.flatnonzero(np.diff(values) > 1e-8) + 1, len(values)]
        spectra[str(n)] = [
            dict(energy=float(np.mean(values[a:b])), multiplicity=int(b - a))
            for a, b in zip(starts[:-1], starts[1:], strict=True)
        ]
    banks = {}
    for closed in (False, True):
        words, h = material(board, dimers, path, cell, internal, closed)
        values, v = eigh(h)
        groups = read_groups(words, sectors, cell)
        for (setting, arm), guide in stocks.items():
            if closed and setting != "native":
                continue
            name = "outlet_cut" if closed else setting
            initial = inject(guide, words, dimers, cell, one_frame[:, :1])
            spectral = v.T @ initial
            rows = []
            for t in WAITS:
                state = v @ (np.exp(-1j * t * values)[:, None] * spectral)
                rows.append(flight_read(state, groups, sectors, response))
            banks.setdefault(name, {})[arm] = rows
    ebit = int(path[0] == dimers[-1][0])
    works = {}
    for arm in ("pulse", "idle"):
        p = abs(stocks["native", arm]) ** 2
        # Each new blank contact has unlike-word energy 8.5 next to the
        # initially empty transit site. Both endpoints have no such hopping
        # coherence across the closed shutter at the switching appointment.
        gate_density = abs(one_frame[cell.index(path[-1]), 0]) ** 2
        works[arm] = dict(
            opening=float(8.5 * (p[:, :, :, ebit].sum() + gate_density)),
            removed_to_D=float(-65 + 32 * p[:, 0].sum() + p[:, 0, 1].sum()),
        )
    return dict(
        protocol_sha256=sha256((HERE / "OUTLET.md").read_bytes()).hexdigest(),
        source_wait=SOURCE_WAIT,
        waits=list(WAITS),
        dimers=dimers,
        path=list(path),
        cell=list(cell),
        inside_edges=internal,
        parent_gap=parent_gap,
        clock_labels=M,
        period=PERIOD,
        cutoff=parent_gap / 10,
        cell_spectra=spectra,
        native_ground_reference=float(one[0]),
        sharp_address_clock_other_light_heavy=calibration.round(13).tolist(),
        false_light_ceiling=false_light,
        switches=works,
        outcomes=list(OUTCOMES),
        panels=banks,
    )


def decode_panels(panels, waits):
    index = waits.index(32)
    for setting, arms in panels.items():
        for rows in arms.values():
            assert len(rows) == len(waits)
            for panel in rows:
                law = np.asarray(panel["core_by_outlet"])
                assert law.shape == (4, 5) and law.min() >= 0
                np.testing.assert_allclose(law.sum(), 1, atol=3e-12, rtol=0)
        pulse, idle = (
            np.asarray(arms[a][index]["core_by_outlet"]) for a in ("pulse", "idle")
        )
        print(
            setting,
            "pulse-minus-idle",
            np.round((pulse - idle).sum(0), 9),
            "pulse lower-core/light",
            pulse[1:3, 2].sum(),
            "idle light",
            idle[:, 2].sum(),
        )


def decode(data):
    decode_panels(data["panels"], data["waits"])
    print(
        "entry sharp clock (other/light/heavy)",
        data["sharp_address_clock_other_light_heavy"][
            data["cell"].index(data["path"][-1])
        ],
    )
    print("maximum clock false light", data["false_light_ceiling"])
    print(
        "switch work pulse-minus-idle",
        {
            k: v - data["switches"]["idle"][k]
            for k, v in data["switches"]["pulse"].items()
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    record = acquire()
    decode(record)
    if args.acquire:
        with args.acquire.open("x") as stream:
            json.dump(record, stream, allow_nan=False, separators=(",", ":"))
            stream.write("\n")
    else:
        compare_records(json.loads((HERE / "outlet.json").read_text()), record)
