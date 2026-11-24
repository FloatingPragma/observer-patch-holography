"""Local contrast, regional voltage and written-head transport on prepared belt walls."""

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo import candidate, currents
from .. import measure
from ...electron_field.flux_string import search
from . import belt


def contact():
    x, c = law.local_term_parts()
    pulse = expm(-1j * np.pi / 64 * (x + candidate.KAPPA * c))
    y = {(n, r): float((1 + pulse[i, i].imag) / 2) for i, (n, r) in enumerate(law.SRC)}
    return {f"{n},{r}": v for (n, r), v in y.items()}, np.array(
        [y[n, 2] - y[n, 0] for n in range(3)]
    )


def region(sites, slope):
    length = 10
    _, _, centre = belt.wall(sites, length)
    n = np.isin(np.arange(length), sites).astype(int)
    unit = np.zeros(length, int)
    unit[int(np.floor(centre))] = 1
    results = []
    for radius in (2, 3):
        centres = np.flatnonzero(
            abs(belt.signed_distance(np.arange(length), centre, length)) <= radius
        )
        local = []
        charges = Counter()
        for j in centres:
            address = (j + np.array([-1, 0, 1])) % length
            # Same unwrapped coordinate for the three different bulk phases at this cut.
            x = j + np.array([-1, 0, 1])
            ref_signals = [slope[(x % 3 == m).astype(int)].sum() for m in range(3)]
            local.append(
                (
                    slope[n[address]].sum() - np.mean(ref_signals),
                    (slope[unit[address]] - slope[0]).sum(),
                )
            )
        # Lift the selected contiguous window about its own centre, including wraparound.
        unwrapped = np.rint(
            belt.signed_distance(centres, centre, length) + centre
        ).astype(int)
        for m in range(3):
            charges[int(n[centres].sum() - (unwrapped % 3 == m).sum())] += 1 / 3
        total, calibration = np.sum(local, axis=0)
        q = total / calibration
        results.append(
            dict(
                radius=radius,
                cut=list(map(int, centres)),
                signed_contrast_sum=total,
                unit_contrast_sum=calibration,
                relative_excess=q,
                regional_voltage_distribution={
                    str(k): v for k, v in sorted(charges.items())
                },
                regional_voltage=measure.voltage_ports(charges),
                sharp_fraction_voltage=measure.voltage_ports({q: 1}),
                regional_full_turn_bright_x=(
                    1 + sum(p * np.cos(2 * np.pi * k) for k, p in charges.items())
                )
                / 2,
                sharp_fraction_full_turn_bright_x=(1 + np.cos(2 * np.pi * q)) / 2,
            )
        )
    return results


def written(launch, order, table, patch):
    s = tuple(launch)
    words = patch.field_component()
    steps = []
    for _ in range(10):
        centre = belt.wall(s, 10)[2]
        (dest,) = [
            d
            for d in table
            if len(set(d) - set(s)) == 1
            and belt.signed_distance(belt.wall(d, 10)[2], centre, 10) == 3
            and abs(
                belt.signed_distance(
                    next(iter(set(d) - set(s))), next(iter(set(s) - set(d))), 10
                )
            )
            <= 2
        ]
        (a,) = set(s) - set(dest)
        (b,) = set(dest) - set(s)
        dist = int(belt.signed_distance(b, a, 10))
        steps.append((a, b, dist, tuple(sorted((order[a], order[b])))))
        s = dest
    assert s == launch and len({e for *_, e in steps}) == 10
    # A nonseparating cross-section between unwrapped sites 9 and 10.
    aligned = patch.seed
    for a, b, dist, e in steps:
        bit = patch.bit[e]
        aligned = (aligned | bit) if order[a] == e[0] else (aligned & ~bit)

    cut_heads = []
    for (u, v), bit in patch.bit.items():
        if u in order and v in order:
            i, j = order.index(u), order.index(v)
            crossings = (i + int(belt.signed_distance(j, i, 10))) // 10
            if crossings:
                cut_heads.append((bit, crossings))

    def tour(seed, write):
        head = seed
        pos = set(launch)
        completed = 0
        material_cut = 0
        for a, b, d, e in steps:
            if a not in pos or b in pos:
                break
            out, newhead = (
                patch.move(order[a], head, e, 1) if write else (order[b], head)
            )
            if out != order[b]:
                break
            crossings = (a + d) // 10
            material_cut += crossings
            completed += 1
            pos = (pos - {a}) | {b}
            head = newhead
        head_cut = sum(
            crossings * (int(bool(seed & bit)) - int(bool(head & bit)))
            for bit, crossings in cut_heads
        )
        return completed, material_cut, head_cut, head

    counts = Counter(tour(w, True)[0] for w in words)
    count, matter, field, end = tour(aligned, True)
    reverse = end
    for a, b, _, e in reversed(steps):
        at, reverse = patch.move(order[b], reverse, e, 1)
        assert at == order[a]
    assert reverse == aligned and count == 10
    no_write = tour(aligned, False)
    return dict(
        scheduled_steps=[
            dict(old=a, new=b, material_step=d, physical_seam=e) for a, b, d, e in steps
        ],
        completed_prefix_counts={str(i): counts[i] for i in range(11)},
        archived_field_words=len(words),
        aligned_head=aligned,
        terminal_head=end,
        changed_heads=(end ^ aligned).bit_count(),
        reverse_restores_head=True,
        winding_laps=3,
        unit_material_crossings=matter,
        independent_head_crossings=field,
        unwritten_control=dict(
            completed=no_write[0],
            material_crossings=no_write[1],
            head_crossings=no_write[2],
        ),
    )


def run():
    cal, slopes = contact()
    patch = search.Patch()
    order = next(iter(belt.belts().values()))
    data = {}
    for n in (3, 4):
        sector = BinarySector(belt.ring(10), n)
        table = [sector.states[i] for i in belt.selected(sector)]
        (launch,) = [s for s in table if belt.wall(s, 10)[0] == 0]
        data[str(n)] = dict(
            launch=launch,
            contrast_and_voltage=region(launch, slopes),
            written_tour=written(launch, order, table, patch),
        )
    a, b = (set(data[str(n)]["launch"]) for n in (3, 4))
    return currents._portable(
        dict(
            calibration_bright_y=cal,
            two_reference_signal=slopes,
            cases=data,
            unit_attachment=dict(removed=sorted(a - b), added=sorted(b - a)),
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    data = run()
    encoded = json.dumps(data, indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as output:
        output.write(encoded)
