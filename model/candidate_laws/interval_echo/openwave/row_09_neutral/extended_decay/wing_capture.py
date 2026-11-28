"""Move a positive material wave coherently onto a distinct physical wing."""

import argparse
from collections import deque
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh, expm
from . import joined_outlet as joined, open_ports, wing_tracer as tracer
from .material import edge, route
from .measure import HERE
from ..passive_decay.measure import compare_records, hamiltonians, physical
from ..passive_decay.readers import circuit
from model.candidate_laws.interval_echo import candidate


def capture_hardware():
    b, dimers, cell, available, _, qualified, _ = open_ports.aperture()
    guide_heads = route((3,), crossing=True)[2]
    busy = (
        {p for ds in (*dimers, *guide_heads) for p in ds} | set(cell) | set(available)
    )
    source_owners = {p // 12 for p in busy}
    for owner in range(b.cells):
        if owner in source_owners:
            continue
        interior = [p for p in range(12 * owner, 12 * owner + 12) if b.degree[p] == 5]
        for p, q, r, s in permutations(interior, 4):
            if (
                q not in b.neighbours[p]
                or r not in b.neighbours[p]
                or s not in b.neighbours[p]
                or s not in b.neighbours[r]
            ):
                continue
            pins = sorted(set(b.neighbours[q]) - {p, q, r, s} - busy)
            if pins:
                break
        else:
            continue
        break
    else:
        raise ValueError("No disjoint copy wing")
    copy = ((p, q), (r, s))
    heads = {edge(p, r): 33, edge(p, s): 32, edge(q, pins[0]): 32}
    wing_h = hamiltonians(4)[1]
    potential = physical(b, copy, heads)
    np.testing.assert_array_equal(potential - potential[0], wing_h.diagonal())
    # A single routed whole-load swap may use only physical blank contacts.
    # All previous guide contacts touching its five dimers were written unless
    # they are its blank dimer or part of the newly qualified open approach.
    guide_sites = set(sum(dimers, ()))
    blank = {edge(*d) for d in dimers} | set(qualified)
    blocked = (
        set(guide_heads)
        | set(heads)
        | {
            edge(a, z)
            for a in guide_sites
            for z in b.neighbours[a]
            if edge(a, z) not in blank
        }
    )

    def transposition(a, z, labels):
        paths = deque([(a,)])
        seen = {a}
        while paths:
            path = paths.popleft()
            if path[-1] == z:
                break
            for j in sorted(b.neighbours[path[-1]]):
                if j not in seen and edge(path[-1], j) not in blocked:
                    seen.add(j)
                    paths.append(path + (j,))
        else:
            raise ValueError((a, z))
        pairs = list(zip(path, path[1:]))
        original = list(labels)
        for x, y in pairs + pairs[-2::-1]:
            labels[x], labels[y] = labels[y], labels[x]
        original[a], original[z] = original[z], original[a]
        assert labels == original
        return list(path)

    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    native = exchange + candidate.KAPPA * contrast
    swap9 = np.eye(9)[[law.SRC.index((y, x)) for x, y in law.SRC]]
    np.testing.assert_allclose(
        expm(-1j * np.pi * native)[:, :9],
        np.vstack((swap9, np.zeros((90, 9)))),
        atol=3e-13,
    )
    swap = circuit(np.eye(4)[[0, 2, 1, 3]], 33)
    np.testing.assert_allclose(swap, np.eye(4)[[0, 2, 1, 3]], atol=1e-13)
    # Hardware comparison is made on the logical one-load registers at the
    # instants of the two logical SWAP gates. Routed intermediates need only be ternary.
    d, e = dimers[-2:]
    labels = list(range(b.n))
    routes = []
    for a, z in zip(e, copy[1], strict=True):
        routes.append(transposition(a, z, labels))
    for a, z in zip(d, e, strict=True):
        labels[a], labels[z] = labels[z], labels[a]
    for a, z in zip(e, copy[0], strict=True):
        routes.append(transposition(a, z, labels))
    for a, z in zip(d, e, strict=True):
        labels[a], labels[z] = labels[z], labels[a]
    target = list(range(b.n))
    copy_sites = copy[0] + copy[1]
    for a, z in zip(d + e, copy_sites, strict=True):
        target[a], target[z] = target[z], target[a]
    assert target == labels and len(set(copy_sites)) == 4
    _, basis = eigh(wing_h)
    frame = circuit(basis, 33)
    raw = np.diag([-1.0, 1.0, 1.0, 1.0])
    phase = circuit(raw, 33)
    # Ramsey control over this real reflection writes one bit, preserving every
    # relative phase within the three-dimensional excited subspace.
    reflected = frame @ phase @ frame.T
    coarse = (np.eye(4) + reflected) / 2
    np.testing.assert_allclose(coarse, frame[:, 1:] @ frame[:, 1:].T, atol=2e-13)
    return dict(
        guide_wing=[list(d), list(e)],
        copy=[list(pair) for pair in copy],
        copy_contacts=[[a, z, w] for (a, z), w in heads.items()],
        blank_transposition_paths=routes,
        logical_swaps=2,
        full_blank_swaps=sum(2 * len(path) - 3 for path in routes),
    )


def acquire():
    hardware = capture_hardware()
    settings = {}
    ticks = [162, 164]
    for name in joined.SETTINGS:
        banks = tracer.traced_flights(
            name,
            np.array(ticks) * 2 * np.pi / 256,
            tags=("any_excited", *tracer.EXTRACTED, "untouched"),
        )
        settings[name] = dict(zip(map(str, ticks), banks, strict=True))
    return dict(
        protocol_sha256=sha256((HERE / "WING_CAPTURE.md").read_bytes()).hexdigest(),
        hardware=hardware,
        settings=settings,
    )


def contrast(pulse, idle, gaps):
    def totals(rows):
        laws = tracer.all_outputs(rows)
        captured = [laws[n] for n in tracer.EXTRACTED]
        retain = laws["any_excited"]
        remove = sum(captured)
        np.testing.assert_allclose(retain.sum(), remove.sum(), atol=2e-11)
        impact = retain - remove
        parked = sum(e * r.sum() for e, r in zip(gaps[1:], captured, strict=True))
        return np.array(
            [
                retain.sum(),
                parked,
                retain[:, 2].sum(),
                remove[:, 2].sum(),
                impact[:, 2].sum(),
                impact[1:3, 2].sum(),
            ]
        )

    return totals(pulse) - totals(idle)


def decode(data):
    for setting, pages in data["settings"].items():
        for tick, page in pages.items():
            for c, arms in page["cases"].items():
                n = "D_phase" if c == "known" else "core_pulse"
                row = contrast(arms[n], arms["idle"], page["isolated_wing_gaps"])
                if setting == "exit_cut":
                    assert abs(row[4]) < 1e-10
                print(
                    setting,
                    tick,
                    c,
                    "tags, captured energy, retain/capture soft, impact and lower-core impact",
                    np.round(row, 7),
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    result = acquire()
    decode(result)
    if args.acquire:
        with args.acquire.open("x") as out:
            json.dump(result, out, allow_nan=False, separators=(",", ":"))
            out.write("\n")
    else:
        compare_records(json.loads((HERE / "wing_capture.json").read_text()), result)
