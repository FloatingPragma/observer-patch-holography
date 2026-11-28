"""Copy a real wing-level record, then retain or extract the actual enrolled excitation."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from . import joined_outlet as joined
from .measure import HERE
from ..passive_decay.measure import compare_records, hamiltonians
from ..passive_decay.readers import Z, circuit

LEVELS = ("wing_0", "wing_1", "wing_2", "wing_3")
SHARP = ("E_empty", "E_two", *LEVELS)
EXTRACTED = ("return_1", "return_2", "return_3")
CASES = (
    ("known", "idle", 0),
    ("known", "D_phase", 0),
    ("decay", "idle", 32),
    ("decay", "core_pulse", 32),
)


def tracer_maps(words, terminal):
    values, v = eigh(hamiltonians(4)[1])
    frame = circuit(v, 33)
    size = 2 * len(words)
    maps = {tag: np.zeros((size, size)) for tag in (*SHARP, *EXTRACTED)}
    lookup = {frozenset(w): i for i, w in enumerate(words)}
    for d in (0, 1):
        for i, w in enumerate(words):
            occ = len(set(w) & set(terminal))
            if occ != 1:
                j = d * len(words) + i
                maps["E_empty" if occ == 0 else "E_two"][j, j] = 1
    outside = sorted(set(p for w in words for p in w) - set(terminal))
    for level, tag in enumerate(LEVELS):
        projector = np.outer(frame[:, level], frame[:, level])
        if level:
            flip = np.eye(4)
            flip[:, [0, level]] = flip[:, [level, 0]]
            reset = circuit(frame @ flip @ frame.T, 33) @ projector
        for q in outside:
            rows = [
                d * len(words) + lookup[frozenset((terminal[1 - e], q))]
                for d in (0, 1)
                for e in (0, 1)
            ]
            maps[tag][np.ix_(rows, rows)] = projector
            if level:
                maps[EXTRACTED[level - 1]][np.ix_(rows, rows)] = reset
    eye = np.eye(size)
    np.testing.assert_allclose(sum(maps[n] for n in SHARP), eye, atol=3e-13)
    for tag in SHARP:
        np.testing.assert_allclose(maps[tag].T @ maps[tag], maps[tag], atol=3e-13)
    for level, name in enumerate(EXTRACTED, 1):
        np.testing.assert_allclose(
            maps[name].T @ maps[name], maps[LEVELS[level]], atol=3e-13
        )
    maps["any_excited"] = sum(maps[n] for n in LEVELS[1:])
    maps["other"] = eye - maps["any_excited"]
    maps["untouched"] = eye
    return values - values[0], {k: sparse.csr_matrix(v) for k, v in maps.items()}


def remote_read(state, apparatus):
    n = 2 * len(apparatus.words)
    ranks = state.shape[1]
    amplitudes = np.einsum(
        "ab,acdr->dbcr", apparatus.energy_read, state.reshape(4, 2, n, ranks)
    )
    waves = amplitudes.reshape(n, -1)
    law = np.zeros((4, 5))
    for count, rows, v in apparatus.groups:
        projected = v.T @ waves[rows]
        lines = (abs(projected.reshape(len(rows), 4, -1)) ** 2).sum(2).T
        if count == 1:
            law[:, 1:4] += lines @ apparatus.response
        else:
            law[:, 0 if count == 0 else 4] += lines.sum(1)
    return law


def traced_flights(name, appointments, tags=None):
    apparatus = joined.apparatus(name)
    h = sparse.csr_matrix(apparatus.h)
    energies, v, stock, preparation = joined.prepare(
        apparatus.h, apparatus.initial, apparatus.cold_filter
    )
    n = 2 * len(apparatus.words)
    gap, maps = tracer_maps(apparatus.words, apparatus.dimers[-1])
    # The same finite free interval follows every read and reset.
    flights = [(v * np.exp(-1j * time * energies)) @ v.T for time in appointments]
    first = (v * np.exp(-32j * energies)) @ v.T
    outputs = [{} for _ in flights]
    if tags is not None:
        maps = {tag: maps[tag] for tag in tags}
    for case, arm, time in CASES:
        launch = stock.reshape(4, 2, n, -1)
        if arm == "core_pulse":
            launch = np.einsum("ba,acdr->bcdr", apparatus.collapse, launch)
        elif arm == "D_phase":
            launch = launch * np.repeat(Z.diagonal(), n // 2)[None, None, :, None]
        incoming = launch.reshape(len(v), -1)
        if time:
            incoming = first @ incoming
        rank = incoming.shape[1]
        blocks = incoming.reshape(8, n, rank).transpose(1, 0, 2).reshape(n, -1)
        rows = [output.setdefault(case, {}).setdefault(arm, {}) for output in outputs]
        for tag, action in maps.items():
            marked = (
                (action @ blocks)
                .reshape(n, 8, rank)
                .transpose(1, 0, 2)
                .reshape(len(v), rank)
            )
            energy = float(np.vdot(marked, h @ marked).real)
            for target, flight in zip(rows, flights, strict=True):
                law = remote_read(flight @ marked, apparatus)
                target[tag] = dict(
                    core_by_outlet=law.round(12).tolist(),
                    postread_weighted_energy=round(energy, 12),
                )
    return [
        dict(
            setting=name,
            preparation=preparation,
            isolated_wing_gaps=gap.tolist(),
            cases=out,
        )
        for out in outputs
    ]


def acquire_setting(name):
    return traced_flights(name, [joined.TICKS[1]])[0]


def all_outputs(rows):
    return {name: np.array(v["core_by_outlet"]) for name, v in rows.items()}


def contrasts(pulse, idle):
    def quantities(rows):
        r = all_outputs(rows)
        retained = sum(r[n] for n in LEVELS[1:])
        extracted = sum(r[n] for n in EXTRACTED)
        effect = retained - extracted
        work = sum(rows[n]["postread_weighted_energy"] for n in EXTRACTED) - sum(
            rows[n]["postread_weighted_energy"] for n in LEVELS[1:]
        )
        sharp = sum(r[n] for n in SHARP)
        coarse = r["any_excited"] + r["other"]
        np.testing.assert_allclose(sharp.sum(), r["untouched"].sum(), atol=2e-11)
        np.testing.assert_allclose(coarse.sum(), r["untouched"].sum(), atol=2e-11)
        np.testing.assert_allclose(retained.sum(), extracted.sum(), atol=2e-11)
        return np.asarray(
            [
                retained.sum(),
                r["untouched"][:, 2].sum(),
                coarse[:, 2].sum(),
                sharp[:, 2].sum(),
                effect[:, 2].sum(),
                effect[1:3, 2].sum(),
                work,
            ]
        )

    return quantities(pulse) - quantities(idle)


def decode(record):
    for name, s in record["settings"].items():
        for case, arms in s["cases"].items():
            p = arms["D_phase" if case == "known" else "core_pulse"]
            result = contrasts(p, arms["idle"])
            if name == "exit_cut":
                assert abs(result[4]) < 1e-10
            print(
                name,
                case,
                "tag, free/coarse/sharp soft, retain-minus-extract soft and joint, reset work",
                np.round(result, 7),
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--setting", choices=joined.SETTINGS)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    if args.setting:
        if not args.acquire:
            parser.error("A single-setting acquisition needs an output")
        result = acquire_setting(args.setting)
    else:
        result = dict(
            protocol_sha256=sha256((HERE / "WING_TRACER.md").read_bytes()).hexdigest(),
            settings={name: acquire_setting(name) for name in joined.SETTINGS},
        )
        decode(result)
    if args.acquire:
        with args.acquire.open("x") as f:
            json.dump(result, f, allow_nan=False, separators=(",", ":"))
            f.write("\n")
    else:
        compare_records(json.loads((HERE / "wing_tracer.json").read_text()), result)
