"""Measure signed repairs, transported doubles and stationary material hands on the shared boards."""

import argparse
from collections import deque
from itertools import combinations, pairwise
import json
from pathlib import Path

import numpy as np

from model.executable_model import federation
from model.executable_model.carrier import seam_digest
from .candidate import KAPPA, NAME, native_law

HERE = Path(__file__).resolve().parent
SEED = 2026112131
RING = (0, 1, 7, 6, 10)
TAG = (5, 4)
SOURCE = (2, 3)
MATTER = 8
VACANCY = 9


def instruments():
    """Read the reduced sharp-sector transition tables from the actual 99-state native operations."""
    echo = native_law.echo()
    h0, h1 = native_law.local_term_parts()
    energies, v = np.linalg.eigh(h0 + KAPPA * h1)
    half = (v * np.exp(-0.5j * np.pi * energies)) @ v.T
    full = (v * np.exp(-1j * np.pi * energies)) @ v.T
    flags = np.arange(len(echo)) % 9 == 4
    basis = np.eye(len(echo))
    flag = np.zeros((3, 3), np.uint8)
    flip = np.zeros((3, 3))
    for a in (0, 2):
        for b in (0, 2):
            col, reverse = 3 * a + b, 3 * b + a
            restored = echo @ (echo[:, col] * (flags if a + b == 2 else ~flags))
            assert np.allclose(restored, basis[:, col], rtol=0, atol=1e-12)
            flag[a, b] = round(np.linalg.norm(echo[flags, col]) ** 2)
            assert abs(abs(full[reverse, col]) - 1) < 1e-12
            assert abs(np.sum(abs(half[[*{col, reverse}], col]) ** 2) - 1) < 1e-12
            if a != b:
                flip[a, b] = abs(half[reverse, col]) ** 2
    # For each prepared tag, keep the sharp public endpoints and test its retained
    # state using the same seam-11 projector, unwrite, endpoint-vacancy read and rearm.
    tips, tag_read = [], []
    for a, b in ((0, 0), (0, 2), (2, 0), (2, 2)):
        col = 3 * a + b
        committed = bool(flag[a, b])
        retained = echo[:, col] if committed else basis[:, col]
        visible = np.linalg.norm(retained[flags]) ** 2
        assert abs(visible - int(committed)) < 1e-12
        undone = echo @ retained if committed else retained
        assert np.allclose(undone, basis[:, col], rtol=0, atol=1e-12)
        left = int(np.argmax(abs(undone))) // 3
        read = int(flag[left, 0])
        if committed:
            assert np.allclose(echo @ undone, retained, rtol=0, atol=1e-12)
        tips.append((1, 1) if committed else (a, b))
        tag_read.append((2 * read - 1) if committed else 0)
    return flag, flip, np.asarray(tips, np.uint8), np.asarray(tag_read, np.int8)


def gateway_words(car):
    """Only use internal seams available on every cell, keeping the two tag ends protected."""
    internal = set(car.seams)
    local = {tuple(sorted((a, b))) for a, b in car.seams if a < 12 and b < 12}
    for c in range(car.cells):
        assert all((12 * c + a, 12 * c + b) in internal for a, b in local)
    needed = [TAG, (TAG[0], VACANCY), SOURCE, (SOURCE[1], VACANCY), (MATTER, VACANCY)]
    needed.extend(pairwise((*RING, RING[0])))
    assert all(tuple(sorted(e)) in local for e in needed)
    words = {}
    for gateway in sorted({p % 12 for edge in car.external for p in edge}):
        queue = deque([[MATTER]])
        seen = {MATTER, *TAG}
        while queue:
            path = queue.popleft()
            if path[-1] == gateway:
                swaps = list(pairwise(path))
                words[gateway] = swaps + swaps[-2::-1]
                break
            for p in range(12):
                if p not in seen and tuple(sorted((path[-1], p))) in local:
                    seen.add(p)
                    queue.append([*path, p])
        if gateway not in words:
            raise ValueError(
                f"Gateway {gateway} cannot be reached with the tag protected"
            )
    return words


def schedules(car):
    neighbors = [set() for _ in range(car.cells)]
    used = [set() for _ in neighbors]
    groups = []
    for a, b in car.external:
        x, y = a // 12, b // 12
        neighbors[x].add(y)
        neighbors[y].add(x)
        k = next(i for i in range(6) if i not in used[x] and i not in used[y])
        if k == len(groups):
            groups.append([])
        groups[k].append((x, y, a % 12, b % 12))
        used[x].add(k)
        used[y].add(k)
    patches = [
        (a, c, b)
        for c, ns in enumerate(neighbors)
        for a, b in combinations(sorted(ns), 2)
    ]
    assert len({tuple(sorted(path)) for path in patches}) == 3 * car.cells
    return [np.asarray(group).T for group in groups], np.asarray(patches, np.int32)


def exchange(loads, cells, a, b):
    loads[cells, a], loads[cells, b] = loads[cells, b].copy(), loads[cells, a].copy()


def route(loads, cells, gateways, words):
    for gate in np.unique(gateways):
        selected = cells[gateways == gate]
        for a, b in words[gate]:
            exchange(loads, selected, a, b)


def public_read(loads, tag_index, flag, sign):
    """Sharp, repeatable local readings under the enrolled probes; sign is an echo read lookup."""
    if np.any(loads[:, VACANCY]):
        raise ValueError("A local probe borrowed its reference without restoring it")
    ring = loads[:, RING]
    masks = flag[ring, np.roll(ring, -1, axis=1)] @ (2 ** np.arange(5))
    one_hot = 2 * np.eye(5, dtype=np.uint8)
    decode = np.full(32, -1, np.int8)
    decode[flag[one_hot, np.roll(one_hot, -1, axis=1)] @ (2 ** np.arange(5))] = (
        np.arange(5)
    )
    phase = decode[masks]
    assert np.all(phase >= 0)
    return np.column_stack(
        (sign[tag_index], flag[loads[:, MATTER], 0], flag[loads[:, 3], 0], phase)
    ).astype(np.int8)


def measure(level, services, seed, tables):
    car = federation(level)
    n = car.cells
    words = gateway_words(car)
    groups, patches = schedules(car)
    prep, coin, move = [
        np.random.default_rng(s) for s in np.random.SeedSequence([seed, level]).spawn(3)
    ]
    flag, flip, tips, sign = tables
    loads = np.zeros((n, 12), np.uint8)
    loads[:, [RING[0], SOURCE[0]]] = 2
    left, right = prep.random((2, n)) < 0.04
    tag_index = 2 * left.astype(np.uint8) + right
    loads[:, TAG] = tips[tag_index]
    loads[:, MATTER] = 2 * (prep.random(n) < 0.16)
    records = [public_read(loads, tag_index, flag, sign)]
    total = int(loads.sum())
    transfers = []
    for step in range(services):
        # Obtain the primary on the two-rail source with a real half contact, and
        # read rail 3 against 9. The stationary rule posts precisely a positive word.
        a, b = SOURCE
        changed = coin.random(n) < flip[loads[:, a], loads[:, b]]
        exchange(loads, changed, a, b)
        positive = flag[loads[:, b], loads[:, VACANCY]].astype(bool)
        for a, b in reversed(tuple(pairwise(RING))):
            exchange(loads, positive, a, b)
        # At most one external contact per cell. Both endpoint transpositions
        # return every borrowed rail before the local occupancy measurement.
        x, y, ga, gb = groups[step % len(groups)]
        route(loads, x, ga, words)
        route(loads, y, gb, words)
        before_x, before_y = loads[x, ga], loads[y, gb]
        swapped = move.random(len(x)) < flip[before_x, before_y]
        loads[x, ga] = np.where(swapped, before_y, before_x)
        loads[y, gb] = np.where(swapped, before_x, before_y)
        route(loads, x, ga, words)
        route(loads, y, gb, words)
        transfers.append(int(swapped.sum()))
        read = public_read(loads, tag_index, flag, sign)
        assert np.array_equal(read[:, 3], (records[-1][:, 3] + positive) % 5)
        assert int(loads.sum()) == total and np.array_equal(
            loads[:, TAG], tips[tag_index]
        )
        records.append(read)
    records = np.stack(records)
    patch_words = public_read(
        loads[patches.ravel()], tag_index[patches.ravel()], flag, sign
    ).reshape(len(patches), 3, 4)
    discrepancies = int(
        np.count_nonzero(np.any(patch_words != records[-1, patches], axis=2))
    )
    if discrepancies:
        raise ValueError("Overlapping observers recovered different checkpoint words")
    headers = dict(
        level=level,
        cells=n,
        seam_sha256=seam_digest(car),
        singletons=n,
        three_cell_patches=len(patches),
        checkpoint_patch_reads=n + 3 * len(patches),
        overlap_disagreements=discrepancies,
        positive_marks=int(records[1:, :, 2].sum()),
        quiet_marks=int(services * n - records[1:, :, 2].sum()),
        signed_negative=int(np.count_nonzero(records[-1, :, 0] == -1)),
        signed_positive=int(np.count_nonzero(records[-1, :, 0] == 1)),
        initial_doubles=int(records[0, :, 1].sum()),
        final_doubles=int(records[-1, :, 1].sum()),
        transfers=sum(transfers),
        final_hands=np.bincount(records[-1, :, 3], minlength=5).tolist(),
    )
    raw = dict(
        singleton_history=records,
        triple_support=patches,
        triple_checkpoint=patch_words,
        final_port_loads=loads,
        initial_tag_input=tag_index,
        transfers_per_service=np.asarray(transfers, np.int32),
    )
    return headers, raw


def run(levels=(0, 1, 2), *, seed=SEED, services=256, output=HERE):
    if (
        services < 1
        or seed < 0
        or not levels
        or len(set(levels)) != len(levels)
        or set(levels) - {0, 1, 2}
    ):
        raise ValueError(
            "Use positive service count, nonnegative seed and distinct shared federation levels 0–2"
        )
    tables = instruments()
    boards, archive = [], {}
    for level in levels:
        summary, raw = measure(level, services, seed, tables)
        boards.append(summary)
        archive.update({f"L{level}_{key}": value for key, value in raw.items()})
    result = dict(
        candidate=NAME,
        seed=seed,
        services=services,
        boards=boards,
        native_half_transfer_probability_02_20=float(tables[1][0, 2]),
        supplied="Independent double preparation on each tag end (0.04), mobile rail (0.16), separate fair sources and five-rail hands, cyclic matching schedule and protected readout windows.",
        record_fields=[
            "signed_repair_in_units_of_two",
            "mobile_double",
            "positive_primary",
            "five_rail_hand",
        ],
        shapes="singleton_history[service, cell, field], initial state at 0. triple_support[patch, site] names three cells with the owner of the hand in the middle. triple_checkpoint[patch, site, field] records all their terminal readings. initial_tag_input encodes two input occupancies as 2*(left occupied)+(right occupied).",
        scope="Stationary branch of the equipped clock during measured transport. Three-cell patches overlap. No distinct long-lived moving-object, intrinsic-abundance or gauge assignment is inferred.",
    )
    if output is not None:
        folder = Path(output)
        folder.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(folder / "patch_records.npz", **archive)
        (folder / "observations.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE)
    parser.add_argument("--services", type=int, default=256)
    parser.add_argument("--seed", type=int, default=SEED)
    options = parser.parse_args()
    print(json.dumps(run(**vars(options)), indent=2))
