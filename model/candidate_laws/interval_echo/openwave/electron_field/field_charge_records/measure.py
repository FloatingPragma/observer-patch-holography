"""Compare independent material and enclosing-head records on the archived gateway."""

import argparse
from collections import Counter, defaultdict, deque
from contextlib import contextmanager
from fractions import Fraction
from hashlib import sha256
from itertools import pairwise
import json
from pathlib import Path
from runpy import run_path
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import numpy as np
from scipy.linalg import eigh


HERE = Path(__file__).resolve().parent
FIELD = HERE.parent
SEARCH = "experiments/interval_echo_collective_spin/search.py"
ARCHIVE_RESULT = "experiments/interval_echo_collective_spin/results.json"
TOL = 3e-12
REGIONS = {
    "left_fan": tuple(range(6)),
    "right_fan": tuple(range(6, 12)),
    "receiving_gateway": (6,),
    "receiving_rim": tuple(range(7, 12)),
    "complete_active_fixture": tuple(range(12)),
}


@contextmanager
def archived_fixture():
    manifest = dict(
        line.split("  ", maxsplit=1)[::-1]
        for line in (FIELD / "records.sha256").read_text().splitlines()
    )
    original_path = sys.path[:]
    with TemporaryDirectory() as temp, ZipFile(FIELD / "records.zip") as archive:
        root = Path(temp)
        assert sorted(archive.namelist()) == sorted(manifest)
        assert all(
            sha256(archive.read(name)).hexdigest() == digest
            for name, digest in manifest.items()
        )
        repo = FIELD.parents[4]
        (root / repo.name).symlink_to(repo, target_is_directory=True)
        archive.extractall(root)
        sys.path[:0] = [str((root / SEARCH).parent), str(repo)]
        try:
            from search import JoinedSpace, fan, vacuum_code

            data = json.loads(archive.read(ARCHIVE_RESULT))
            best = data["headline"]["best_worst_qubit_success"]
            yield fan(4), JoinedSpace, vacuum_code, best, manifest
        finally:
            sys.path[:] = original_path


def head_records(positions, words, edges):
    """Read the selected endpoint of each seam and the load at each port."""
    ends = np.asarray(edges, dtype=int)
    values = np.asarray(words, dtype=np.int64)
    heads = np.where(
        (values[:, None] >> np.arange(len(edges))) & 1, ends[:, 0], ends[:, 1]
    )
    at_port = np.arange(12)[None, :] == np.asarray(positions)[:, None]
    return heads, at_port


def record_calibration(space, sign):
    """Compute separate load and incoming-head reads for each enclosing cut."""
    positions, words = zip(*space.basis, strict=True)
    heads, at_port = head_records(positions, words, space.edges)
    q = sign * at_port.astype(int)
    incoming = (heads[:, :, None] == np.arange(12)).sum(axis=1)
    initial = space.portal[0][0, 0]
    background = incoming[initial] - q[initial]
    np.testing.assert_array_equal(incoming - q, np.broadcast_to(background, q.shape))
    assert int(background.sum()) == len(space.edges) - sign

    # Each row is a geometric cut, each column a possible sharp detector event.
    masks = ((np.arange(1 << 12)[:, None] >> np.arange(12)) & 1).astype(bool)
    inward = np.zeros((len(masks), len(positions)), dtype=int)
    ends = np.asarray(space.edges)
    for seam, (u, v) in enumerate(ends):
        cut = masks[:, u] ^ masks[:, v]
        inward[cut] += masks[cut][:, heads[:, seam]]
    delta_head = inward - inward[:, [initial]]
    load = masks @ q.T
    delta_load = load - load[:, [initial]]
    np.testing.assert_array_equal(delta_head, delta_load)
    region_ids = [sum(1 << p for p in ports) for ports in REGIONS.values()]
    np.testing.assert_array_equal(delta_head[-1], 0)
    np.testing.assert_array_equal(load[-1], sign)
    # Partition identity includes heads on opposite sides of inner seams, which
    # contribute one constant to the unsubtracted sum.
    np.testing.assert_array_equal(
        load[region_ids[2]] + load[region_ids[3]], load[region_ids[1]]
    )
    np.testing.assert_array_equal(
        delta_head[region_ids[2]] + delta_head[region_ids[3]], delta_head[region_ids[1]]
    )
    return {
        "head_by_seam": heads,
        "signed_loads_by_port": q,
        "signed_loads_by_region": load[region_ids].T,
        "signed_inward_head_changes_by_region": delta_head[region_ids].T,
        "reference_inward_head_counts": inward[region_ids, initial].tolist(),
        "reference_signed_loads": load[region_ids, initial].tolist(),
        "incoming_head_stock": len(space.edges),
        "neutral_background_summed_heads": int(background.sum()),
        "checked_cut_count": len(masks),
    }


def verify_steps(space, field, sign, degrees, h, records):
    """Assemble the generator independently from each dated contact rule."""
    heads = records["head_by_seam"]
    loads = records["signed_loads_by_port"]
    edges = space.edges
    rings = [tuple(x + off for x in ring) for off in (0, 6) for ring in field.rings]
    ring_pairs = [list(pairwise((*ring, ring[0]))) for ring in rings]
    twice = np.zeros(h.shape, dtype=int)
    counts = Counter()
    for i, (x, word) in enumerate(space.basis):
        twice[i, i] = 16 * (6 - degrees[x])
        for seam, (u, v) in enumerate(edges):
            if x not in (u, v) or ((heads[i, seam] == x) != (sign == 1)):
                continue
            dest = v if x == u else u
            j = space.ids[dest, word ^ (1 << seam)]
            change = np.bincount(heads[j], minlength=12) - np.bincount(
                heads[i], minlength=12
            )
            np.testing.assert_array_equal(change, loads[j] - loads[i])
            assert np.count_nonzero(heads[j] != heads[i]) == 1
            twice[i, i] += 8
            twice[j, i] -= 8
            counts["material_steps"] += 1
            if (u, v) == (0, 6):
                counts["gateway_material_steps"] += 1
        for path in ring_pairs:
            seams = [edges.index(tuple(sorted(edge))) for edge in path]
            orientation = [heads[i, seam] == u for (u, _), seam in zip(path, seams)]
            if any(orientation) and not all(orientation):
                continue
            j = space.ids[x, word ^ sum(1 << seam for seam in seams)]
            np.testing.assert_array_equal(
                np.bincount(heads[j], minlength=12), np.bincount(heads[i], minlength=12)
            )
            np.testing.assert_array_equal(loads[j], loads[i])
            twice[i, i] += 1
            twice[j, i] -= 1
            counts["ring_steps"] += 1
    np.testing.assert_array_equal(twice, 2 * h)
    initial = set(space.portal[0].ravel())
    reachable = set(initial)
    queue = deque(initial)
    while queue:
        index = queue.popleft()
        for j in np.flatnonzero(twice[index] < 0):
            if j not in reachable:
                reachable.add(j)
                queue.append(j)
    assert len(reachable) == len(h)
    return counts


def controlled_half_contact(space, source, sign, records):
    """Use Gaussian-integer amplitudes for the two half-contact comparisons."""
    bridge_index = space.edges.index((0, 6))
    head = records["head_by_seam"][:, bridge_index]
    a, b = [side.ravel() for side in space.portal]
    assert np.all(head[a] != head[b]) and len(set(head[a])) == 1
    np.testing.assert_array_equal(
        records["signed_loads_by_region"][b, 1]
        - records["signed_loads_by_region"][a, 1],
        sign,
    )
    basis = [(0, 0), (0, 1), (1, 0), (1, 1)]
    half_numerator = np.array([[1, 1j], [1j, 1]])
    joint_numerator = np.zeros((len(basis), len(a), 2), dtype=complex)
    separated_numerator = np.zeros_like(joint_numerator)
    input_numerator = source[a] * 4
    np.testing.assert_array_equal(np.abs(input_numerator), 1)
    joint_numerator[[0, 3]] = half_numerator[:, :1, None] * input_numerator
    separated_numerator[:] = (
        np.kron(half_numerator, half_numerator)[:, :1, None] * input_numerator
    )

    def probabilities(value, denominator):
        counts = (value.real**2 + value.imag**2).sum(axis=1)
        np.testing.assert_array_equal(counts, np.rint(counts))
        return [
            [Fraction(int(weight), denominator) for weight in row] for row in counts
        ]

    native = probabilities(joint_numerator, 32)
    separated = probabilities(separated_numerator, 64)
    assert native == [[Fraction(1, 2)] * 2, [0] * 2, [0] * 2, [Fraction(1, 2)] * 2]
    assert separated == [[Fraction(1, 4)] * 2] * 4
    output = []
    for (placed, flipped), coupled, independent in zip(
        basis, native, separated, strict=True
    ):
        output.append(
            {
                "receiver_signed_load": placed * sign,
                "receiver_signed_inward_head_change": flipped * sign,
                "coupled_probability_by_input": [str(x) for x in coupled],
                "independent_probability_by_input": [str(x) for x in independent],
            }
        )
    return output


def release_probabilities(h, source, receiver, records, times):
    energies, vectors = eigh(h)
    amplitudes = vectors.T @ source
    groups = defaultdict(list)
    for i, (loads, heads) in enumerate(
        zip(
            records["signed_loads_by_region"],
            records["signed_inward_head_changes_by_region"],
            strict=True,
        )
    ):
        groups[tuple(loads), tuple(heads)].append(i)
    arrival = records["signed_loads_by_region"][:, 1] != 0
    np.testing.assert_array_equal(
        arrival, records["signed_inward_head_changes_by_region"][:, 1] != 0
    )
    # The whole-cell boundary read is a scalar on every accepted receiver mode.
    np.testing.assert_array_equal(receiver[~arrival], 0)
    np.testing.assert_array_equal(receiver.T @ receiver, np.eye(receiver.shape[1]))
    rows = []
    for time in times:
        evolved = vectors @ (np.exp(-1j * time * energies[:, None]) * amplitudes)
        probabilities = abs(evolved) ** 2
        np.testing.assert_allclose(probabilities.sum(axis=0), 1, atol=TOL, rtol=0)
        joint_law = [
            {
                "signed_loads": list(map(int, loads)),
                "signed_inward_head_changes": list(map(int, heads)),
                "probability_by_input": probabilities[indices].sum(axis=0).tolist(),
            }
            for (loads, heads), indices in sorted(groups.items())
        ]
        whole = probabilities[arrival].sum(axis=0)
        receiver_doublet = (abs(receiver.T @ evolved) ** 2).sum(axis=0)
        assert np.all(receiver_doublet <= whole + TOL)
        rows.append(
            {
                "archived_exposure": time,
                "all_starts_joint_load_and_flux": joint_law,
                "right_cell_load_and_signed_flux_probability_by_input": whole.tolist(),
                "right_cell_flux_and_gateway_load_and_receiver_doublet_probability_by_input": receiver_doublet.tolist(),
            }
        )
    return rows


def sector(field, joined_space, vacuum_code, heads, degrees, sign, times):
    space = joined_space(field, (heads, heads), sign)
    degree = [6, *degrees[0], 6, *degrees[1]]
    h = space.hamiltonian(degree)
    records = record_calibration(space, sign)
    counts = verify_steps(space, field, sign, degree, h, records)
    vacuum, code = vacuum_code(field, heads)
    source = np.zeros((len(h), 2))
    receiver = np.zeros((len(h), 4, 2))
    source[space.portal[0]] = np.einsum("ik,j->ijk", code, vacuum)
    receiver[space.portal[1]] = np.einsum("ia,jk->ijak", np.eye(4), code)
    np.testing.assert_array_equal(source.T @ source, np.eye(2))
    receiver = receiver.reshape(len(h), 8)

    parked = h.copy()
    for a, b in zip(*(side.ravel() for side in space.portal), strict=True):
        assert space.material_step(space.basis[a], (0, 6)) == space.basis[b]
        assert parked[a, b] == parked[b, a] == -4
        parked[a, b] = parked[b, a] = 0
        parked[a, a] -= 4
        parked[b, b] -= 4
    right = records["signed_loads_by_region"][:, 1] != 0
    np.testing.assert_array_equal(parked[np.ix_(right, ~right)], 0)
    np.testing.assert_array_equal(source[right], 0)

    return {
        "sign": sign,
        "word_count": len(h),
        "local_step_counts_directed": dict(counts),
        "cut_count": records["checked_cut_count"],
        "active_head_stock": records["incoming_head_stock"],
        "putative_neutral_background_head_total": records[
            "neutral_background_summed_heads"
        ],
        "reference_region_loads": records["reference_signed_loads"],
        "reference_region_inward_heads": records["reference_inward_head_counts"],
        "half_contact_exact": controlled_half_contact(space, source, sign, records),
        "parked_right_charge_and_flux_probability_both_inputs_all_times": [0, 0],
        "release": release_probabilities(h, source, receiver, records, times),
    }


def run():
    with archived_fixture() as (field, joined, source_code, best, manifest):
        release = best["release"]
        assert release["sign"] == 1
        assert release["head_sectors"][0] == release["head_sectors"][1]
        heads = tuple(release["head_sectors"][0])
        conjugate = tuple(
            sum(i in edge for edge in field.edges) - h for i, h in enumerate(heads)
        )
        times = [release["best_receiver_only_time"], release["best_worst_time"]]
        assert times == [27.0625, 30.8125]
        sectors = [
            sector(field, joined, source_code, labels, best["degrees"], sign, times)
            for sign, labels in ((1, heads), (-1, conjugate))
        ]
        for time_record, conjugate_record in zip(
            *(section["release"] for section in sectors), strict=True
        ):
            key = "right_cell_flux_and_gateway_load_and_receiver_doublet_probability_by_input"
            np.testing.assert_allclose(
                time_record[key], conjugate_record[key], atol=TOL
            )
        np.testing.assert_allclose(
            np.mean(sectors[0]["release"][0][key]),
            release["receiver_only_mean_joint_success_at_that_time"],
            atol=TOL,
            rtol=0,
        )
        return {
            "design_sha256": sha256((HERE / "DESIGN.md").read_bytes()).hexdigest(),
            "archived_search_sha256": manifest[SEARCH],
            "archived_result_sha256": manifest[ARCHIVE_RESULT],
            "physical_ports_by_cell": best["actual_ports"],
            "port_regions_local_order": REGIONS,
            "sectors": sectors,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as out:
            out.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(FIELD / "replay.py"))["compare"]
        # Match the JSON's sequence types before comparing the replay.
        compare(
            json.loads(json.dumps(result)),
            json.loads((HERE / "results.json").read_text()),
            "field_charge_records",
        )
        for item in result["sectors"]:
            print(
                f"q={item['sign']:+}, {item['word_count']} admitted words, "
                f"{item['cut_count']} cut identities, equal-marginal mismatch 0 vs 1/2"
            )


if __name__ == "__main__":
    main()
