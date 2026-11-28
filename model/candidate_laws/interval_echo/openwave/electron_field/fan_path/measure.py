"""Read full, energy-admitted and separated charge/ring records on physical paths."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import eigsh, expm_multiply

from .channel_bound import certify_channel
from .fixture import FanPath, archived_fixture
from .spectral import separated_reference, verified_inertia

HERE = Path(__file__).resolve().parent
TIMES = (27.0625, 30.8125)


def probability_range(effect):
    e = np.linalg.eigvalsh((effect + effect.conj().T) / 2)
    assert e[0] >= -5e-12 and e[-1] <= 1 + 5e-12
    return np.clip(e[[0, -1]], 0, 1).tolist()


def separated_effect(path, columns, distance):
    result = np.zeros((columns.shape[1],) * 2, dtype=complex)
    for block in path.blocks.values():
        axes = tuple(i for i in range(path.size) if abs(i - block.cell) >= distance)
        if not axes:
            continue
        part = columns[block.rows]
        tensor = part.reshape(*block.shape, columns.shape[1])
        # Each unoccupied ring has the unique normalized, uniform zero mode.
        vacuum = tensor.sum(axis=axes) / np.sqrt(
            np.prod([block.shape[i] for i in axes])
        )
        projected = vacuum.reshape(-1, columns.shape[1])
        result += part.conj().T @ part - projected.conj().T @ projected
    return result


def reference_read(kraus):
    flat = kraus.reshape(-1, 4)
    joint = flat.T @ flat.conj() / 2
    partial = joint.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    eig = np.linalg.eigvalsh(partial)
    success = np.trace(joint).real
    return {
        "reference_joint_probability": float(success),
        "reference_partial_transpose_minimum_per_launch": float(eig[0]),
        "conditional_reference_negativity": float(-np.minimum(eig, 0).sum() / success)
        if success > 1e-25
        else None,
    }


def terminal_law(path, columns):
    n = path.size - 1
    assert columns.shape == (len(path.basis), 2)
    kraus = path.decode(columns)
    vacuum = kraus.sum(axis=0) / np.sqrt(len(kraus))
    vacuum_effect = vacuum.conj().T @ vacuum
    accepted = np.sum(abs(columns) ** 2, axis=0)
    in_cell = np.sum(abs(columns[path.owner == n]) ** 2, axis=0)
    at_hub = np.sum(abs(columns[path.blocks[n, 0].rows]) ** 2, axis=0)
    effect = np.einsum("eos,eot->st", kraus.conj(), kraus)
    decoded = effect.diagonal().real
    categories = np.array(
        [1 - accepted, accepted - in_cell, in_cell - at_hub, at_hub - decoded, decoded]
    )
    np.testing.assert_allclose(categories.sum(axis=0), 1, rtol=0, atol=5e-12)
    assert categories.min() > -5e-12
    shape = path.blocks[n, 0].shape[:-1]
    word_tensor = kraus.reshape(*shape, 2, 2)
    separated_at_capture = []
    for d in range(1, path.size):
        axes = tuple(range(n - d + 1))
        quiet = word_tensor.sum(axis=axes) / np.sqrt(np.prod([shape[i] for i in axes]))
        joint = (np.vdot(kraus, kraus) - np.vdot(quiet, quiet)).real / 2
        rate = effect.trace().real / 2
        assert -5e-12 <= joint <= rate + 5e-12
        separated_at_capture.append(
            {
                "distance_at_least": d,
                "per_launch": float(np.clip(joint, 0, rate)),
                "conditional_on_received_charge_and_code": float(
                    np.clip(joint / rate, 0, 1)
                )
                if rate > 1e-25
                else None,
            }
        )
    return {
        "per_input_outcomes": np.clip(categories, 0, 1).T.tolist(),
        "joint_decoder_probability_range": probability_range(effect),
        **reference_read(kraus),
        "remote_positive_ring_and_decoded": separated_at_capture,
        "all_other_fans_vacuum_and_decoded": {
            **reference_read(vacuum[None, :, :]),
            "probability_range": probability_range(vacuum_effect),
            "receiver_accepted_but_spectator_vacuum_refused_per_input": np.clip(
                (effect - vacuum_effect).diagonal().real, 0, 1
            ).tolist(),
        },
        "separated_positive_ring_by_distance": [
            {
                "distance_at_least": d,
                "mean": float(np.clip(np.trace(e).real / 2, 0, 1)),
                "input_range": probability_range(e),
            }
            for d in range(1, path.size)
            for e in [separated_effect(path, columns, d)]
        ],
    }


def history(path, h, times):
    wave = path.source().astype(complex)
    last = 0
    result = []
    for t in times:
        dt = t - last
        if dt:
            op = (-1j * dt) * h
            wave = expm_multiply(op, wave, traceA=-1j * dt * h.diagonal().sum())
        np.testing.assert_allclose(wave.conj().T @ wave, np.eye(2), atol=7e-12, rtol=0)
        result.append({"appointment": t, **terminal_law(path, wave)})
        last = t
    return result


def spectral_read(path, h, threshold, bracket, sectors, times):
    cert = [verified_inertia(h, cut) for cut in bracket]
    count = cert[0]["count_below"]
    assert count == cert[1]["count_below"]
    rng = np.random.default_rng(923)
    e, v = eigsh(
        h,
        k=max(count + 2, 4),
        which="SA",
        v0=rng.standard_normal(h.shape[0]),
        tol=1e-13,
    )
    assert np.count_nonzero(e < threshold) == count
    residuals = np.linalg.norm(h @ v - v * e, axis=0)
    assert max(residuals) < 2e-10
    selected = v[:, :count]
    admitted_source = selected @ (selected.T @ path.source())
    records = {
        "threshold": threshold,
        "threshold_bracket": [str(cut) for cut in bracket],
        "inertia_certificates": cert,
        "energies_below_and_next": e[: count + 1].tolist(),
        "ritz_residual_norms": residuals[: count + 1].tolist(),
        "energy_admission_range": probability_range(
            admitted_source.T @ admitted_source
        ),
    }
    min_ground = min(
        Fraction(row["charged_ground_dyadic_bracket"][0]) for row in sectors
    )
    if count == 1:
        # Variational upper bound using an exactly representable integer trial vector.
        trial = np.rint(v[:, 0] * 2**12).astype(np.int64)
        twice_h = (h * 2).astype(np.int64)
        assert (twice_h != 2 * h).nnz == 0
        max_word = int(np.max(abs(trial)))
        assert len(trial) * max_word**2 * int(np.max(abs(twice_h).sum(axis=1))) < 2**62
        twice_energy = int(trial @ (twice_h @ trial))
        twice_norm = 2 * int(trial @ trial)
        upper = Fraction(twice_energy, twice_norm)
        separated_ceiling = upper - min_ground
        assert separated_ceiling < 1
        records.update(
            {
                "rank_one_time_independent_arrival": terminal_law(
                    path, admitted_source
                ),
                "eigenstate_separated_probability_by_distance": [
                    probability_range(separated_effect(path, selected, d))[0]
                    for d in range(1, path.size)
                ],
                "ground_energy_variational_upper": str(upper),
                "separated_probability_upper_from_energy": str(separated_ceiling),
            }
        )
    else:
        coefficients = selected.T @ path.source()
        records["energy_selected_appointments"] = [
            {
                "appointment": t,
                **terminal_law(
                    path, selected @ (np.exp(-1j * t * e[:count, None]) * coefficients)
                ),
            }
            for t in times
        ]
    if count == 2:
        records["two_level_channel_certificate"] = certify_channel(
            path, h, e[:count], selected, bracket[0], min_ground
        )
        environment_size = len(path.decode(selected))
        maps = [
            path.decode(selected[:, j : j + 1] @ coefficients[j : j + 1]).sum(axis=0)
            / np.sqrt(environment_size)
            for j in range(count)
        ]
        records["fine_energy_record_retained_vacuum_herald"] = reference_read(
            np.stack(maps)
        )
    return records


def conjugacy(path, negative):
    mask = sum(path.field.bits.values())
    perm = [negative.ids[p, tuple(mask ^ w for w in words)] for p, words in path.basis]
    assert (negative.parked[perm][:, perm] != path.parked).nnz == 0
    assert (negative.bridge_term[perm][:, perm] != path.bridge_term).nnz == 0
    gauges = []
    for c in (0, path.size - 1):
        a = path.source(c)
        b = negative.source(c)[perm]
        # Unoccupied terminal/source words may have a normalization involving sqrt(2).
        gauges.append(np.rint(a.T @ b).astype(int))
        np.testing.assert_array_equal(abs(gauges[-1]), np.eye(2))
        np.testing.assert_array_equal(a @ gauges[-1], b)
    return [g.tolist() for g in gauges]


def run():
    records = []
    with archived_fixture() as archive:
        for asymmetric in (False, True):
            for n in (2, 3, 4):
                path = FanPath(archive, n, 1, asymmetric=asymmetric)
                negative = FanPath(archive, n, -1, asymmetric=asymmetric)
                gauge = conjugacy(path, negative)
                threshold, bracket, sectors = separated_reference(path)
                times = sorted({0.0, *TIMES, *((n - 1) * t for t in TIMES)})
                cases = []
                for rate in (4, 2):
                    h = path.parked + (rate / 4) * path.bridge_term
                    cases.append(
                        {
                            "material_bridge_rate": rate,
                            "spectral": spectral_read(
                                path, h, threshold, bracket, sectors, times
                            ),
                            "appointments": history(path, h, times),
                        }
                    )
                records.append(
                    {
                        "board_cells": path.board.cells,
                        "fan_count": n,
                        "physical_fans": path.cells,
                        "physical_bridges": path.bridges,
                        "configuration_count": len(path.basis),
                        "source_and_receiver_conjugation": gauge,
                        "parked_reference": sectors,
                        "cases": cases,
                    }
                )
                print(
                    f"{path.board.cells} cells, {n} fans: "
                    + ", ".join(
                        f"rate {c['material_bridge_rate']} count {c['spectral']['inertia_certificates'][0]['count_below']}"
                        for c in cases
                    ),
                    flush=True,
                )
        census = path.census
    return {
        "protocol_sha256": {
            p: sha256((HERE / p).read_bytes()).hexdigest()
            for p in ("PROTOCOL.md", "TERMINAL_CONTROL.md")
        },
        "geometry_census_board_two_bridges_three_bridges": census,
        "outcome_order": [
            "energy_filter_refused",
            "accepted_charge_outside_receiver_cell",
            "receiver_cell_off_hub",
            "receiver_hub_doublet_refused",
            "receiver_hub_and_doublet",
        ],
        "paths": records,
    }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record", type=Path)
    args = cli.parse_args()
    data = run()
    if args.record:
        with args.record.open("x", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2, allow_nan=False)
            handle.write("\n")
    else:
        from runpy import run_path

        compare = run_path(str(HERE.parent / "replay.py"))["compare"]
        compare(
            json.loads(json.dumps(data)),
            json.loads((HERE / "results.json").read_text()),
            "fan_path",
        )
        print("Reproduced the path records.")


if __name__ == "__main__":
    main()
