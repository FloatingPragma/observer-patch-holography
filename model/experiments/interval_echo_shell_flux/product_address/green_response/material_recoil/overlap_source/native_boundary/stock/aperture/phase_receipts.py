"""Source-calibrated earliest-receipt time for two one-use ordinary signals."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import expm

from . import thermal_radar as survey
from . import source_dwell as dwell

contact = dwell.contact
HERE = contact.HERE
SHOTS = 2**25


def written_receiver():
    # Start receiver at its port. Equal paths for an absent sender must produce
    # a dark port even when the written head has its own vacuum evolution.
    blank = contact.stock.kernels()[2][:9, :9]
    split = expm(-1j * np.pi * blank[np.ix_([3, 1], [3, 1])] / 8)
    probs = []
    for mag, (base, extra) in contact.blocks()["full"].items():
        ref = expm(-1j * dwell.DWELL * base / 4)
        for present, block in enumerate((base, base + extra)):
            arm = expm(-1j * dwell.DWELL * block / 4)
            dark = split[0, 0] ** 2 * arm + split[0, 1] * split[1, 0] * ref
            p = (abs(dark) ** 2).sum(0)
            if present and mag == 2:
                np.testing.assert_allclose(p, p.mean(), atol=1e-13)
                probs.append(float(p.mean()))
            else:
                assert p.max() < 1e-24
    return probs[0]


def free_probe(board, edges):
    """Reduce each native column on an evacuated chart with one test load."""
    _, _, full, _ = contact.stock.kernels()
    cut = set(map(tuple, edges))
    blank = full[:9, :9]
    pairs = contact.blocks()["full"]
    check = {}
    for load in (1, 2):
        # On each written seam (0,0), (load,0) and (0,load) have the
        # same head block; this equality includes every outgoing column.
        checked = 0
        for mag, (base, _) in pairs.items():
            heads = 9 * np.array(
                [
                    k + 1
                    for k, (_, m) in enumerate(contact.stock.LAW.LABELS)
                    if abs(m) == mag
                ]
            )
            for word in (0, load, 3 * load):
                allowed = heads + word
                expect = np.zeros_like(full[:, allowed])
                expect[allowed] = base
                np.testing.assert_array_equal(full[:, allowed], expect)
                checked += len(allowed)
        reduced = np.zeros((board.n, board.n), dtype=int)
        for a, b in board.seams:
            if (a, b) not in cut:
                reduced[np.ix_([a, b], [a, b])] += blank[
                    np.ix_([3 * load, load], [3 * load, load])
                ]
        # Six-slot completion: the contribution of missing blank neighbors
        # is local to each occupied site and independent of the head values.
        reduced[np.diag_indices(board.n)] += (6 - board.degree) * 4 * 8 * load**2
        vacuum_offset = 6 * 4 * 8 * load**2
        h = survey.lab.motion.cold.motion(board, cut, load, 8)[0]
        np.testing.assert_array_equal(
            reduced - vacuum_offset * np.eye(board.n), 4 * h.toarray()
        )
        check[str(load)] = dict(
            material_dimension=board.n,
            written_columns_per_contact=checked,
            relative_generator_sha256=sha256(
                (4 * h.toarray()).astype("<i8").tobytes()
            ).hexdigest(),
        )
    return check


def count_flags(rng, source, success):
    # Sample fresh joint heads and then the two physically separate detector
    # flags on each sharp head. Keep every failure; a head is gated only once.
    heads = rng.multinomial(SHOTS, source)
    bank = np.zeros(4, dtype=np.int64)
    for h, n in enumerate(heads):
        a, b = np.array([h // 2, h % 2]) * success
        bank += rng.multinomial(n, [(1 - a) * (1 - b), (1 - a) * b, a * (1 - b), a * b])
    return bank


def run():
    old = json.loads((HERE / "pair_contact.json").read_bytes())
    e = np.array(list(old["decoded"]["ordinary_energy_bands"].values()))
    # |cos(x)-cos(y)| <= |x-y| bounds a half-turn interrogation.
    t = dwell.DWELL
    center = (1 - np.cos(t * e.mean(axis=1)).mean()) / 2
    width = t * (e[:, 1] - e[:, 0]).mean() / 4
    calibration = survey.band(center, width)
    success = written_receiver()
    assert calibration[0] <= success <= calibration[1]
    np.testing.assert_allclose(success, (1 + np.cos(np.pi / 65)) / 4)
    delta = survey.lab.probe.log_factor(old["source"])
    strength = survey.lab.probe.mix_bands(delta, calibration**2)
    edges = old["physical_seams"]
    board = survey.lab.motion.cold.federation(0)
    _, live, origin = survey.chart(board, edges)
    paths = []
    for u, v in edges:
        paths.append(
            np.minimum(live[:, u, None] + live[v, :], live[:, v, None] + live[u, :])
        )
    # Equal clock triples give the same terminal observable. Acquire a fresh
    # flag bank per triple, with the address multiplicities preserved.
    start, stop = np.where(~np.eye(board.n, dtype=bool))
    triples = np.stack([a[start, stop] for a in (live, *paths)], axis=1)
    distinct, which, multiplicity = np.unique(
        triples, axis=0, return_inverse=True, return_counts=True
    )
    table = np.array(old["accepted_weight"], float)
    p = table.ravel() / table.sum()
    q = np.outer(table.sum(1), table.sum(0)).ravel() / table.sum() ** 2
    small = np.sqrt(np.log(8 * len(distinct) / 0.01) / SHOTS)
    rng = np.random.default_rng(636519)
    rows = []
    resolved = []
    gains = []
    for clock in distinct:
        direct, left, right = clock
        arrivals = np.array([left, right]) + 1 + t / np.pi
        slack = np.maximum(0, direct - arrivals)
        common = float(min(slack))
        times = np.array(
            [
                min(
                    [direct]
                    + [
                        a
                        for bit, a in zip((h // 2, h % 2), arrivals, strict=True)
                        if bit
                    ]
                )
                for h in range(4)
            ]
        )
        np.testing.assert_allclose(times @ [1, -1, -1, 1], common, atol=2e-14)
        banks = np.array([count_flags(rng, s, success) for s in (p, q)])
        assert np.all(banks.sum(1) == SHOTS)
        forecast = common * strength
        difference = (banks[1] - banks[0]) / SHOTS
        held = survey.band(float(difference @ times), float(np.ptp(times) * small))
        assert held[0] <= forecast[1] + 1e-13 and forecast[0] <= held[1] + 1e-13
        # Losing all assistance or using just one written route carries no
        # connected change at equal one-head marginals.
        for bit, saving in (([0, 0, 1, 1], slack[0]), ([0, 1, 0, 1], slack[1])):
            single = direct - saving * np.array(bit)
            assert abs(difference @ single) <= small * np.ptp(single) + 1e-12
        rows.append(
            dict(
                clock=clock.tolist(),
                banks=banks.tolist(),
                predicted=forecast.tolist(),
                measured=held.tolist(),
            )
        )
        resolved.append(held[0] > 0)
        gains.append(common)
    gains = np.array(gains)[which]
    ids = np.flatnonzero(gains == np.max(gains))
    i = int(ids[0])
    sel = start == origin
    pc = old["decoded"]["connected_charge_nats"]
    return dict(
        sha256={
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in ("PHASE_RECEIPTS.md", "pair_contact.json", "thermal_radar.json")
        },
        evacuated_free_probe=free_probe(board, edges),
        time_unit="one blank swap (pi)",
        accepted_per_arm_per_clock=SHOTS,
        native_dark_if_sender_and_positive_head=success,
        ordinary_calibration=calibration.tolist(),
        positive_charge_nats=pc,
        connected_joint_flag_prediction=strength.tolist(),
        joint_miss_threshold=float(strength[0] / 2),
        ideal_accepted_per_arm_at_one_percent=int(
            np.ceil(4 * np.log(100) / strength[0] ** 2)
        ),
        ordered_pairs=len(start),
        clock_triples=len(rows),
        pairs_sensitive=int((gains > 0).sum()),
        pairs_resolved=int(np.sum(np.asarray(resolved)[which])),
        survey_center=origin,
        center_destinations_sensitive=int((gains[sel] > 0).sum()),
        mean_predicted_delay_all=(gains.mean() * strength).tolist(),
        largest=dict(
            sender=int(start[i]),
            receiver=int(stop[i]),
            common_saving=float(gains[i]),
            **rows[which[i]],
        ),
        multiplicities=multiplicity.tolist(),
        strata=rows,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    out = run()
    path = HERE / "phase_receipts.json"
    if opt.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        contact.timing.lab.lab.cold.apparatus.check(json.loads(path.read_bytes()), out)
    print(
        {
            k: v
            for k, v in out.items()
            if k not in ("strata", "multiplicities", "sha256")
        }
    )
