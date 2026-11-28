"""Assemble full accepted ordinary-matched registers from fresh native returns."""

import argparse
from bisect import bisect_right
from fractions import Fraction as Q
from hashlib import sha256
from itertools import accumulate
import json
import math
import random

import numpy as np

from . import peer_acquisition as acquisition
from . import peer_population as population

fibers = population.fibers
HERE = population.HERE


def run(*, observe=None, pulse_seed=177813, register_seed=440781):
    k, shots = 16, 4096
    w, c, groups, cat, types, _, _, _, _ = fibers.apparatus()
    bank = []
    theory, ensembles = population.histogram(k, draw_bank=bank)
    size = theory["size"]
    paired = fibers.BLANK_PAIRS
    protected = np.flatnonzero(cat < 0)
    native_cat = cat.copy()
    native_cat[protected] = len(types) + np.arange(len(protected))
    word_quota, remainder = divmod(size, len(w))
    assert remainder == 0
    full, initial, _, service, _ = acquisition.returns.retention.setup()
    live = initial > 0
    subset = np.flatnonzero(live)
    lookup = {tuple(word): i for i, word in enumerate(w)}
    enrolled = np.full(len(full), -1)
    enrolled[subset] = [lookup[tuple(word)] for word in full[subset]]
    tables = acquisition.returns.ticket_tables([m for m, _ in service])
    rng = np.random.default_rng(pulse_seed)
    tickets = random.Random(register_seed)
    cumulative = list(accumulate(number for number, _ in bank))
    blank = np.array([i[-1] for i in types])
    first_cell = np.array([i[0] == 0 for i in types])
    flags = [w[:, -1] == 0, c == 0]
    counts = np.empty((shots, len(types) + len(protected)), int)
    # Choose a full register in proportion to the exact number of native
    # proposal sequences and peer permutations it represents. No measured
    # blank frequency or receiver outcome enters the choice.
    for row in counts:
        _, v = bank[bisect_right(cumulative, tickets.randrange(cumulative[-1]))]
        free = v.copy()
        totals = np.array([v[i] + v[j] for i, j in paired])
        num_blank = int(sum(v[j] for _, j in paired))
        parts = rng.multivariate_hypergeometric(totals, num_blank)
        for (i, j), n, q in zip(paired, totals, parts, strict=True):
            free[i], free[j] = n - q, q
        row[: len(types)] = free
        row[len(types) :] = word_quota
    assert np.all(counts.sum(1) == size)
    # Materialize every requested member on actual independent first
    # returns. Fill slots of each class in arrival order. Surplus returns
    # never become members. Uniform ordering within the batch is immaterial
    # for its counts and for the independently random released slot.
    need = counts.sum(0)
    stored = [np.empty(x, np.int16) for x in need]
    taken = np.zeros_like(need)
    attempts = longest = 0
    while np.any(taken < need):
        fresh, steps = acquisition.acquire(rng, 65536, subset, live, tables)
        longest = max(longest, steps)
        words = enrolled[fresh]
        labels = native_cat[words]
        for j in range(len(need)):
            keep = words[labels == j][: need[j] - taken[j]]
            stored[j][taken[j] : taken[j] + len(keep)] = keep
            taken[j] += len(keep)
        attempts += len(fresh)
    members = np.empty((shots, size), np.int16)
    offsets = np.zeros(len(need), int)
    for i, row in enumerate(counts):
        members[i] = np.concatenate(
            [stored[j][offsets[j] : offsets[j] + n] for j, n in enumerate(row)]
        )
        offsets += row
    assert np.array_equal(offsets, need)
    # A conditioned peer permutes right labels within the fixed left blocks
    # after allocating the specified number of zeros to each block.
    rights_shuffled = 0
    for row in members:
        left = c[row] // 2
        right = c[row] % 2
        peer = right.copy()
        for x in (0, 1):
            peer[left == x] = rng.permutation(right[left == x])
        assert np.array_equal(
            np.bincount(2 * left + peer, minlength=4), np.bincount(c[row], minlength=4)
        )
        rights_shuffled += int(np.sum(peer != right))
    reference, refsteps = acquisition.acquire(rng, 2**20, subset, live, tables)
    reference = enrolled[reference]
    if observe is not None:
        observe(members, reference)
    released = members[np.arange(shots), rng.integers(size, size=shots)]
    original = fibers.source.records()[1]["original"]
    z = ensembles["agreed"]
    # Bernstein on independent full batches. All members of one batch
    # are correlated. The two conserved statistics have the same sum in
    # every triple split represented by one aggregated type.
    log = math.log(2 * 6 / 0.006)
    readings = {}
    for name, flag, kind in zip(
        ("blank", "signed_square"), flags, (blank, first_cell), strict=True
    ):
        assert np.array_equal(flags[1], (w[:, 1] - 1) ** 2)
        avg = z @ flag
        var = (
            sum(n * (Q(int(v @ kind), size) - avg) ** 2 for n, v in bank)
            / cumulative[-1]
        )
        width = math.sqrt(2 * float(var) * log / shots) + 2 * log / (3 * shots)
        x = int(flag[members].sum())
        base = int(flag[reference].sum())
        err = math.sqrt(log / (2 * len(reference)))
        obs = x / (shots * size)
        assert abs(obs - float(avg)) <= width
        assert abs(base / len(reference) - float(original @ flag)) <= err
        delta = obs - base / len(reference)
        interval = [delta - width - err, delta + width + err]
        assert interval[0] > 0
        nrel = int(flag[released].sum())
        assert abs(nrel / shots - float(avg)) <= math.sqrt(log / (2 * shots))
        readings[name] = dict(
            register_counts=x,
            reference_counts=base,
            released_counts=nrel,
            exact_change=str(avg - original @ flag),
            observed_change=delta,
            change_bounds=interval,
            exact_batch_variance=str(var),
        )
    _, _, _, _, codes, h, *_ = fibers.source.tidal.lab.energy.pair.setup()
    energy = np.vstack(
        (
            w[:, :4].T,
            h.diagonal()[codes],
            fibers.source.tidal.lab.energy.returned(w)[-1],
        )
    )
    assert np.array_equal(energy @ z, energy @ original)
    # Check all eight component-calibrated ordinary totals in each register.
    # The empirical sharp-word values fluctuate within their exchange group.
    lev = np.empty((len(native_cat), 8), object)
    for g in np.unique(groups):
        mask = groups == g
        lev[mask] = [Q(int(x), int(mask.sum())) for x in energy[:, mask].sum(1)]
    cat_lev = [lev[np.flatnonzero(native_cat == i)[0]] for i in range(counts.shape[1])]
    assert np.array_equal(np.array(cat_lev)[native_cat], lev)
    for j in range(8):
        vals = [row[j] for row in cat_lev]
        unit = math.lcm(*(q.denominator for q in vals))
        exact = np.array([int(v * unit) for v in vals])
        assert np.all(counts @ exact == int((energy @ original)[j] * size * unit))
    selective = np.isin(w[:, -1], fibers.source.PARTIAL)
    assert np.all(z[selective] == original[selective])
    x = json.loads((HERE / "landmark_arrival.json").read_bytes())
    panel = x["panels"]["25/36"]
    exact_delta = Q(panel["matched_partial_fixed_stations"]["104"]) / Q(
        x["ordinary_energy_matched"]["preparations"]["entropy_stationary"][
            "blank_change"
        ]
    )
    transfer = panel["held_transfer"]["open_minus_closed"]
    prediction = [u * v for u in readings["blank"]["change_bounds"] for v in transfer]
    ds = (z - original) @ ((w[:, :4] == 1).sum(1))
    coeff = (4 * (z - original) @ flags[0] / (original @ selective), ds, -ds)
    return dict(
        sha256={
            p: sha256((HERE / p).read_bytes()).hexdigest()
            for p in ("PEER_RETURN.md", "peer_population.json", "landmark_arrival.json")
        },
        size=size,
        batches=shots,
        reference_returns=len(reference),
        members=members.size,
        preparation_return_attempts=attempts,
        surplus=attempts - members.size,
        longest_first_return=max(longest, refsteps),
        unresolved=0,
        changed_right_positions=rights_shuffled,
        copies_per_selective_word_in_every_register=word_quota,
        readings=readings,
        exact_ordinary_home_and_energy=list(map(str, energy @ z)),
        empirical_ordinary_home_and_energy=energy[:, members].mean((1, 2)).tolist(),
        joint_local_coverage=0.994,
        exact_transfer=str(exact_delta),
        exact_arrival_change=str(exact_delta * Q(readings["blank"]["exact_change"])),
        shared_arrival_change_bounds=[min(prediction), max(prediction)],
        # Held transfer uses a separate 0.003 error budget.
        joint_source_and_transfer_coverage_at_least=0.991,
        ordinary_coefficients=list(map(str, coeff)),
        passive_probability=theory["acceptance_given_free_register"],
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    raw = json.dumps(out, indent=2, allow_nan=False) + "\n"
    path = HERE / "peer_register_lab.json"
    if args.write:
        path.write_text(raw)
    else:
        assert raw == path.read_text()
    print(
        out["members"],
        out["readings"]["blank"]["change_bounds"],
        out["shared_arrival_change_bounds"],
    )
