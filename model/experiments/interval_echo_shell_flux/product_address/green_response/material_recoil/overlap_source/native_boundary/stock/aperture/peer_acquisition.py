"""Release native first-return words after two independent local reports agree."""

import argparse
from hashlib import sha256
import json
import math

import numpy as np

from . import collective_source as source
from . import landmark_returns as returns

HERE = returns.HERE


def acquire(rng, size, subset, live, tables):
    # Each proposal has its own initial native word and independent pulse tickets.
    # Retain trial identities in the order in which their first returns arrive.
    out = np.full(size, -1)
    sites = rng.choice(subset, size)
    pending = np.arange(size)
    for step in range(1, 4097):
        for table in tables:
            sites = table[sites, rng.integers(128, size=len(sites))]
        done = live[sites]
        out[pending[done]] = sites[done]
        sites, pending = sites[~done], pending[~done]
        if not len(pending):
            return out, step
    raise RuntimeError(
        f"{len(pending)} unresolved native first returns after {step} sweeps"
    )


def run():
    full, initial, _, service, _ = returns.retention.setup()
    live = initial > 0
    subset = np.flatnonzero(live)
    _, w, c, _, codes, h, *_ = source.tidal.lab.energy.pair.setup()
    index = {tuple(word): i for i, word in enumerate(w)}
    native = np.full(len(full), -1)
    native[subset] = [index[tuple(word)] for word in full[subset]]
    pi = source.records()[1]["original"]
    joint = np.array([pi @ (c == j) for j in range(4)], object).reshape(2, 2)
    paired = np.outer(joint.sum(1), joint.sum(0)).ravel()
    shots = 2**18
    targets = paired[c]
    agree = pi @ targets
    once = pi * targets / agree
    decline = pi * (1 - targets) / (1 - agree)
    renewed = pi * targets / joint.ravel()[c]
    assert np.all(once * agree + decline * (1 - agree) == pi)
    assert np.array_equal(renewed, source.records()[1]["source"])
    expected = sum(paired / joint.ravel())
    tables = returns.ticket_tables([m for m, _ in service])
    rng = np.random.default_rng(857241)
    longest = 0

    def bank(size):
        nonlocal longest
        words, duration = acquire(rng, size, subset, live, tables)
        longest = max(longest, duration)
        return native[words]

    left, right = bank(shots), bank(shots)
    requested = 2 * (c[left] // 2) + c[right] % 2
    first = bank(shots)
    yes = c[first] == requested
    members = first.copy()
    pending = np.flatnonzero(~yes)
    attempts = shots
    rounds = 1
    while len(pending):
        if rounds >= 200:
            raise RuntimeError(f"{len(pending)} unmatched reporter pairs")
        batch = bank(len(pending))
        good = c[batch] == requested[pending]
        members[pending[good]] = batch[good]
        attempts += len(pending)
        pending = pending[~good]
        rounds += 1
    assert np.array_equal(c[members], requested)
    scores = {
        "blank": w[:, -1] == 0,
        "selective": np.isin(w[:, -1], source.PARTIAL),
        "signed_square": (w[:, 1] - 1) ** 2,
        **{f"head_cell_{j}": c == j for j in range(4)},
    }
    energy = np.vstack(
        (w[:, :4].T, h.diagonal()[codes], source.tidal.lab.energy.returned(w)[-1])
    )
    result = {}
    # Each single-attempt arm is iid conditional on its observed size.
    # The union bound also covers correlated readings in different arms.
    for arm, proposal, law in (
        ("ledger_only", first, pi),
        ("first_matches", first[yes], once),
        ("first_declines", first[~yes], decline),
        ("renewed_material", members, renewed),
    ):
        eps = math.sqrt(
            math.log(2 * (4 * len(scores) + 4) / 0.01) / (2 * len(proposal))
        )
        row = {}
        for name, flag in scores.items():
            theoretical = law @ flag
            observed = int(flag[proposal].sum())
            assert abs(observed / len(proposal) - float(theoretical)) <= eps
            row[name] = dict(count=observed, exact=str(theoretical))
        result[arm] = dict(
            trials=len(proposal),
            simultaneous_error=eps,
            readings=row,
            ordinary_home_and_energy=list(map(str, energy @ law)),
        )
    counts = np.bincount(requested, minlength=4)
    eps = math.sqrt(math.log(2 * (4 * len(scores) + 4) / 0.01) / (2 * shots))
    assert max(abs(counts / shots - np.array(paired, float))) <= eps
    return dict(
        sha256={
            p: sha256((HERE / p).read_bytes()).hexdigest()
            for p in ("PEER_RETURN.md", "landmark_returns.json")
        },
        paired_report_law=list(map(str, paired)),
        paired_report_counts=counts.tolist(),
        paired_report_error=eps,
        first_agreement=str(agree),
        expected_material_attempts_per_pair=str(expected),
        total_material_attempts=attempts,
        declined_material_attempts=attempts - shots,
        total_reporter_returns=2 * shots,
        maximum_renewals_including_first=rounds,
        maximum_sweeps_for_one_return=longest,
        unresolved=0,
        joint_coverage=0.99,
        arms=result,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    raw = json.dumps(result, indent=2, allow_nan=False) + "\n"
    path = HERE / "peer_acquisition.json"
    if args.write:
        path.write_text(raw)
    else:
        assert raw == path.read_text()
    print(
        result["total_material_attempts"],
        result["unresolved"],
        result["first_agreement"],
    )
