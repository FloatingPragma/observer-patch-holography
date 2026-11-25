"""Mass-blind search of all directed two-loop head parities and charge detectors."""

import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag, expm

from ..two_winding_loops import measure as loop
from .certificate import certify_poles, commutant_dimension


def hermitian_blocks():
    # On one filled-background insertion the magnitude-one odd sector is
    # exactly 9/8 times the even sector after the common filled head subtraction.
    blocks = loop.blocks()
    even, odd = blocks[1, 1], blocks[1, -1]
    for h in [even, odd]:
        assert h.shape == (18, 18)

    def relative(h):
        return h - np.kron(h[4::9, 4::9], np.eye(9, dtype=int))

    np.testing.assert_array_equal(8 * relative(odd), 9 * relative(even))
    return loop.build(False)


def inventory():
    words, index, terms = hermitian_blocks()
    ins = [i for i, (site, heads) in enumerate(words) if site < 3]
    outs = sorted({(site, heads // 8) for site, heads in words if site >= 3})
    receivers = [[index[s, 8 * b + a] for s, b in outs] for a in (0, 7)]
    addresses = [w[0] for w in words]
    rows = []
    for bits in product(range(2), repeat=6):
        h = terms[-1] + sum(
            x + (x // 8) * b for x, b in zip(terms[:6], bits, strict=True)
        )
        multiplicity, identities = certify_poles(h, ins, receivers)
        u = expm(-2j * h)
        received = loop.choi(loop.kraus(u, index))
        rows.append(
            dict(
                parities="".join(map(str, bits)),
                positive_poles=multiplicity,
                zero_minor_identities=identities,
                full_port_commutant=commutant_dimension(h, addresses),
                fixed_port_acceptance=received.trace().real,
                broadband_pt_min=np.linalg.eigvalsh(loop.pt(received))[0],
            )
        )
    return rows, addresses


def native_meter():
    words, index, terms = loop.build(True)
    h8 = sum(terms)
    levels, v = np.linalg.eigh(h8 / 8)
    site = np.array([s for s, _ in words])
    positive = levels > 0
    transition = sum(abs(v.T @ ((site == i)[:, None] * v)) ** 2 for i in range(6))
    np.testing.assert_allclose(transition.sum(axis=0), 1, atol=2e-13)
    sub = transition[np.ix_(positive, positive)]
    costs = []
    for extra in product(range(3), repeat=len(sub) - 1):
        labels = (0,) + extra
        if set(labels) != {0, 1, 2} or labels.index(1) > labels.index(2):
            continue
        lab = np.array(labels)
        costs.append((sub[lab[:, None] != lab[None, :]].sum() / len(sub), labels))
    # Direct free evolution and two ordinary material measurements. Neither
    # energy projection nor a coherent external clock enters these banks.
    rng = np.random.default_rng(174901)
    u = expm(-2j * h8)
    after = expm(-0.5j * h8)
    banks = []
    for flag in [False, True]:
        # Different emitted head words on the same gateway. Both have net -1.
        initial = index[0, 7 * int(flag)]
        mid = u[:, initial]
        for meter in ["none", "owner", "port"]:
            groups = {"none": np.zeros_like(site), "owner": site // 3, "port": site}[
                meter
            ]
            dist = np.vstack(
                [
                    np.bincount(
                        site,
                        weights=abs(after @ (mid * (groups == a))) ** 2,
                        minlength=6,
                    )
                    for a in np.unique(groups)
                ]
            )
            np.testing.assert_allclose(dist.sum(), 1, atol=5e-13)
            counts = rng.multinomial(2**20, dist.ravel() / dist.sum()).reshape(
                dist.shape
            )
            banks.append(
                dict(
                    source_word=7 * int(flag),
                    first_meter=meter,
                    joint_counts=counts.tolist(),
                    last_owner=int(counts[:, 3:].sum()),
                    last_gateway=int(counts[:, 3].sum()),
                )
            )
    dims = {
        str(odd): {
            key: commutant_dimension(sum(loop.build(odd)[2]), group)
            for key, group in [("owner", site // 3), ("port", site)]
        }
        for odd in [False, True]
    }
    return dict(
        positive_transitions=sub.tolist(),
        positive_to_negative=float(
            transition[np.ix_(~positive, positive)].sum() / len(sub)
        ),
        partitions=len(costs),
        minimum_positive_crossing=min(costs)[0],
        minimizing_labels=min(costs)[1],
        dimensions=dims,
        banks=banks,
        starts_per_bank=2**20,
        seed=174901,
    )


def known_controls():
    # A known spin-1/2 particle on a connected three-site orbital path.
    # Site-resolved scalar charge counters act identically on its spin.
    orbital = np.array([[1, -1, 0], [-1, 3, -1], [0, -1, 2]])
    identity = np.eye(2, dtype=int)
    x = np.array([[0, 1], [1, 0]])
    y = np.array([[0, -1j], [1j, 0]])
    z = np.diag([1, -1])
    sites = np.repeat(np.arange(3), 2)
    h = np.kron(orbital, identity)
    spin_rot = [-1j * x, -1j * y]
    a, b = spin_rot
    np.testing.assert_array_equal(a @ b, -b @ a)
    np.testing.assert_array_equal(a @ a, -identity)
    split = h + np.kron(np.eye(3, dtype=int), 4 * z)
    flavors = block_diag(
        h, h + 20 * np.eye(6, dtype=int), h + 50 * np.eye(6, dtype=int)
    )
    # Explicit exact species projectors for three imposed distinct rest offsets.
    family = np.repeat(np.arange(3), 6)
    for f in range(3):
        tag = np.diag((family == f).astype(int))
        assert not (tag @ flavors - flavors @ tag).any()

    def energy_blind_result(hamiltonian):
        values, columns = np.linalg.eigh(hamiltonian)
        components = []
        for energy in np.unique(np.round(values, 10)):
            band = columns[:, abs(values - energy) < 1e-9]
            projector = band @ band.T
            # Same left emitter and right detector on every spin component.
            components.append(projector[4:6, :2])
        choi = loop.choi(components)
        return dict(
            mean_arrival=float(choi.trace().real),
            pt_min=float(np.linalg.eigvalsh(loop.pt(choi))[0]),
        )

    return dict(
        unsplit_energy_read=energy_blind_result(h),
        split_energy_read=energy_blind_result(split),
        scalar_commutant=commutant_dimension(orbital, list(range(3))),
        unsplit_spin_commutant=commutant_dimension(h, sites),
        split_spin_commutant=commutant_dimension(split, sites),
        split_energy_unique=len(np.unique(np.round(np.linalg.eigvalsh(split), 10)))
        == 6,
        spin_order_sign=-1,
        scalar_order_sign=1,
        three_flavor_commutant=commutant_dimension(flavors, np.tile(sites, 3)),
    )


def run():
    rows, _ = inventory()
    totals = Counter()
    for r in rows:
        totals.update(r["positive_poles"])
    return dict(
        preparations=rows,
        positive_roots_by_multiplicity=dict(totals),
        broadband_negative_at_16=sum(
            bool(r["broadband_pt_min"] < -1e-10) for r in rows
        ),
        benches_with_repeated_positive=sum(
            2 in r["positive_poles"] and r["positive_poles"][2] > 0 for r in rows
        ),
        native_readout=native_meter(),
        comparison=known_controls(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = json.dumps(run(), indent=2) + "\n"
    if parser.parse_args().write:
        Path(__file__).with_name("results.json").write_text(data)
    else:
        print(data, end="")
