"""Record protected four-home acknowledgements on independently evolved native registers."""

import argparse
import json
from hashlib import sha256

import numpy as np

from . import sector_register_pairs as pairs

HERE = pairs.collective.HERE


def read_words():
    words, _, ports, _, _, uv, loads, periods = pairs.chart.apparatus()
    homes = np.searchsorted(ports, [8, 20, 32, 44])
    internal = np.isin(uv, homes).sum(1) == 1
    edges = [np.flatnonzero((uv == h).any(1) & internal)[0] for h in homes]
    readings = np.array([periods[tuple(w > 0)][edges] for w in loads])
    classes, ids = np.unique(readings, axis=0, return_inverse=True)
    np.testing.assert_array_equal(readings - 1, words[:, :4] > 0)
    return ports[uv[edges]], classes, ids


def moments(counts):
    return {
        "registers": len(counts),
        "total_counts": counts.sum(0).tolist(),
        "sum_of_count_products": (counts.T @ counts).tolist(),
    }


def acquire():
    coll = pairs.collective
    links, classes, ids = read_words()
    _, _, subset, full, mats = pairs.chart.source.sources((8,))
    native = np.full(len(full), -1)
    native[subset] = np.arange(len(subset))
    tickets = coll.peer_register_lab.acquisition.returns.ticket_tables(mats)
    rng = np.random.default_rng(67839521)
    result = {
        "links": links.tolist(),
        "reported_periods": classes.tolist(),
        "protocol_sha256": sha256(
            (HERE / "SECTOR_RETURN_CLOCKS.md").read_bytes()
        ).hexdigest(),
    }

    def observe(members, ordinary):
        size = members.shape[1]
        result["size"] = size
        for name, starts in (
            ("peer", members),
            ("ordinary", ordinary[: len(ordinary) // size * size].reshape(-1, size)),
        ):
            batches = []
            batches_in_arm = 64 if name == "peer" else 65
            assert len(starts) % batches_in_arm == 0
            for part in np.array_split(starts, batches_in_arm):
                end, _ = coll.sector_sampling.report_from(
                    tickets, subset, native, rng, 8, part.ravel()
                )
                readings = ids[end].reshape(len(part), size)
                count = np.zeros((len(part), len(classes)), dtype=np.int64)
                np.add.at(count, (np.arange(len(part))[:, None], readings), 1)
                assert np.all(count.sum(1) == size)
                batches.append(moments(count))
            result[name] = batches

    coll.peer_register_lab.run(observe=observe, pulse_seed=907537, register_seed=814273)
    return result


def pair_law(nreg, total, products, size):
    single = total / (nreg * size)
    distinct = (products - np.diag(total)) / (nreg * size * (size - 1))
    independent = (np.outer(total, total) - products) / (nreg * (nreg - 1) * size**2)
    change = distinct - independent
    np.testing.assert_allclose(change.sum(0), 0, atol=1e-15)
    return change, single


def estimates(blocks, size):
    # Preserve independent whole registers within each equal-sized block.
    nr = np.array([x["registers"] for x in blocks])
    total = np.array([x["total_counts"] for x in blocks])
    prod = np.array([x["sum_of_count_products"] for x in blocks])
    n, t, p = nr.sum(), total.sum(0), prod.sum(0)
    change, single = pair_law(n, t, p, size)
    assert np.ptp(nr) == 0
    jk = np.array(
        [
            pair_law(n - x, t - y, p - z, size)[0]
            for x, y, z in zip(nr, total, prod, strict=True)
        ]
    )
    return change, single, jk


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    dst = HERE / "sector_return_clock_lab.json"
    data = acquire()
    if cli.parse_args().write:
        dst.write_text(json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n")
    else:
        assert data == json.loads(dst.read_text())
    print({k: sum(x["registers"] for x in data[k]) for k in ("peer", "ordinary")})
