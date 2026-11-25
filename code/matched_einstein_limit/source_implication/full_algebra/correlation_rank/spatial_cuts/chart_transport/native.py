"""Directed native four-tick source effects of the three partial translations."""

import argparse
import json

import numpy as np
import sympy as sp
from scipy.linalg import qr

from ..radial_profiles import effects as fixed
from . import mounts

HERE = mounts.HERE


def choose():
    """Select the 31 response coordinates on each mount using material data only."""
    selected = []
    for records in mounts.float_effects():
        keys, values = list(records), np.array(list(records.values()))
        norms = np.linalg.norm(values, axis=1)
        assert min(norms) > 0
        _, _, ids = qr(
            (values / norms[:, None] ** 0.3).T, mode="economic", pivoting=True
        )
        selected.append([keys[i] for i in sorted(ids[:31])])
    return selected


def exact_rows(selected, apparatus=mounts.apparatus):
    reg = fixed.cuts.regional
    board, route, launch, patches = apparatus()
    edge_sets = [fixed.cuts.native.flux_on_edges(board, p) for p in patches]
    coefficient = fixed.coefficients()
    joint = [reg.cert.prior.sym_matrix(b) for b in reg.cert.iso.basis(True)]
    local = mounts.full_source.local_basis()
    effects = [{} for _ in patches]
    for q in (-1, 1):
        for odd in (False, True):
            h = reg.bridge.native.generator(board, route, q, odd)
            enrolled = []
            for patch, edges, keys in zip(patches, edge_sets, selected):
                mat, ops = mounts.channels(board, route, patch, edges, q, odd)
                assert (mat - h).count_nonzero() == 0
                enrolled.append(
                    [((sign, name), *ops[name]) for sign, name in keys if q == sign]
                )
            idx = np.unique(
                np.concatenate(
                    [op.nonzero()[0] for arm in enrolled for _, op, _ in arm]
                )
            )
            w = fixed.columns(h, idx, launch, coefficient)
            for arm, rows in zip(enrolled, effects):
                for k, op, unit in arm:
                    ids = np.unique(op.nonzero()[0])
                    sites = len(set(ids // 4))
                    norm = int(np.ceil(abs(op).sum(axis=1).max() / abs(unit)))
                    assert sites <= 960 and norm <= 200
                    compressed = reg.region.compress(
                        w[:, np.searchsorted(idx, ids)], op[ids, :][:, ids], unit
                    )
                    compressed *= sp.Rational(
                        reg.bridge.native.prior.FRAC**2, fixed.UNIT**2
                    )
                    j = [reg.cert.inner(compressed, b) for b in joint]
                    m = [reg.cert.inner(compressed, b) for b in local]
                    err = sp.Rational(5 * sites * max(1, norm), 10**56)
                    rows.setdefault(k, []).append((j, m, compressed.trace(), err))
    data = []
    for entries, keys, patch in zip(effects, selected, patches):
        rows = []
        for sign, name in keys:
            (j, m, tr, e), (jj, mm, ttr, ee) = entries[sign, name]
            centers = j + jj + m + mm + [tr - ttr]
            errors = [e] * 9 + [ee] * 9 + [e] * 6 + [ee] * 6 + [2 * (e + ee)]
            row = [
                [
                    int(sp.floor((c - r) * fixed.SCALE)),
                    int(sp.ceiling((c + r) * fixed.SCALE)),
                ]
                for c, r in zip(centers, errors)
            ]
            rows.append(dict(charge=sign, reading=name, row=row))
        data.append(
            dict(
                center=patch["center"],
                target=patch["target"],
                orientation=patch["orientation"],
                cells=patch["cells"],
                images=patch["images"],
                rows=rows,
            )
        )
    return data


def store(data, path=HERE / "effects.json"):
    # One physical detector per line, with even/odd Pauli order specified by mounts.source_components.
    arms = []
    for item in data:
        header = json.dumps({k: v for k, v in item.items() if k != "rows"})[:-1]
        rows = ",\n".join(json.dumps(r, separators=(",", ":")) for r in item["rows"])
        arms.append(header + ', "rows": [\n' + rows + "\n]}")
    path.write_text(
        '{"scale": '
        + str(fixed.SCALE)
        + ', "charts": [\n'
        + ",\n".join(arms)
        + "\n]}\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    record = HERE / "effects.json"
    if record.exists():
        previous = json.loads(record.read_text())
        keys = [
            [(v["charge"], v["reading"]) for v in c["rows"]] for c in previous["charts"]
        ]
    elif args.write:
        keys = choose()
    else:
        parser.error("the native effect archive is missing; use --write to generate it")
    result = exact_rows(keys)
    if args.write:
        store(result)
    else:
        assert result == previous["charts"] and previous["scale"] == fixed.SCALE
    print("Verified the native finite effects of all three material charts.")
