"""Lower bounds on missing source dimensions from independently checkable interval minors."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import qr
import sympy as sp

from .six_meter_pointer import HERE, bounds, components, positive


def pivots(matrix, count):
    return qr(matrix, mode="economic", pivoting=True)[2][:count]


def certificate(rows, columns, names, histories, bands):
    mid, rad = components(
        [[bands[histories[c]][names[r]] for c in columns] for r in rows]
    )
    try:
        proposed = np.linalg.inv(np.array(mid, dtype=float))
    except np.linalg.LinAlgError:
        return None
    k = len(rows)
    # Store each binary64 inverse entry exactly using its round-trip JSON representation.
    K = sp.Matrix(proposed).applyfunc(sp.Rational)
    H = positive(sp.eye(k) - K * mid) + positive(K) * rad
    try:
        scale = np.linalg.solve(np.eye(k) - np.array(H, dtype=float), np.ones(k))
    except np.linalg.LinAlgError:
        return None
    if not np.all(np.isfinite(scale)) or min(scale) <= 0:
        return None
    v = sp.Matrix(scale).applyfunc(sp.Rational)
    alpha = max(z / w for z, w in zip(H * v, v, strict=True))
    if alpha >= 1:
        return None
    return {
        "rank_lower": k,
        "meters": [names[r] for r in rows],
        "histories": [histories[c] for c in columns],
        "preconditioner": proposed.tolist(),
        "norm_weights": scale.tolist(),
        "contraction_upper": bounds(alpha, 0)[1],
    }


def rank(names, ordinary):
    bands = {h: z["observed_minus_prediction"] for h, z in ordinary.items()}
    histories = list(bands)
    M = np.array([[np.mean(bands[h][n]) for h in histories] for n in names])
    M = M / np.maximum(np.max(abs(M), axis=1, keepdims=True), np.finfo(float).tiny)
    M = M / np.maximum(np.max(abs(M), axis=0, keepdims=True), np.finfo(float).tiny)
    U, _, V = np.linalg.svd(M, full_matrices=False)
    for k in range(min(M.shape), 0, -1):
        c = pivots(M, k)
        r = pivots(M.T, k)
        for rows, cols in [
            (pivots(U[:, :k].T, k), pivots(V[:k], k)),
            (pivots(M[:, c].T, k), c),
            (r, pivots(M[r], k)),
        ]:
            # Quick floating point rejection saves rational work on bad proposals.
            sub = np.array(
                [[bands[histories[j]][names[i]] for j in cols] for i in rows]
            )
            a = sub.mean(axis=2)
            try:
                h = abs(np.linalg.inv(a)) @ ((sub[:, :, 1] - sub[:, :, 0]) / 2)
                if max(abs(np.linalg.eigvals(h))) >= 1:
                    continue
            except np.linalg.LinAlgError:
                continue
            proof = certificate(rows, cols, names, histories, bands)
            if proof is not None:
                return proof
    return {"rank_lower": 0}


def run():
    p = HERE / "regional_source_calibration.json"
    data = json.loads(p.read_bytes())
    labels = data["labels"]
    groups = {
        "shape": [n for n in labels if n.startswith("ratio/")],
        "material": [n for n in labels if not n.startswith("ratio/")],
        "both": labels,
    }
    return {
        "results": {
            k: {g: rank(ns, z["ordinary"]) for g, ns in groups.items()}
            for k, z in data["calibrations"].items()
            if z["status"] == "verified inverse"
        },
        "inputs": {
            f.name: sha256(f.read_bytes()).hexdigest()
            for f in (p, HERE / "REGIONAL_RANK.md")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    path = HERE / "regional_rank.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    print(
        {k: {g: v["rank_lower"] for g, v in z.items()} for k, z in r["results"].items()}
    )
