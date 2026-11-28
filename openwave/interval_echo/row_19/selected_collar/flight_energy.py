"""Compare a source-only hopping-work measurement to a phase-referenced regional energy read."""

import argparse
from hashlib import sha256
import json
from math import comb

import numpy as np
from scipy.sparse.linalg import expm_multiply

from . import directed_service as flow

HERE, lab, iv, reg, region = flow.HERE, flow.lab, flow.iv, flow.reg, flow.region
STEP = 1 / 128
ORDER = 6
SEED = 28112123


def partition(values, ports):
    selected = values[:, ports]
    return np.column_stack([selected, values.sum(axis=1) - selected.sum(axis=1)])


def acquire(controls=None, *, seed=SEED):
    rng = np.random.default_rng(seed)
    native = lab.six.survey.native
    board = native.iso.pairing.formation.federation(1)
    route = native.iso.pairing.formation.region(board)
    ports = flow.dist.read("distributed_source")["ports"]
    if controls is None:
        controls = flow.dist.read("six_meter_pointer")["central_settings"]
    knobs = {int(t): complex(*v) for t, v in controls.items()}
    _, stocks, *_ = region.c.trials()
    history = {}
    norms = {}
    for sign in (-1, 1):
        Hs = [native.generator(board, route, sign, bool(b)).tocsr() / 4 for b in (0, 1)]
        norms[str(sign)] = max(float(abs(H).sum(axis=1).max()) for H in Hs)
        for home in ports:
            ramsey = {tau: {} for tau in knobs}
            bonds = {tau: [] for tau in knobs}
            for b, H in enumerate(Hs):
                gen = -1j * H
                init = np.zeros((H.shape[0], 4))
                init[4 * home : 4 * home + 4] = np.eye(4)
                tube = region.density.flight_columns(gen, init, 2 * max(knobs)).reshape(
                    2 * max(knobs) + 1, -1, 4, 4
                )
                neighbors = sorted(
                    {
                        int(j)
                        for j in H[4 * home : 4 * home + 4].nonzero()[1] // 4
                        if j != home
                    }
                )
                for tau, gamma in knobs.items():
                    mid = tube[tau]
                    for j in neighbors:
                        block = H[4 * j : 4 * j + 4, 4 * home : 4 * home + 4].toarray()
                        integer = np.rint(4 * block).astype(np.int64)
                        np.testing.assert_array_equal(4 * block, integer)
                        v = round(np.linalg.norm(integer, 2)) / 4
                        square = integer @ integer.T
                        np.testing.assert_array_equal(
                            square @ square, int((4 * v) ** 2) * square
                        )
                        X = stocks["reference"][b]
                        carrier = mid[home]
                        remote = mid[j]
                        co = np.trace(block @ carrier @ X @ remote.conj().T)
                        right = block @ carrier
                        left = block.conj().T @ remote
                        hit = (
                            np.trace(
                                right @ X @ right.conj().T + left @ X @ left.conj().T
                            ).real
                            / v**2
                        )
                        record = {"parity": b, "other": j, "unit": v}
                        for labname, mu in (
                            ("work", 2 * co.real),
                            ("incoming_current", 2 * co.imag),
                        ):
                            record[labname] = lab.counts(
                                np.array(
                                    [(hit + mu / v) / 2, (hit - mu / v) / 2, 1 - hit]
                                ),
                                rng,
                            )
                        bonds[tau].append(record)
                    click = mid @ mid[home]
                    miss = tube[2 * tau] - click
                    initial = np.concatenate(
                        [miss.reshape(-1, 4), click.reshape(-1, 4)], axis=1
                    )
                    evolved = expm_multiply(
                        gen,
                        initial,
                        start=0,
                        stop=ORDER * STEP,
                        num=ORDER + 1,
                        traceA=gen.diagonal().sum(),
                    ).reshape(ORDER + 1, -1, 4, 8)
                    m, c = np.split(evolved, 2, axis=-1)

                    def panels(amp, X):
                        weighted = amp @ X
                        norm = np.einsum(
                            "snij,snij->sn", weighted, amp.conj(), optimize=True
                        ).real
                        cross = np.einsum(
                            "snij,nij->sn", weighted, amp[0].conj(), optimize=True
                        )
                        return partition(norm, ports), partition(cross, ports)

                    for arm, g, X in (
                        ("dense", 0, stocks["reference"][b]),
                        ("ordinary", 0, stocks["ordinary"][b]),
                        ("overlap", gamma, stocks["reference"][b]),
                        ("bare", 1, stocks["reference"][b]),
                    ):
                        pop, co = panels(m + g * c, X)
                        if abs(g) < 1:
                            more, extra = panels(c, X)
                            pop += (1 - abs(g) ** 2) * more
                            co += (1 - abs(g) ** 2) * extra
                        total = pop + pop[0]
                        signed = -2 * co.imag
                        prob = np.stack([total + signed, total - signed], axis=-1) / 4
                        if arm not in ramsey[tau]:
                            ramsey[tau][arm] = np.zeros_like(prob)
                        ramsey[tau][arm] += prob
            for tau in knobs:
                history[f"{tau}/{home}/{sign}"] = {
                    "bonds": bonds[tau],
                    "ramsey": {
                        arm: [
                            np.reshape(lab.counts(row.ravel(), rng), (7, 2)).tolist()
                            for row in prob[1:]
                        ]
                        for arm, prob in ramsey[tau].items()
                    },
                }
    return {
        "starts": lab.SHOTS,
        "seed": seed,
        "step": STEP,
        "ports": ports,
        "knobs": controls,
        "generator_norms": norms,
        "histories": history,
    }


def signed(counts, n, unit):
    plus, minus = counts
    alpha = np.log(2 * reg.src.FAMILY / 0.01) / n
    hit = (plus + minus) / n
    width = unit * (
        np.sqrt(2 * alpha * (hit + reg.src.probability_error(hit, n)))
        + 4 * alpha / 3
        + 3e-13
    )
    mid = unit * ((plus - minus) / n)
    return reg.band(mid - width, mid + width)


def decode(data):
    iv.dps = 35
    dose = region.density.DOSE
    n = data["starts"]
    rows = data["histories"]
    source = {}
    for h, z in rows.items():
        k = complex(*data["knobs"][h.split("/")[0]])
        work = sum(signed(v["work"][:2], n, v["unit"]) for v in z["bonds"])
        current = sum(
            signed(v["incoming_current"][:2], n, v["unit"]) for v in z["bonds"]
        )
        source[h] = {
            "bare": dose * work,
            "overlap": dose * (k.real * work - k.imag * current),
        }

    def read(keys):
        total = n * len(keys)
        maxnorm = max(data["generator_norms"][k.rsplit("/", 1)[-1]] for k in keys)
        rem = STEP**ORDER * maxnorm ** (ORDER + 1) / (ORDER + 1)
        readings = {}
        for arm in ("dense", "ordinary", "overlap", "bare"):
            sites = {}
            for ix, name in enumerate(
                [*map(str, data["ports"]), "complement", "total"]
            ):
                val = iv.mpf(0)
                for j in range(1, ORDER + 1):
                    counts = [
                        sum(
                            tape[i][col]
                            for k in keys
                            for tape in [rows[k]["ramsey"][arm][j - 1]]
                            for i in (range(7) if name == "total" else [ix])
                        )
                        for col in range(2)
                    ]
                    coeff = (-1) ** (j + 1) * comb(ORDER, j)
                    val += coeff * signed(counts, total, 1) / (iv.mpf(j) * STEP)
                sites[name] = dose * (val + reg.band(-rem, rem))
            readings[arm] = sites
        forecasts = {
            arm: sum(source[k][arm] for k in keys) / len(keys)
            for arm in ("bare", "overlap")
        }
        changes = {
            arm: {
                name: lab.six.bounds(v - readings["dense"][name])
                for name, v in sites.items()
            }
            for arm, sites in readings.items()
            if arm != "dense"
        }
        return {
            "forecasts": {k: lab.six.bounds(v) for k, v in forecasts.items()},
            "late_minus_dense": changes,
            "dense_energy": {
                k: lab.six.bounds(v) for k, v in readings["dense"].items()
            },
        }

    result = {k: read([k]) for k in rows}
    groups = {t: read([k for k in rows if k.split("/")[0] == t]) for t in data["knobs"]}
    groups["mixture"] = read(list(rows))
    mix = groups["mixture"]
    cap = flow.dist.read("local_cap")
    caps = []
    for row in cap["balls"]:
        p = str(row["port"])
        a = iv.mpf(row["readings"]["overlap"]["area"])
        entropy = iv.mpf(row["readings"]["overlap"]["modular"])
        energy = iv.mpf(mix["late_minus_dense"]["overlap"][p])
        base = iv.mpf(mix["dense_energy"][p])
        fixed = entropy - iv.pi * energy / 2
        moved = fixed + iv.pi / 2 * base * a / iv.mpf(row["base_area"])
        caps.append(
            {
                "triangle": row["triangle"],
                "port": row["port"],
                "fixed_radius_remainder": lab.six.bounds(fixed),
                "fixed_volume_remainder": lab.six.bounds(moved),
            }
        )
    return {"histories": result, "groups": groups, "cap_remainders": caps}


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = acquire()
    r["result"] = decode(r)
    r["inputs"] = {
        n: sha256((HERE / n).read_bytes()).hexdigest()
        for n in (
            "FLIGHT_ENERGY_PROTOCOL.md",
            "six_meter_pointer.json",
            "local_cap.json",
            "distributed_source.json",
        )
    }
    path = HERE / "flight_energy.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    for k, z in r["result"]["groups"].items():
        print(
            k, z["forecasts"], {a: v["total"] for a, v in z["late_minus_dense"].items()}
        )
